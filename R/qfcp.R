#' @importFrom MASS ginv
#' @importFrom Matrix bdiag
#' @importFrom stats optim pnorm dnorm
#' @export

gaussian<-function(x,kernel_sigma){
  exp(-x^2/(2*kernel_sigma^2))
}

getgram<-function(s,kernel_sigma){
    kernel_sigma = 2*kernel_sigma^2
    Gram = outer(s,s,FUN=function(x,y) exp(-(x-y)^2/kernel_sigma))
    return(Gram)
}

smooth_function <- function(x, d = 0, type = "norm_cdf"){
    if(d==0){
        if(type == "indicator"){
            g = as.numeric(x>0)
        }
        else if (type == "norm_cdf"){
            g = pnorm(x)
        }
        else if (type == "sigmoid") {
            g = 1/(1+exp(-x)) 
        } 
        else { 
        g = pnorm(x) + x * dnorm(x)
        }
    }
    else if(d==1){
         if(type == "indicator"){
            stop("Indicator function is un-smoothing function!")
        }
        else if (type == "norm_cdf"){
            g = dnorm(x)
        }
        else if (type == "sigmoid") {
            g = 1/(1+exp(-x))^2 
        } 
        else { 
        g = (1 - x - x^2) * dnorm(x)
        }        
    }
    else if(d == 2){
        if(type == "indicator"){
            stop("Indicator function is un-smoothing function!")
        }
        else if(type == "norm_cdf"){
            g = -x * dnorm(x)
        }

    }

    return(g)
}

qfcp<- function(x, y, z, sub_x, m, kernel_sigma, sigma, lambda, tau) {
    n         = dim(x)[1]
    p         = dim(x)[2]
    d         = dim(sub_x)[2]
    q         = dim(z)[2]
    o         = -2/3
    s         = seq(0, 1, length.out = m )   
    kern      = getgram(s, kernel_sigma)  
    f         = matrix(rep(0, n*m), n, m)   
    hat_u     = matrix(rep(0, n * m), n, m)                                
    hat_w     = matrix(rep(0,n*m), n, m)                       
    hat_gamma = rep(0, dim(z)[2]-1)
    omega     = kronecker(as.matrix(bdiag(diag(1,p,p), diag(1,d,d))), kern) 
    pk        = kronecker(rep(1,(p+d)),kern)
    ch        = 1
    iters     = 30

    iter      = 0
    while(iter<iters){
        
        iter  = iter + 1
        h = ch*sd(z[,1]+z[,-1] %*% hat_gamma)*n^(o)
        # h  = log(n)*n^(-1)
        g  = z[,1]+z[,-1] %*% hat_gamma 
        G  = as.vector(smooth_function(g/h))
        xx = cbind(x, (sub_x * G))
        ginv_xx   =  ginv(m * (t(xx)  %*% xx))

        temp      = ((-f + y + hat_w)>(tau/sigma) + 0) * tau/sigma + ((-f + y + hat_w)<(tau-1)/sigma + 0) * (tau-1)/sigma
        hat_u     = ( -f + y + hat_w) * (((-f + y + hat_w)>(tau/sigma) | (-f + y + hat_w)<(tau-1)/sigma ) + 0) - temp 
        hat_u     = y - hat_u 
        rr_1      = hat_u  + hat_w
        xxx = 0
        kkk = 0
        for(i in 1:n){
            xxx = xxx + (xx[i,]) %*% t(xx[i,])
            }
        for(j in 1:m){
            kkk = kkk + (kern[,j]) %*% t(kern[,j])
            }
        xxkk  = kronecker(xxx, kkk)
        xnp        = kronecker(t(xx),matrix(1,m,m)) * kronecker(matrix(1,(p+d),n),kern)
        ginv_xxkk  = ginv(xxkk + 2*lambda/sigma * omega)
        uxw    = hat_u  + hat_w
        xu     = xnp %*% as.vector(t(uxw))
        hat_d  = ginv_xxkk %*% as.vector(xu)
        fn <- function(gamma){
            # h  = ch*sd(z[,1]+z[,-1] %*% gamma)*n^(o)
            h = log(n)*n^(-1)
            g  = z[,1]+z[,-1] %*% gamma 
            G  = as.vector(smooth_function(g/h))
            xx = cbind(x, sub_x * G) 
            f1 = kronecker(xx, t(rep(1,m))) %*% (pk * t(kronecker(t(hat_d), rep(1,m))))	
            rss = sum((hat_u - f1 + hat_w)^2)*sigma/2 
            return(rss) 
        } 

        R = optim(par=hat_gamma, fn=fn,gr=NULL)     
        hat_gamma = R$par 
        f    = kronecker(xx, t(rep(1,m))) %*% (pk * t(kronecker(t(hat_d), rep(1,m))))	
        hat_w  = hat_w + hat_u - f 
    }

    return (list(gamma=hat_gamma, d=hat_d , w=hat_w))
}

qfcp0 <- function(x, y, m, kernel_sigma, sigma, lambda, tau) {
    
    n         = dim(x)[1]
    p         = dim(x)[2]
    s         = seq(0, 1, length.out = m )   
    kern      = getgram(s, kernel_sigma)      # m x m
    
    hat_u     = matrix(rep(0, n * m), n, m)   # (p+d)m x 1       
    hat_w     = matrix(rep(0,n*m), n, m)      # n x m
    iter      = 0
    hat_b     = rep(0,m*p)
    f         = matrix(0,n,m)
    omega     = kronecker(diag(1,p,p),kern)

    xxkk      = matrix(rep(0,m*p*m*p), m*p, m*p)
    ginv_xx   = ginv(m * (t(x)  %*% x)) 
	xm		  = kronecker(x, t(rep(1,m)))
	pk        = kronecker(rep(1,p),kern)
    iters     = 30

    xxx = 0
    kkk = 0
    for(i in 1:n){
        xxx = xxx + (x[i,]) %*% t(x[i,])
        }
    for(j in 1:m){
        kkk = kkk + (kern[,j]) %*% t(kern[,j])
        }
    xxkk  = kronecker(xxx, kkk)

    ginv_xxkk  = ginv(xxkk + 2*lambda/sigma * omega)
    xnp        = kronecker(t(x),matrix(1,m,m)) * kronecker(matrix(1,p,n),kern)

    while(iter<iters){
 
        iter   = iter + 1
        temp   = ((-f + y + hat_w)>(tau/sigma) + 0) * tau/sigma + ((-f + y + hat_w)<(tau-1)/sigma + 0) * (tau-1)/sigma
        hat_u  = (-f + y + hat_w) * (((-f + y + hat_w)>(tau/sigma) | (-f + y + hat_w)<(tau-1)/sigma ) + 0) - temp 
        hat_u  = y - hat_u 
        uxw    = hat_u  + hat_w
        xu     = xnp %*% as.vector(t(uxw))
        hat_b  = ginv_xxkk %*% as.vector(xu)
		f    = xm %*% (pk * t(kronecker(t(hat_b), rep(1,m))))	
        hat_w  = hat_w + hat_u - f 
        df	   = sum((hat_w>(tau-1))*(hat_w<tau))
        res     = y - f
    }
    return (list(d=hat_b,res=res,haty=f,Df = df))

}

