#' @importFrom MASS ginv
#' @importFrom stats cor
#' @export

ADMM_quantile_CP <- function(x, y, m, kernel_sigma, sigma, iters, tau, ginv_xxkk, kern) {
    n         = dim(x)[1]
    p         = dim(x)[2]
    
    hat_u     = matrix(rep(0, n * m), n, m)   # (p+d)m x 1       
    hat_w     = matrix(rep(0,n*m), n, m)      # n x m
    iter      = 0

    hat_b     = rep(0,m*p)
    hat_b0    = rep(0,p)  
    f         = matrix(0,n,m)
    f_1       = matrix(0,n,m)
    omega     = kronecker(diag(1,p,p),kern)
	xm		  = kronecker(x, t(rep(1,m)))
	pk        = kronecker(rep(1,p),kern)
    xnp       = kronecker(t(x),matrix(1,m,m)) * kronecker(matrix(1,p,n), kern)

    while(iter<iters){
 
        iter   = iter + 1
	
        # update u
        temp   = ((-f + y + hat_w)>(tau/sigma) + 0) * tau/sigma + ((-f + y + hat_w)<(tau-1)/sigma + 0) * (tau-1)/sigma
        hat_u  = (-f + y + hat_w) * (((-f + y + hat_w)>(tau/sigma) | (-f + y + hat_w)<(tau-1)/sigma ) + 0) - temp 
        hat_u  = y - hat_u 
        
        # update b
        uxw    = hat_u + hat_w
        xu     = xnp %*% as.vector(t(uxw))
        hat_b  = ginv_xxkk %*% as.vector(xu)

        # update w
		f_1    = xm %*% (pk * t(kronecker(t(hat_b), rep(1,m))))	
		f      = f_1 
        hat_w  = hat_w + hat_u - f 
		df	   = sum((hat_w>(tau-1))*(hat_w<tau))
		# print(df)

        res     = y - f

    }
	# print(df)

    return (list(d=hat_b,ep_nu=hat_b0,res=res,haty=f,Df = df))

}

quantile_loss<- function(u, tau){
	u * (tau - (u<0))
}

getgram<-function(s,kernel_sigma){
    kernel_sigma = 2*kernel_sigma^2
    Gram = outer(s,s,FUN=function(x,y) exp(-(x-y)^2/kernel_sigma))
    return(Gram)
}

wast_pval <- function(x, y, z, m, tx, tau = 0.5, B=1000, kernel_sigma=0.2){
	n  = nrow(x)
	p  = ncol(x)
	p2 = ncol(tx)
	p3 = ncol(z)

	rho 		= cor(t(z)) - diag(n)
	Omega 		= 0.25 + 0.5*atan(rho/sqrt(1-rho^2))/pi 
	diag(Omega) = 0

	sigma   = 1
	maxIter = 40

	s         = seq(0, 1, length.out = m )
	kern      = getgram(s, kernel_sigma)
	omega1    = kronecker(diag(1,p,p),kern)
	xxkk      = matrix(rep(0,m*p*m*p), m*p, m*p)
    ginv_xx   = ginv(m * (t(x)  %*% x)) 
    for(i in 1:n){
        for (j in 1:m){ 
            xxkk = xxkk + kronecker(((x[i,]) %*% t(x[i,])), ((kern[,j]) %*% t(kern[,j]) ) ) 
        }
    }

	nlamb    = 50
	lamb_a   = seq(2, 7, length.out=nlamb)
	gacv     = rep(0, nlamb)
	for (k in 1:nlamb ){
		lambda     = lamb_a[k]
		ginv_xxkk  = ginv(xxkk + 2*lambda/sigma * omega1 )
		result1    = ADMM_quantile_CP(x, y, m, kernel_sigma, sigma=1, iters=maxIter, tau, ginv_xxkk, kern)
		resids0    = result1$res
		gacv[k]    = sum(quantile_loss(resids0,tau))/(n*m - result1$Df )
	}

	lam_min    = lamb_a[which.min(gacv)]
	ginv_xxkk  = ginv(xxkk + 2*lam_min/sigma * omega1 )
	result     = ADMM_quantile_CP(x, y, m, kernel_sigma, sigma=1, iters=maxIter, tau, ginv_xxkk, kern)
	muhat      = result$haty
	resids0    = result$res

	
	teststat = rep(0,m)
    for(i in 1:m){
        spsi 		= tx*as.numeric((resids0[,i]<0) - tau )
		teststat[i] = sum(diag(t(spsi)%*%Omega%*%spsi))/(n*(n-1))
    }
	Tn0 = sum(as.numeric(teststat))

	#--------------- Bootstrap ---------------------------------------------
	yb 		    = matrix(0, n, B)
	Tn		    = rep(0, B)
	for(k in 1:B){
		set.seed(k)
		wb 		= matrix(sample(c(-2*tau, 2*(1-tau)), size=n*m, prob=c(tau,1-tau), replace=TRUE), n, m)
		yb		= muhat + abs(resids0)*wb

		resultb  = ADMM_quantile_CP(x, yb, m, kernel_sigma, sigma=1,  iters=20, tau, ginv_xxkk, kern)

		for(i in 1:m){
			spsi 		= tx*as.numeric((resultb$res[,i]<0) - tau )
			teststat[i] = sum(diag(t(spsi)%*%Omega%*%spsi))/(n*(n-1))
		}
		Tn[k] = sum(as.numeric(teststat))
	}

	p_value  = mean(Tn > Tn0)
	return(p_value)
}