#' @importFrom MASS mvrnorm 
#' @importFrom stats rchisq rnorm rweibull
#' @importFrom VGAM rlaplace
#' @export

x_covariance <-function(s){
  covar = matrix(0,s,s)
  for (i in 1:s){
    for (j in 1:s) {
      covar[i,j] = 0.5^(abs(i-j))
      }
  }
  return(covar)
}

e_covariance <-function(s){
  a = 1
  h = 0.8
  covar = matrix(0,s,s)
  for (i in 1:s){
    for (j in 1:s) {
      covar[i,j] = a*exp(-((i-j)/h)^2)
    }
  }
  return(covar)
}

rmv_t <- function(n,df,mu,Sigma){
	WW = sqrt(df/rchisq(n,df))
	LL = t(chol(Sigma))
	Z  = matrix(rnorm(n*length(mu)),ncol=n)
	return(mu+ t(WW*t(LL%*%Z)))
}

rmv_lapla <- function(n,mu,Sigma){
	LL = t(chol(Sigma))
	Z  = matrix(rlaplace(n*length(mu),0,1),ncol=n)
	return(t(LL%*%(Z)))
}

paramater_function <- function(kind,s){

    if(kind == 1){
        paramater = (1-s)^3
    }
    else if(kind == 2){
        paramater =  exp(-3*s)
    }
    else if(kind == 3){
        paramater =  sin(pi*s) + s^3
    }
    else if(kind == 4){
        paramater =  4*cos(pi*s/2)+3*s^3
    }
    else if(kind == 5){
        paramater =  3*s^2 + 3
    }
    else if(kind == 6){
        paramater =  sqrt(2)*cos(2*pi*s)
    }
    else if(kind == 0){
        paramater =  0*s
    }
    else              {
         stop("wrong!")
        # paramater = sqrt(2)*cos(2*pi*s)
    }

    return (paramater)
}

gendata <- function(seed, n, m, p, d, q, param_kind, gamma, is_sig, err_kind, tau, istest=0, ch=1) {

    set.seed(seed)
    s = seq(0, 1, length.out = m )
    x   =  mvrnorm(n, rep(0, (p-1)), x_covariance(p-1))
    x       = cbind(rep(1,n), x)

	err = switch(err_kind,
				'gaussian'	= mvrnorm(n, rep(0, m), e_covariance(m)),
				't3'        = t(rmv_t(n,3,rep(0,m),e_covariance(m))),
                'laplace'   = rmv_lapla(n,rep(0,m),e_covariance(m)),
				'weibull'   = matrix(rweibull(n*m, shape = 0.3, scale = 0.5),n,m)
	)

    beta  = matrix(0,m,p)
    if(istest == 0){
        z_1   = rnorm(n, 0, 1)
        z_2   = cbind(rep(1,n), matrix(rnorm(n*(q-1),1,1),n,(q-1)))
        z     = cbind(z_1, z_2)
        zg    = z_2 %*% gamma
        ind   = ((z_1 + zg )>0)           
    }
    else if (istest == 1) {
        z1    = matrix(rnorm(n*q, 0, 1),n,q)
        zg    = z1 %*% gamma
        z1g   = quantile(zg, 0.35)
        ind   = (zg > z1g) * ch
        z     = cbind(1,z1)
    }
    for (i in 1:p){    
        beta[,i] = paramater_function(param_kind[i],s)
    }

    if (is_sig == 1){
        sub_x = x[ , (p-d+1):p] 
        # sub_x = x[ , 2:3] 
        # sub_x = matrix(rnorm(n*d, mean=0, sd=1), n, d)
        deta  = cbind(paramater_function(param_kind[p+1], s), paramater_function(param_kind[p+2], s))
        sub   = sub_x * as.vector(ind)
        y_nm  = (x %*% t(beta)) + (sub %*% t(deta)) + err
        # cat("n",n,"sum",sum(ind),"\n")

    }
    else{
        sub_x = matrix(rep(0,n*d),n,d) 
        deta  = matrix(rep(0,m*d),m,d)        
        y_nm  = (x %*% t(beta))  +  err
    }
    # loss = sum(err*(tau-(err<0)))/(n*m)

	return( list(y=y_nm, x=x, z=z, sub_x=sub_x, deta = deta,  beta = beta))
	
}
