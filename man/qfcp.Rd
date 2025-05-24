\name{qfcp}
\alias{qfcp}
\title{Functional Change-plane Quantile Regression via ADMM}
\description{
  Implements Alternating Direction Method of Multipliers (ADMM) for functional change-plane quantile regression models. This method handles  quantile regression with subgroup structures  via RKHS.
}

\usage{
qfcp(x, y, z, sub_x, m, kernel_sigma, sigma, lambda, tau)
}

\arguments{
  \item{x}{Design matrix (n x p) of scalar covariates.}
  \item{y}{Functional response matrix (n x m).}
  \item{z}{Grouping variable matrix (n x q), second column is intercept.}
  \item{sub_x}{Threshold covariate subset matrix (n x d).}
  \item{m}{Number of functional domain points. }
  \item{kernel_sigma}{Bandwidth parameter for Gaussian kernel construction.}
  \item{sigma}{ADMM penalty parameter. }
  \item{lambda}{Regularization parameter. }
  \item{tau}{Quantile level in (0,1).}
}

\value{
\item{gamma}{Estimated subgroup parameters (length q-1).}
\item{d}{Coefficient functions (m(p+d) x 1 vector).}
\item{w}{Dual variables matrix (n x m).}
}

\references{
Guan X, Li Y, Liu X, You J (2025). 
"Change-plane analysis in functional response quantile regression." 
\emph{arXiv preprint} arXiv:2503.07332. 

Boyd S, Parikh N, Chu E, Peleato B, Eckstein J (2011). 
"Distributed Optimization and Statistical Learning via the Alternating Direction Method of Multipliers." 
\emph{Foundations and Trends in Machine Learning} 3(1):1-122. 

Yuan M, Cai TT (2010). 
"A Reproducing Kernel Hilbert Space Approach to Functional Linear Regression." 
\emph{Annals of Statistics} 38(6):3412-3444. 
}

\examples{
seed = 10086
n    = 200
m    = 30   
p    = 3
q    = 2
d    = 2    
tau  = 0.5      
t_gamma    = c(-1, 1)               
param_kind = c(3,1,2,4,5)                
issig      = 1 
istest     = 0                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
theta      = c(-0.5, 0.5)

data = gendata(seed, n, m, p, d, q,
                  param_kind, t_gamma, issig, 't3', tau, istest)

fit1 <- qfcp(data$x, data$y, data$z, data$sub_x, m, kernel_sigma=0.5,
             sigma=1, lambda=0.01, tau=0.5)


}
