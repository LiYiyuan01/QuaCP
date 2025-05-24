\name{wast_pval}
\alias{wast_pval}
\title{Subgroup Testing for Functional Change-plane Quantile Regression via WAST Method}
\description{
A test method based on Weighted Squared Score Statistic (WAST) is proposed to detect the subgroup 
effect in functional quantile regression. This method significantly reduces the computational 
complexity through the closed-loop test statistic formula.
}
\usage{
wast_pval(x, y, z, m, tx, tau = 0.5, B = 1000, 
     kernel_sigma = 0.2)
}

\arguments{
  \item{x}{Design matrix of scalar covariates}
  \item{y}{Functional response matrix}
  \item{z}{Grouping variable matrix}
  \item{m}{Functional domain grid points vector.}
  \item{tx}{Covariate matrix for testing}
  \item{tau}{Quantile level in (0,1), default=0.5 (median)}
  \item{B}{Number of bootstrap samples, default=500}
  \item{kernel_sigma}{Bandwidth for Gaussian kernel smoothing (>0), default=0.2}
}

\value{
\item{p_value}{P-value for the hypothesis test}
}

\references{
Guan X, Li Y, Liu X, You J (2025). 
"Change-plane analysis in functional response quantile regression." 
\emph{arXiv preprint} arXiv:2503.07332. 

Liu X, Huang J, Zhou Y, Zhang F, Ren P (2024). 
"Efficient subgroup testing in change-plane models." 
\emph{arXiv preprint} arXiv:2408.00602. 
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
issig      = 0 
istest     = 1     
ch         = 1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
theta      = c(-0.5, 0.5)

data = gendata(seed, n, m, p, d, q,
                  param_kind, t_gamma, issig, 't3', tau, istest, ch)

result <- wast_pval(data$x, data$y, data$z, m, data$sub_x, 
              tau = 0.5, B = 500,
              kernel_sigma = 0.2)

}
