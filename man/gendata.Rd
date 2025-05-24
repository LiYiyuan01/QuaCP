\name{gendata}
\alias{gendata}
\title{
  Generate Simulated Data for Functional Change-plane Quantile Regression Model.
}

\description{
  Generates simulated datasets for functional change-plane quantile regression models. 
  This function creates synthetic data with customizable parameters including sample size, 
  covariate dimensions, and error distributions.
}

\usage{
gendata(seed, n, m, p, d, q, param_kind, gamma, is_sig, err_kind, tau, istest=0, ch=1)
}

\arguments{
  \item{seed}{Random seed for reproducibility (integer).}
  \item{n}{Sample size (positive integer).}
  \item{m}{Number of observation points in functional domain (positive integer).}
  \item{p}{Dimension of covariates x, the second colum is 1 (positive integer, here p must be less than 6).}
  \item{d}{Dimension of covariates sub_x which the subset of x (positive integer).}
  \item{q}{Dimension of grouping variable (positive integer).}
  \item{param_kind}{Kind of slope function(There are 1 to 5 function options available, and each number contains one function).}
  \item{gamma}{Grouping parameter vector (length q-1).}
  \item{is_sig}{Threshold effect activation flag (0=disable, 1=enable, default=0).}
  \item{err_kind}{Error distribution: "gaussian", "t3", "weibull" or "laplace" .}
  \item{tau}{Quantile level.}
  \item{istest}{Simulaton for estimation (\eqn{0}) or testing (\eqn{1}).}
  \item{ch}{Signal strength for testing.}
  
}

\value{
  A list containing:
  \item{y}{Function response matrix (n x m).}
  \item{x}{Scalar covariate matrix (n x p).}
  \item{sub_x}{Scalar covariate matrix, the subset of X (n x d).}
  \item{z}{Grouping variable (n x q).}
  \item{beta}{True slope functions of x. }
  \item{deta}{True slope functions of sub_x.}
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

sim_data = gendata(seed, n, m, p, d, q,
                  param_kind, t_gamma, issig, 't3', tau, istest)

# Inspect structure
str(sim_data)

}
