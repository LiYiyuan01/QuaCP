# QuaCP
 A Package of "Change-plane Analysis in Functional Response Quantile Regression". In the paper, 
 we consider a change-plane model within the framework of functional response quantile regression, 
 capable of identifying and testing subgroups in non-Gaussian functional responses with scalar predictors. 
The alternating direction method of multipliers algorithm is designed to estimate the function coefficients and grouping parameters in the proposed model, thereby dividing the population into different subgroups.
To further test the existence of subgroups, we develop a weighted average of the squared score test statistic, which has a closed form and reduces computational burden. 
         
# Installation

    #install.packages("devtools")
    library(devtools)
    install_github("LiYiyuan01/QuaCP")

    # or
    #install.packages("remotes")
    library(remotes)
    remotes::install_github("LiYiyuan01/QuaCP") 

# Usage
 - [x] [QuaCP-manual.pdf](https://github.com/LiYiyuan01/QuaCP/blob/master/inst/QuaCP-manual.pdf) ---------- Details of the usage of the package.


# Reference

Guan X, Li Y, Liu X, You J (2025). Change-plane analysis in functional response quantile regression. arXiv:2503.07332. 

Liu X, Huang J, Zhou Y, Zhang F, Ren P (2024). Efficient subgroup testing in change-plane models. arXiv:2408.00602. 

Boyd S, Parikh N, Chu E, Peleato B, Eckstein J (2011). Distributed Optimization and Statistical Learning via the Alternating Direction Method of Multipliers. Foundations and Trends in Machine Learning. 3(1):1-122. 

Yuan M, Cai TT (2010). A Reproducing Kernel Hilbert Space Approach to Functional Linear Regression. Annals of Statistics. 38(6):3412-3444. 

# Development
This R package is developed by Yiyuan Li (liyiyuan@stu.sufe.edu.cn ).




