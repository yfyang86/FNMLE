# FNMLE

| Contents | Logs |
|:---------|:-----|
| Date | 2024-3-05 | 
| Version | v-0.1 |

This package provides an numerical approach to derive the maximum likelihood estimation (MLE) of multidimensional folded normal distribution (FN).

Will add the definition of multidimensional FN here....

Note: We assume $\mu_i \ge 0$ in this package. 

# A Simple Walk through

In this section we will illustrate the two dimensional case of the folded normal distribution. 

Will illustrate Sigma matrix here...

## Simulation Data Generation


```r
library(MASS)
set.seed(123)
n <- 100
u1 <- 4
u2 <- 6
sig1 <- 1
sig2 <- 2
rho <- 0.2
mean <- c(u1, u2)
sigma2 <- matrix(c(
                sig1^2,             rho * sig1 * sig2,
                rho * sig1 * sig2,  sig2^2),
                2, 2)
dat <- abs(mvrnorm(n, mean, sigma2))
```
## The MLE of multidimensional FN distribution
Applying the function "FN_MLE" in our R package "FNMLE", the MLE of multidimensional FN is derived.
```r
## MLE
result <- FN_MLE(dat)
print(result$par)
names(result$par)=c("mu1","mu2","sigma11","sigma21","sigma22")
```
