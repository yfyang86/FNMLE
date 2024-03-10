# FNMLE

| Contents | Logs |
|:---------|:-----|
| Date | 2024-3-05 | 
| Version | v-0.1 |

This package provides an numerical approach to derive the maximum likelihood estimation (MLE) of multidimensional folded normal distribution (FN).

# An introduction to a multidimensional folded normal distribution
First, let's recall the definition of the $n$-dimensional folded normal random vector. More details can be found in Chakraborty and Chatterjee (2013), Xi Liu, Yiqiao Jin, Yifan Yang and Xiaoqing Pan (2023).

A random vector ${\bf X} = (X_1,\cdots,X_n)^{'}$ is said to have a multivariate folded normal distribution with a real vector ${\bf\mu}\in\mathbb{R}^n$ and a symmetric positive definite matrix $\bf \Sigma_{n \times n}$,
if its probability density function is given by
$$f_{\bf X}(\bf x; \bf \mu, \bf \Sigma)=\sum_{\bf s \in \bf S(n)} (2 \pi)^{-\frac{n}{2}} |{\bf \Sigma}|^{-\frac{1}{2}}\exp[{-\frac{1}{2}(\bf\Lambda_{\bf s}^{(n)}{\bf x} - {\bf\mu})^{'} {\bf\Sigma}^{-1}(\bf\Lambda_{\bf s}^{(n)}\bf x - {\bf\mu})}],  \  {\bf x \geq \bf 0},$$
where 
$${\bf s}=(s_1, \cdots, s_n)\in {\bf S}(n) = [{(s_1, \cdots, s_n): s_i = \pm 1, i = 1,\cdots, n }]$$ 
represents a possible sign vector, and the diagonal sign matrix is ${\bf \Lambda}_{\bf s}^{(n)} = diag(s_1, \cdots, s_n)$. 
We further denote ${\bf X} \sim FN_n({\bf \mu}, {\bf \Sigma})$ for simplicity, where the parameters ${\bf \mu}$ and ${\bf \Sigma}$ are the mean vector and variance matrix of the corresponding $n$-dimensional random vector with multivariate normal distribution $N_n({\bf \mu}, {\bf \Sigma})$.

Note: We assume $\bf \mu \ge \bf 0$ in this package. 


# An example of the bivariate case
In this section, we will introduce an example of the bivariate case $FN_2(\bf \mu, \bf \Sigma)$, where $\bf \mu = (\mu_1, \mu_2)'$ and 
$$\bf \Sigma = \begin{array}{cc} \sigma_{11} & \sigma_{21} \\ \sigma_{21} & \sigma_{22} \end{array}.$$



## Simulation Data Generation

We began by generating 100 samples from a bivariate folded normal distribution with mean parameters $\bf \mu=(4, 6)$ and the covariance matrix $\bf \Sigma$ given by:
$\begin{bmatrix}
	1 & 0.4 \\
	0.4 & 4
\end{bmatrix}.$


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
sigma2 <- matrix(c(sig1^2,  rho * sig1 * sig2,
                rho * sig1 * sig2,  sig2^2), nrow = 2)
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

The estimated covariance matrix is:
```r
### The estimated covariance matrix
est_helper.G2(result$par)[["sigma2"]]
```

These results highlight the accuracy of our estimation method in capturing the underlying parameters of the bivariate folded normal distribution based on the given samples.

## The estimated coverage probability of the 95% confidence intervals for the parameters
Additional simulations are performed to assess the estimated coverage probability of the $95\%$ confidence intervals for the parameters.
Following the approach of Tsagris et al. (2014), $95\%$ confidence intervals are calculated using the normal approximation, where the variance matrix is estimated using the inverse of the Hessian matrix.
```r
set.seed(123456)
simu_n <- 1000
dim_p <- 2
result_df <- data.frame()
n <- 20
mu <- c(2.5, 2.5)
ss <- matrix(c(25, 5, 5, 25), 2, 2)
lower_idx <- f_lower_idx(dim_p)
for (i in seq_len(simu_n)) {
    dat <- mvrnorm(n, mu, ss)
    inputx <- abs(dat)
    fit <- FN_MLE(inputx)
    fisher_info <- solve(fit$hessian)
    prop_sigma <- sqrt(diag(fisher_info))
    prop_sigma <- diag(prop_sigma)
    a <- diag(prop_sigma)
    prop_sigma <- a
    upper <- fit$par + 1.96 * prop_sigma
    lower <- fit$par - 1.96 * prop_sigma
    true_para <- c(mu, ss[lower_idx])
    p <- c()
    for (k in 1:5) {
        p[k] <- (lower[k] <= true_para[k] & true_para[k] <= upper[k])
    }
    result_df <- rbind(result_df, p)
}
## calculate coverage rate of parameters ##
# Be sure to handle NA's before the
# coverage rate calculation
result_df[is.na(result_df)] <- FALSE
p <- apply(result_df, 2, mean)
names(p)=c("mu1","mu2","sigma11","sigma21","sigma22")
round(p, 2)
```

Then the estimated coverage rate the parameters are:
| n | simu_n   |   mu1   |     mu2     |   sigma11   | sigma21   | sigma22 |
|--:|---------:|--------:|------------:|------------:|----------:|--------:|
|20 |1000      |0.68     |0.67         |0.71         |0.75       |0.69     |

Here:
- n: samples size
- simu_n: simulation times



