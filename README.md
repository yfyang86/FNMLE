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
$$f_{\bf X}(\bf x; \bf \mu, \bf \Sigma)=\sum_{\bf s \in \bf S(n)} (2 \pi)^{-\frac{n}{2}} |{\bf \Sigma}|^{-\frac{1}{2}}\exp\{-\frac{1}{2}(\bf\Lambda_{\bf s}^{(n)}{\bf x} - {\bf\mu})^{'} {\bf\Sigma}^{-1}(\bf\Lambda_{\bf s}^{(n)}\bf x - {\bf\mu})\},$$
where 
$${\bf s}=(s_1, \cdots, s_n)\in {\bf S}(n) = \{(s_1, \cdots, s_n): s_i = \pm 1, i = 1,\cdots, n \}$$ 
represents a possible sign vector, and the diagonal sign matrix is ${\bf \Lambda}_{\bf s}^{(n)} = diag(s_1, \cdots, s_n)$. 
We further denote ${\bf X} \sim FN_n({\bf \mu}, {\bf \Sigma})$ for simplicity, where the parameters ${\bf \mu}$ and ${\bf \Sigma}$ are the mean vector and variance matrix of the corresponding $n$-dimensional random vector with multivariate normal distribution $N_n({\bf \mu}, {\bf \Sigma})$.

Note: We assume $\bf \mu \ge \bf 0$ in this package. 



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
