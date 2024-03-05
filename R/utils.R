# License: Apache 2
# Author: Yifan Yang <yfyang.86@hotmail.com>

library(mvtnorm)

#' @export
#' @title Generate combination of {1, -1} of length n
#' @param n length of the combination
#' @return a matrix of all possible combinations of {1, -1} of length n
permsign <- function(n) {
    if (n < 2) stop("n>=2")
    if (n == 2) {
        return(
            matrix(
                c(1,  1,
                  1, -1,
                 -1,  1,
                 -1, -1),
                byrow = TRUE, ncol = 2));
    }
    ex <- permsign(n - 1)
    return(rbind(cbind(1, ex), cbind(-1, ex)));
}
#' @export
#' 
#' @title pointwise mult-op by row
#' @param mat matrix
#' @param vec vector
#' @return a matrix of pointwise multiplication of each row of mat by vec
prodbyrow <- function(mat, vec) {
        t(t(mat) * vec)
}

#' @export
#' @title Density of multivariate folded normal distribution
#' @param xvec a vector of length n
#' @param mu mean, a vector of length n
#' @param sigma2 variance matrix, a matrix of dimension n x n
#' @return density of multivariate folded normal distribution
permute_f <- function(xvec, mu, sigma2){
    n <- length(mu)
    if (n < 2) stop("dim >= 2")
    sign_vec <- permsign_c(n)
    m <- nrow(sign_vec)
    pdf_val <- 0.
    i_row <- 1
    while (i_row <= m) {
        pdf_val <- pdf_val + dmvnorm(prodbyrow_c(xvec, sign_vec[i_row, ]),
                              mean = mu, sigma = sigma2)
        i_row = i_row + 1
        }
    return(pdf_val);
}


#' mu, covmat for 2 dim only
est_helper <- function(params) {
    u1 <- params[1]
    u2 <- params[2]
    s1 <- params[3]
    s2 <- params[4]
    rho <- params[5]
    list(mu = c(u1, u2),
         sigma2 = matrix(
        c(s1^2, s1 * s2 * rho, s1 * s2 * rho, s2^2),
        2, 2))
}

#' @export
#' @title log-likelihood; mu, covmat for 2 dim only
#' @param params parameters for 2 dim case. params = c(mean1, mean2, sig11, sig12, sig22),
#'               hence length(params) = 5
#' @return log-likelihood
loglikdim2 <- function(params) {
    u1 <- params[1]
    u2 <- params[2]
    s1 <- exp(params[3])
    s2 <- exp(params[4])
    rho <- params[5]
    ttmmpp = est_helper(c(u1, u2, s1, s2, rho))
    return(-sum(log(permute_f(matrix(c(x, y), ncol = 2), ttmmpp$mu, ttmmpp$sigma2))));
}

#' @export
#' @title log-likelihood; mu, covmat by row (upper triangle),  dim >= 2
#' @description log-likelihood; mu, covmat by row (upper triangle),  dim >= 2
#' mu, covmat by row (upper triangle), 100 >= dim >= 2
#' 
#' Generally, method == "chole" is more stable (Non-negative Matrix),
#' while method == "tri" only takes garantee that the matrix is symtric.
#' 
#' @param params length(params) = p + (p+1)*p/2
#' params[1:p] stores the mean vector:
#' 
#' -    mu_1 mu_2 ... mu_p
#' 
#' params[-(1:p)] stores the triangle elements of the covariance matrix 
#' (or the Cholesky Decomposition upper triangle matrix v s.t. Cov = v' v)
#'  
#'     sig_{1, 1}, sig_{1, 2}, ..., sig_{1, p}
#'      
#'     0         , sig_{2,2} , ..., sig_{2, p}
#'      
#'     0 ...
#'      
#'     0 ...                      , sig_{p, p}
#' 
#' @param method "tri" or "chole"
est_helper.G2 <- function(params, method = 'tri') {
    pp = length(params)

    df_f <- Vectorize(function(p){p^2/2 + 3*p/2 - pp})
    p = which(df_f(1:100) == 0)
    mu = params[1:p]
    sigma2 = matrix(0, p, p)
    s = p + 1
    for(i in 1:p){
        for(j in i:p){
            if(method == 'tri'){
                if (i == j){
                    sigma2[i, j] = params[s]/2
                }else{
                    sigma2[i, j] = params[s]
                }
            }else{
                # Cholesky
                 sigma2[i, j] = params[s]
            }
            s = s + 1
        }
    }
     if(method == 'tri'){
            sigmamat = sigma2 + t(sigma2)
            }else{
            sigmamat = t(sigma2) %*% (sigma2)
            }
    
    return(list(p = p, mu = mu, sigma2 = sigmamat));
}

# log-likelihood; mu, covmat by row (upper triangle),  dim >= 2
loglik.G2 <- function(params, inputx) {
    ttmmpp = est_helper.G2(params)
    return(-sum(log(permute_f(matrix(inputx, ncol = ttmmpp$p), ttmmpp$mu, ttmmpp$sigma2))));
}

# log-likelihood; mu, covmat by row (upper triangle),  dim >= 2
loglik.G2cholesky <- function(params, inputx) {
    ttmmpp = est_helper.G2(params, 'chol')
    return(-sum(log(permute_f(matrix(inputx, ncol = ttmmpp$p), ttmmpp$mu, ttmmpp$sigma2))));
}

loglik.G2cholesky <- function(params, inputx) {
  ttmmpp = est_helper.G2(params, 'chol')
  return(-sum(log(permute_f(matrix(inputx, ncol = ttmmpp$p), ttmmpp$mu, ttmmpp$sigma2))));
}

#' @export
#' @title Lower triangle index of a matrix
#' @param dim_p dimension of the matrix
#' @return lower triangle index of a matrix
f_lower_idx <- function(dim_p){
    return(matrix(1:(dim_p^2),
                    byrow = TRUE, dim_p)[lower.tri(
                                        matrix(0, dim_p, dim_p), diag = TRUE)])
}