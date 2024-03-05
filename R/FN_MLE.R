#' @useDynLib FNMLE, .registration=TRUE

#' @export
#' @title Estimation of Multivariate Folded normal distribution
#'
#' @param data: the matrix of normal distribution, A multidimensional normal distribution matrix with n rows and
#'              p columns, n is the sample size, and p is the dimension of the distribution.
#' @param hessian: If return Hessian matrix. Default is TRUE.
#'
#' @return par:    The best set of parameters found.
#'         value:  The value of folded normal density function corresponding to par.
#'         counts: A two-element integer vector giving the number of calls to fn and gr respectively. This excludes
#'         those calls needed to compute the Hessian, if requested, and any calls to fn to compute a finite-difference
#'         approximation to the gradient.
#'         convergence: An integer code. 0 indicates successful completion, 1 indicates that the iteration limit
#'         maxit had been reached.
#'         message: A character string giving any additional information returned by the optimizer, or NULL.
#'         hessian: Only if argument hessian is true. A symmetric matrix giving an estimate of the Hessian at
#'         the solution found. Note that this is the Hessian of the unconstrained problem even if the box
#'         constraints are active.
#'
#' @examples{
#' data <- matrix(abs(rnorm(100)), ncol = 4)
#' FN_MLE(data)
#' }
FN_MLE <- function(data, hessian = TRUE) {
  dim_p <- ncol(data)
  lower_idx <- f_lower_idx(dim_p)
  inputx <- abs(data)
  sss <- cov(data)
  init <- c(
    colMeans(data),
    sss[lower_idx]
  )

  loglik_hack <- function(params) {
    return(loglik.G2(params, inputx))
  }

  fit <- optim(init, loglik_hack, method = "BFGS", hessian = hessian)
  fit
}
