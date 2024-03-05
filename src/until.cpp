#include <Rcpp.h>

using namespace Rcpp;

//' @export
// [[Rcpp::export]]
NumericMatrix prodbyrow_c(NumericMatrix mat, NumericVector vec){
    /**
    example:
    prodbyrow_c(cbind(runif(10), -runif(10)), c(-1, 1))
    */
    NumericMatrix result(mat.nrow(), mat.ncol());
    int i = 0;
    int j = 0;
    for(; i < mat.nrow(); ++i){
        for(; j < mat.ncol(); ++j){
            result(i, j) = (vec(j) < 0 ? -mat(i, j) : mat(i, j));
        }
        j = 0;
    }
    return result;
}

//' @export
// [[Rcpp::export]]
NumericMatrix permsign_c(int n){
    if(n == 2){
        NumericMatrix result(4, 2);
        result(0, 0) = 1; result(0, 1) = 1;
        result(1, 0) = 1; result(1, 1) = -1;
        result(2, 0) = -1; result(2, 1) = 1;
        result(3, 0) = -1; result(3, 1) = -1;
        return (result);
    }
    NumericMatrix tmp = permsign_c(n-1);
    int v = (2<<(n-2));
    int u = (2<<(n-1));
    NumericMatrix result(u, n);
    for(int i = 0; i < v; i++){
        result(i, 0) = 1;
        for(int j = 1; j < n; j++){
            result(i, j) = tmp(i, j-1);
        }
    }
    
    for(int i = v; i < u; i++){
        result(i, 0) = -1;
        for(int j = 1; j < n; j++){
            result(i, j) = tmp(i-v, j-1);
        }
    }

    return result;
}