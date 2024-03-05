UNITEST_FN_MLE <- function() {
    data <- matrix(c(
        1.96424772, 0.05716957,
        1.14657676, 1.01913463,
        0.08995317, 0.59953389,
        0.35224735, 1.34653609,
        0.02756375, 0.50781523,
        0.38535959, 1.32675296,
        0.13432975, 1.59026667,
        0.38934246, 1.22333517,
        0.20953879, 0.18159422,
        1.02054088, 1.89612705
    ), ncol = 2)
    result <- FN_MLE(data)
    return(sprintf("%0.3f", result[["par"]][1]))
}

# Compare this snippet from R/FN_MLE.R:
test_that("FN_MLE works", {
    expect_equal(UNITEST_FN_MLE(), "0.294")
})
