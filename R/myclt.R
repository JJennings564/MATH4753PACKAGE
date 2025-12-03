#' myclt
#'
#' @param n sample size
#' @param iter number of iterations
#'
#' @returns a histogram of sample means.
#' @export
#'
#' @examples myclt(10, 10000)
myclt <- function(n, iter) {
  y <- runif(n * iter, 0, 5)
  data <- matrix(y, nr = n, nc = iter, byrow = TRUE)
  sm <- apply(data, 2, mean)       # Change sum -> mean
  hist(sm, main = "Histogram of Sample Means")  # Plot means instead
  sm
}
