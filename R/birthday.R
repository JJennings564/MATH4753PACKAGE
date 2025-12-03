#' Birthday Problem Probability
#'
#' Computes the probability of at least two people in a group with size n sharing the same birthday.
#'
#' @param x An integer or numeric vector giving the group size(s).
#'
#' @return A numeric vector of probabilities.
#'
#' @examples
#' birthday(20:25)
#' @export
birthday <- function(x){

  1 - exp(lchoose(365,x) + lfactorial(x) - x * log(365))

}
