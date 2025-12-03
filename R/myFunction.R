#' Predict spruce height from breast diameter
#'
#' This function uses coefficients from a quadratic regression model
#' to predict spruce tree height based on breast diameter.
#'
#' @param x Numeric vector of breast height diameters in cm.
#' @param coefs Numeric vector of coefficients (default = NULL).
#' If NULL, the function will look for quad.lm in the global environment.
#' @return Numeric vector of predicted heights in m.
#' @examples
#' x <- c(10, 12, 14, 16, 18)
#' y <- 2 + 0.5 * x + 0.01 * x^2
#' quad.lm <- lm(y ~ poly(x, 2, raw = TRUE))

#' predict_spruce(c(15, 18, 20), coefs = coef(quad.lm))
#'
#' @export
predict_spruce <- function(x, coefs = NULL) {
  if (is.null(coefs)) {
    if (!exists("quad.lm", envir = .GlobalEnv)) {
      stop("You must supply coefficients or have quad.lm fitted in the global environment.")
    }
    coefs <- stats::coef(get("quad.lm", envir = .GlobalEnv))
  }
  coefs[1] + coefs[2]*x + coefs[3]*x^2
}
