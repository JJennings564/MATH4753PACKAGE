#' myncurve Function
#'
#' This function draws a normal curve and shades the area
#' between two x-values. It also prints the probability value.
#'
#' @param mu The mean of the normal distribution.
#' @param sigma The standard deviation of the normal distribution.
#' @param a The starting x-value for shading.
#'
#' @return The probability between a and b.
#' @export
#'
#' @examples
#' myncurve(mu = 10, sigma = 4, a = 10)

myncurve <- function(mu, sigma, a) {
  # Set up the x and y values for the normal curve
  x <- seq(mu - 4*sigma, mu + 4*sigma, length = 1000)
  y <- dnorm(x, mean = mu, sd = sigma)

  # Plot the curve
  plot(x, y, type = "l", lwd = 2,
       main = paste("Normal Curve: mu =", mu, ", sigma =", sigma))

  # Define shaded region from left tail to x=a
  xshade <- seq(mu - 4*sigma, a, length = 100)
  yshade <- dnorm(xshade, mean = mu, sd = sigma)

  # Fill the shaded area
  polygon(c(min(xshade), xshade, a), c(0, yshade, 0), col = "lightblue")

  # Compute probability
  prob <- pnorm(a, mean = mu, sd = sigma)
  prob <- round(prob, 4)

  # Display area on graph
  text(a, max(y)/4, paste("P(X <=", a, ") =", prob))

  # Return a list with results
  return(list(mu = mu, sigma = sigma, area = prob))
}
