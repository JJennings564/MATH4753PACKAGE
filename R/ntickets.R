#' Compute number of tickets to sell given seat count, show probability, and target overbook probability
#'
#' @param N integer number of seats on the plane
#' @param gamma probability that the plane is truly overbooked
#' @param p probability that a single purchaser shows up
#' @return A named list with elements:
#'   \item{nd}{integer found using exact discrete binomial calculation}
#'   \item{nc}{integer found using normal approximation (continuity corrected)}
#'   \item{N}{input N}
#'   \item{p}{input p}
#'   \item{gamma}{input gamma}
#' @details
#' The function creates two separate plots of the objective function:
#' The solution is the smallest integer n where the objective crosses zero.
#' @examples
#' ntickets(N = 400, gamma = 0.02, p = 0.95)
#' @export
#'
ntickets <- function(N, gamma, p) {
  # N = number of seats on the flight
  # gamma = probability of overbooking (more people show than seats)
  # p = probability that a passenger shows up

  # 1: Discrete distribution Binomial
  n_start <- N
  n_max <- (N * 1.1)
  nd <- NA
  for (n in n_start:n_max) {
    objective <- 1 - gamma - pbinom(N, n, p)
    if (objective > 0) {
      nd <- n  # first n where objective > 0
      break
    }
  }

  # 2: Normal approximation
  q <- 1 - p
  z <- qnorm(1 - gamma)


  nc <- N / p  # Initial guess
  for (i in 1:10) {
    nc <- (N + 0.5 - z * sqrt(nc * p * q)) / p
  }

  if (is.na(nd)) nd <- n_max


  # plots of objective functions
  n_vals <- seq(N, ceiling(max(nd, nc, na.rm = TRUE)) + 15, by = 0.1)

  # Discrete objective function
  n_discrete <- seq(N, ceiling(max(nd, nc, na.rm = TRUE)) + 15, by = 1)
  obj_discrete <- sapply(n_discrete, function(n) {
    1 - gamma - pbinom(N, n, p)
  })

  # Continuous objective function
  obj_continuous <- sapply(n_vals, function(n) {
    mu <- n * p
    sigma <- sqrt(n * p * (1 - p))
    1 - gamma - pnorm(N + 0.5, mu, sigma)
  })

  par(mfrow = c(2, 1), mar = c(4, 4, 3, 2))


  plot(n_discrete, obj_discrete, type = "b",
       xlab = "n",
       ylab = "Objective",
       main = sprintf("Objective Vs n to find optimal tickets sold\n(%d) gamma= %.2f N=%d discrete",
                      nd, gamma, N),
       col = "blue", pch = 19, cex = 0.8,
       ylim = c(0, 1))
  abline(h = 0, col = "red", lty = 1, lwd = 2)
  abline(v = nd, col = "blue", lty = 1, lwd = 2)
  grid()


  plot(n_vals, obj_continuous, type = "l",
       xlab = "n",
       ylab = "Objective",
       main = sprintf("Objective Vs n to find optimal tickets sold\n(%.3f) gamma= %.2f N=%d continuous",
                      nc, gamma, N),
       col = "black", lwd = 1.5,
       ylim = c(0, 1))
  abline(h = 0, col = "red", lty = 1, lwd = 2)
  abline(v = nc, col = "blue", lty = 1, lwd = 2)
  grid()

  par(mfrow = c(1, 1))

  result <- list(nd = nd, nc = nc, N = N, p = p, gamma = gamma)
  print(result)

  return(invisible(result))
}
