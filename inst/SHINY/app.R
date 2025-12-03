# app.R
library(shiny)
library(ggplot2)

# Helper: transform/untransform parameters and negative log-likelihoods.
# We optimize on unconstrained theta; map to actual params via exp or plogis.

# Negative log-likelihood for each distribution in terms of unconstrained theta
nll_normal <- function(theta, x) {
  mu <- theta[1]            # unconstrained
  sigma <- exp(theta[2])    # positive
  n <- length(x)
  -sum(dnorm(x, mean = mu, sd = sigma, log = TRUE))
}
nll_exponential <- function(theta, x) {
  rate <- exp(theta[1])
  -sum(dexp(x, rate = rate, log = TRUE))
}
nll_poisson <- function(theta, x) {
  lambda <- exp(theta[1])
  -sum(dpois(x, lambda = lambda, log = TRUE))
}
nll_binomial <- function(theta, x, size) {
  p <- plogis(theta[1]) # in (0,1)
  -sum(dbinom(x, size = size, prob = p, log = TRUE))
}
nll_gamma <- function(theta, x) {
  shape <- exp(theta[1])
  rate  <- exp(theta[2])
  -sum(dgamma(x, shape = shape, rate = rate, log = TRUE))
}
nll_beta <- function(theta, x) {
  a <- exp(theta[1])
  b <- exp(theta[2])
  -sum(dbeta(x, shape1 = a, shape2 = b, log = TRUE))
}

# Function to run an optim-based MLE and return estimates + se via Hessian + delta method
mle_optim <- function(nll_fun, theta_start, names_out, transform_fun, x, extra = list()) {
  # nll_fun: function(theta, x, ...) returns negative log-likelihood
  # theta_start: numeric vector initial values on unconstrained scale
  res <- tryCatch(
    optim(theta_start, fn = function(t) nll_fun(t, x, !!!extra), method = "BFGS", hessian = TRUE),
    error = function(e) list(par = NA, value = NA, convergence = 99, message = e$message),
    warning = function(w) invokeRestart("muffleWarning")
  )
  if (!is.list(res) || is.null(res$par)) return(list(success = FALSE, message = "optim failed"))
  theta_hat <- res$par
  # map to original parameter scale and compute Jacobian for delta method
  par_hat <- transform_fun(theta_hat)
  # compute approximate covariance on theta scale
  if (!is.null(res$hessian) && !any(is.na(res$hessian))) {
    cov_theta <- tryCatch(solve(res$hessian), error = function(e) NULL)
  } else cov_theta <- NULL
  se <- rep(NA, length(par_hat))
  if (!is.null(cov_theta)) {
    # numeric Jacobian matrix J_ij = d g_i(theta) / d theta_j
    eps <- 1e-6
    J <- matrix(0, nrow = length(par_hat), ncol = length(theta_hat))
    for (j in seq_along(theta_hat)) {
      thp <- theta_hat
      thp[j] <- thp[j] + eps
      g1 <- transform_fun(thp)
      thm <- theta_hat
      thm[j] <- thm[j] - eps
      g2 <- transform_fun(thm)
      J[, j] <- (g1 - g2) / (2 * eps)
    }
    cov_par <- J %*% cov_theta %*% t(J)
    se <- sqrt(pmax(0, diag(cov_par)))
  }
  
  list(success = TRUE,
       par = par_hat,
       se = se,
       logLik = -res$value,
       convergence = res$convergence,
       message = ifelse(is.null(res$message), "", res$message))
}

# transform functions for each distribution
transform_normal <- function(theta) c(mu = theta[1], sigma = exp(theta[2]))
transform_exponential <- function(theta) c(rate = exp(theta[1]))
transform_poisson <- function(theta) c(lambda = exp(theta[1]))
transform_binomial <- function(theta) c(p = plogis(theta[1]))
transform_gamma <- function(theta) c(shape = exp(theta[1]), rate = exp(theta[2]))
transform_beta  <- function(theta) c(alpha = exp(theta[1]), beta = exp(theta[2]))

# UI
ui <- fluidPage(
  titlePanel("MLE demonstration for univariate distributions"),
  sidebarLayout(
    sidebarPanel(
      helpText("Simulate data or upload your own. Then choose a distribution and click 'Run MLE'."),
      selectInput("dist", "Distribution",
                  choices = c("Normal", "Exponential", "Poisson", "Binomial", "Gamma", "Beta"),
                  selected = "Normal"),
      numericInput("n", "Sample size (for simulation)", value = 200, min = 2),
      hr(),
      h4("Simulation parameters"),
      conditionalPanel(
        "input.dist == 'Normal'",
        numericInput("norm_mu", "True mu", value = 0),
        numericInput("norm_sigma", "True sigma", value = 1, min = 1e-6)
      ),
      conditionalPanel(
        "input.dist == 'Exponential'",
        numericInput("exp_rate", "True rate", value = 1, min = 1e-6)
      ),
      conditionalPanel(
        "input.dist == 'Poisson'",
        numericInput("pois_lambda", "True lambda", value = 3, min = 0)
      ),
      conditionalPanel(
        "input.dist == 'Binomial'",
        numericInput("binom_size", "Number of trials (n) for each observation", value = 10, min = 1),
        numericInput("binom_p", "True p", value = 0.3, min = 0, max = 1)
      ),
      conditionalPanel(
        "input.dist == 'Gamma'",
        numericInput("gamma_shape", "True shape (k)", value = 2, min = 1e-6),
        numericInput("gamma_rate", "True rate (theta)", value = 1, min = 1e-6)
      ),
      conditionalPanel(
        "input.dist == 'Beta'",
        numericInput("beta_a", "True alpha", value = 2, min = 1e-6),
        numericInput("beta_b", "True beta", value = 5, min = 1e-6)
      ),
      actionButton("sim", "Simulate data"),
      hr(),
      fileInput("file", "Upload CSV (single column numeric data, or two columns for binomial: successes, trials)", accept = ".csv"),
      checkboxInput("header", "Has header (if uploading CSV)", TRUE),
      hr(),
      actionButton("run", "Run MLE"),
      width = 3
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Plot & Fit", plotOutput("fitPlot", height = "450px"), br(),
                 verbatimTextOutput("fitSummary")),
        tabPanel("Data", tableOutput("dataTable")),
        tabPanel("About", 
                 h4("What this app does"),
                 p("This app demonstrates maximum likelihood estimation (MLE) for several univariate distributions."),
                 p("MLE is done by numeric optimization of the negative log-likelihood using optim()."),
                 p("Constrained parameters (e.g., rates > 0 or probabilities in (0,1)) are optimized on an unconstrained scale and mapped back with exp() or plogis().")
        )
      )
    )
  )
)

# Server
server <- function(input, output, session) {
  # Reactive stored data
  data_r <- reactiveVal(NULL)
  
  observeEvent(input$sim, {
    n <- as.integer(input$n)
    dist <- input$dist
    x <- switch(dist,
                "Normal" = rnorm(n, mean = input$norm_mu, sd = input$norm_sigma),
                "Exponential" = rexp(n, rate = input$exp_rate),
                "Poisson" = rpois(n, lambda = input$pois_lambda),
                "Binomial" = rbinom(n, size = input$binom_size, prob = input$binom_p),
                "Gamma" = rgamma(n, shape = input$gamma_shape, rate = input$gamma_rate),
                "Beta"  = rbeta(n, shape1 = input$beta_a, shape2 = input$beta_b)
    )
    data_r(x)
  })
  
  observeEvent(input$file, {
    req(input$file)
    df <- read.csv(input$file$datapath, header = input$header, stringsAsFactors = FALSE)
    # If single column numeric -> vector
    if (is.data.frame(df) && ncol(df) == 1) {
      if (!is.numeric(df[[1]])) {
        showNotification("Uploaded data is not numeric", type = "error")
      } else data_r(df[[1]])
    } else if (is.data.frame(df) && ncol(df) >= 2 && input$dist == "Binomial") {
      # for binomial allow two columns: successes, trials OR three columns? We'll assume first two
      x_success <- as.numeric(df[[1]])
      x_trials  <- as.numeric(df[[2]])
      if (any(is.na(x_success)) || any(is.na(x_trials))) {
        showNotification("Binomial upload: numeric columns required", type = "error")
      } else {
        # store as matrix-like with successes and trials as list
        data_r(list(success = x_success, trials = x_trials))
      }
    } else {
      # fallback: try to extract numeric columns
      nums <- sapply(df, is.numeric)
      if (any(nums)) {
        data_r(unlist(df[, nums, drop = FALSE][1]))
      } else {
        showNotification("Uploaded CSV didn't contain usable numeric column(s).", type = "error")
      }
    }
  })
  
  output$dataTable <- renderTable({
    x <- data_r()
    if (is.null(x)) return(NULL)
    if (is.list(x) && !is.null(x$success)) {
      data.frame(successes = x$success, trials = x$trials)
    } else data.frame(x = x)
  })
  
  # Run MLE when user clicks
  fit_res <- eventReactive(input$run, {
    x_raw <- data_r()
    if (is.null(x_raw)) {
      showNotification("No data: simulate or upload data first", type = "error")
      return(NULL)
    }
    dist <- input$dist
    # prepare x appropriately for each dist
    if (dist == "Binomial" && is.list(x_raw) && !is.null(x_raw$success)) {
      x <- x_raw$success
      size <- x_raw$trials
    } else if (dist == "Binomial" && is.numeric(x_raw)) {
      x <- x_raw
      size <- input$binom_size
    } else if (is.list(x_raw) && !is.null(x_raw$success) && dist != "Binomial") {
      # user uploaded successes/trials but selected other dist
      showNotification("Uploaded two-column data interpreted as binomial; switch to Binomial or upload single column", type = "error")
      return(NULL)
    } else {
      x <- as.numeric(x_raw)
      size <- NULL
    }
    # preliminary checks per distribution
    if (dist %in% c("Exponential", "Gamma", "Beta") && any(x <= 0) && dist != "Beta") {
      if (dist == "Exponential" || dist == "Gamma") {
        if (any(x <= 0)) {
          showNotification("Data must be strictly positive for Exponential/Gamma", type = "error")
          return(NULL)
        }
      }
    }
    if (dist == "Beta" && any(x <= 0 | x >= 1)) {
      showNotification("Data for Beta must be in (0,1)", type = "error")
      return(NULL)
    }
    
    # run appropriate MLE
    if (dist == "Normal") {
      theta0 <- c(mean(x), log(sd(x) + 1e-6))
      res <- mle_optim(nll_normal, theta0, c("mu", "sigma"), transform_normal, x)
    } else if (dist == "Exponential") {
      theta0 <- c(log(1 / (mean(x) + 1e-6)))
      res <- mle_optim(nll_exponential, theta0, c("rate"), transform_exponential, x)
    } else if (dist == "Poisson") {
      theta0 <- log(mean(x) + 1e-6)
      res <- mle_optim(nll_poisson, c(theta0), "lambda", transform_poisson, x)
    } else if (dist == "Binomial") {
      # x are successes, size is known (could be scalar)
      if (is.null(size)) {
        showNotification("Binomial: supply trials (size)", type = "error")
        return(NULL)
      }
      # if size is a vector of length 1, replicate
      if (length(size) == 1) size <- rep(size, length(x))
      theta0 <- qlogis((sum(x) + 0.5) / (sum(size) + 1)) # logit of observed proportion
      # we need a wrapper nll that passes size
      nll_binom_wrap <- function(theta, x, size) nll_binomial(theta, x, size)
      res <- mle_optim(nll_binom_wrap, c(theta0), "p", transform_binomial, x, extra = list(size = size))
    } else if (dist == "Gamma") {
      # use method-of-moments start
      m <- mean(x); v <- var(x)
      shape0 <- ifelse(v > 0, m^2 / v, 1)
      rate0  <- ifelse(v > 0, m / v, 1)
      theta0 <- c(log(shape0 + 1e-6), log(rate0 + 1e-6))
      res <- mle_optim(nll_gamma, theta0, c("shape","rate"), transform_gamma, x)
    } else if (dist == "Beta") {
      # method of moments start
      m <- mean(x); v <- var(x)
      tmp <- max(1e-6, (m*(1-m)/v - 1))
      a0 <- max(1e-3, m * tmp)
      b0 <- max(1e-3, (1-m) * tmp)
      theta0 <- c(log(a0), log(b0))
      res <- mle_optim(nll_beta, theta0, c("alpha","beta"), transform_beta, x)
    } else {
      res <- NULL
    }
    res$dist <- dist
    res$x <- x
    res$size <- if (!is.null(size)) size else NA
    res
  })
  
  output$fitSummary <- renderPrint({
    res <- fit_res()
    if (is.null(res)) return("No fit yet.")
    if (!res$success) {
      cat("Fit failed:", res$message, "\n")
      return()
    }
    cat("Distribution:", res$dist, "\n")
    cat("Log-likelihood:", formatC(res$logLik, digits = 6), "\n")
    cat("Convergence code (0 = success):", res$convergence, "\n\n")
    params <- res$par
    ses <- res$se
    df <- data.frame(Estimate = params, StdErr = ses, row.names = names(params))
    print(df, digits = 6)
    cat("\nNotes: Standard errors are approximate (delta method using optim() Hessian on transformed scale).\n")
  })
  
  output$fitPlot <- renderPlot({
    res <- fit_res()
    if (is.null(res) || !res$success) {
      # show data only if exists
      x <- data_r()
      if (is.null(x)) return(NULL)
      if (is.list(x) && !is.null(x$success)) {
        plot_df <- data.frame(x = x$success)
      } else plot_df <- data.frame(x = as.numeric(x))
      ggplot(plot_df, aes(x = x)) + geom_histogram(bins = 30, fill = "grey80", color = "black") +
        ggtitle("Data (no fit yet)")
    } else {
      x <- res$x
      dist <- res$dist
      par <- res$par
      # plotting for discrete vs continuous distributions
      if (dist %in% c("Poisson", "Binomial")) {
        # discrete: bar chart of observed frequencies + fitted pmf
        xvec <- as.integer(x)
        obs_tab <- as.data.frame(table(xvec))
        names(obs_tab) <- c("x", "count")
        obs_tab$x <- as.integer(as.character(obs_tab$x))
        obs_tab$prop <- obs_tab$count / sum(obs_tab$count)
        # fitted pmf over range
        xmin <- min(obs_tab$x)
        xmax <- max(obs_tab$x)
        xs <- seq(xmin, xmax)
        if (dist == "Poisson") {
          pmf <- dpois(xs, lambda = par["lambda"])
          fit_df <- data.frame(x = xs, pmf = pmf)
        } else { # Binomial
          size <- if (!is.na(res$size[1])) res$size[1] else input$binom_size
          pmf <- dbinom(xs, size = size, prob = par["p"])
          fit_df <- data.frame(x = xs, pmf = pmf)
        }
        ggplot() +
          geom_bar(data = obs_tab, aes(x = x, y = prop), stat = "identity", fill = "grey80", color = "black") +
          geom_point(data = fit_df, aes(x = x, y = pmf), color = "red", size = 2) +
          geom_line(data = fit_df, aes(x = x, y = pmf), color = "red") +
          ylab("Proportion / PMF") + xlab("x") +
          ggtitle(paste0("Observed frequencies and fitted PMF (", dist, ")"))
      } else {
        # continuous: histogram + fitted density curve
        plot_df <- data.frame(x = as.numeric(x))
        gg <- ggplot(plot_df, aes(x = x)) +
          geom_histogram(aes(y = ..density..), bins = 30, fill = "grey90", color = "black")
        xs <- seq(min(plot_df$x), max(plot_df$x), length.out = 400)
        dens <- switch(dist,
                       "Normal" = dnorm(xs, mean = par["mu"], sd = par["sigma"]),
                       "Exponential" = dexp(xs, rate = par["rate"]),
                       "Gamma" = dgamma(xs, shape = par["shape"], rate = par["rate"]),
                       "Beta" = dbeta(xs, shape1 = par["alpha"], shape2 = par["beta"]),
                       rep(0, length(xs))
        )
        fit_df <- data.frame(x = xs, y = dens)
        gg + geom_line(data = fit_df, aes(x = x, y = y), size = 1, color = "red") +
          ggtitle(paste0("Histogram and fitted density (", dist, ")"))
      }
    }
  })
}

# Run the app
shinyApp(ui, server)
