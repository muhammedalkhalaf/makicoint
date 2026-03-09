#' @title Maki Cointegration Test with Multiple Structural Breaks
#'
#' @description Performs the Maki (2012) cointegration test that allows for an
#'   unknown number of structural breaks in the cointegrating relationship.
#'   The test extends the Gregory-Hansen (1996) framework to multiple breaks.
#'
#' @param y Numeric vector. The dependent variable (must be I(1)).
#' @param x Numeric matrix or vector. The independent variable(s) (must be I(1)).
#'   Maximum 4 regressors supported due to critical value availability.
#' @param model Integer (0-3). The model specification:
#'   \itemize{
#'     \item 0: Level shift - breaks only in intercept
#'     \item 1: Level shift with trend - breaks in intercept, linear trend included
#'     \item 2: Regime shift - breaks in intercept and slope coefficients
#'     \item 3: Regime shift with trend - breaks in intercept, slope, and trend
#'   }
#' @param max_breaks Integer (1-5). Maximum number of structural breaks to consider.
#' @param trimming Numeric (0-0.5). Trimming parameter defining the minimum
#'   segment length as a fraction of sample size. Default is 0.10.
#' @param max_lags Integer. Maximum number of lags for ADF test. Default is 12.
#' @param lag_method Character. Lag selection method:
#'   \itemize{
#'     \item "t-sig": General-to-specific based on t-statistic significance (default)
#'     \item "fixed": Use \code{max_lags} for all tests
#'     \item "aic": Minimize Akaike Information Criterion
#'     \item "bic": Minimize Bayesian Information Criterion
#'   }
#'
#' @return An object of class "makicoint" containing:
#'   \item{statistic}{The test statistic (minimum ADF t-statistic)}
#'   \item{critical_values}{Named vector with 1%, 5%, and 10% critical values}
#'   \item{breakpoints}{Vector of estimated break point locations (observation indices)}
#'   \item{break_fractions}{Vector of break point locations as fractions of sample}
#'   \item{nobs}{Number of observations}
#'   \item{model}{Model specification used}
#'   \item{max_breaks}{Maximum breaks tested}
#'   \item{lags_used}{Number of lags used in ADF test}
#'   \item{trimming}{Trimming parameter}
#'   \item{lag_method}{Lag selection method}
#'   \item{reject}{Logical indicating whether null is rejected at 5% level}
#'
#' @details
#' The Maki (2012) test is based on residual-based cointegration testing with

#' structural breaks. The procedure:
#' \enumerate{
#'   \item Estimates the cointegrating regression with break dummies for each
#'     candidate break point
#'   \item Computes the ADF test statistic on the residuals
#'   \item Selects the break points that minimize the ADF statistic
#'   \item The minimum ADF statistic is compared to critical values
#' }
#'
#' The null hypothesis is no cointegration. Rejection indicates evidence of
#' cointegration with structural break(s). The test is left-tailed; more
#' negative values provide stronger evidence against the null.
#'
#' @references
#' Maki, D. (2012). Tests for cointegration allowing for an unknown number of
#' breaks. \emph{Economic Modelling}, 29(5), 2011-2015.
#' \doi{10.1016/j.econmod.2012.04.022}
#'
#' Gregory, A. W., & Hansen, B. E. (1996). Residual-based tests for
#' cointegration in models with regime shifts. \emph{Journal of Econometrics},
#' 70(1), 99-126. \doi{10.1016/0304-4076(69)41685-7}
#'
#' @examples
#' # Generate cointegrated series with a structural break
#' set.seed(123)
#' n <- 200
#' e <- rnorm(n)
#' x <- cumsum(rnorm(n))
#'
#' # Create break in intercept at t=100
#' y <- 1 + 2*x + 3*(1:n > 100) + e
#'
#' # Test for cointegration with breaks
#' result <- makicoint(y, x, model = 0, max_breaks = 2)
#' print(result)
#' summary(result)
#'
#' @export
makicoint <- function(y, x, model = 2, max_breaks = 1, trimming = 0.10,
                      max_lags = 12, lag_method = "t-sig") {

  # Input validation
  if (!is.numeric(y) || !is.vector(y)) {
    stop("'y' must be a numeric vector.")
  }

  if (is.vector(x)) {
    x <- matrix(x, ncol = 1)
  } else if (is.data.frame(x)) {
    x <- as.matrix(x)
  } else if (!is.matrix(x)) {
    stop("'x' must be a numeric vector, matrix, or data frame.")
  }

  if (!is.numeric(x)) {
    stop("'x' must contain numeric values.")
  }

  n <- length(y)
  k <- ncol(x)

  if (nrow(x) != n) {
    stop("'y' and 'x' must have the same number of observations.")
  }

  if (any(is.na(y)) || any(is.na(x))) {
    stop("Missing values are not allowed. Please remove or impute NA values.")
  }

  if (n < 30) {
    stop("Sample size too small. Minimum 30 observations required.")
  }

  if (!model %in% 0:3) {
    stop("'model' must be 0 (level shift), 1 (level shift with trend), ",
         "2 (regime shift), or 3 (regime shift with trend).")
  }

  if (!max_breaks %in% 1:5) {
    stop("'max_breaks' must be between 1 and 5.")
  }

  if (k < 1 || k > 4) {
    stop("Number of independent variables must be between 1 and 4. ",
         "Critical values are not available for more than 4 regressors.")
  }

  if (trimming <= 0 || trimming >= 0.5) {
    stop("'trimming' must be between 0 and 0.5 (exclusive).")
  }

  if (max_lags < 0) {
    stop("'max_lags' must be non-negative.")
  }

  lag_method <- tolower(lag_method)
  lag_method <- gsub("-", "", lag_method)  # Handle "t-sig" vs "tsig"

  if (!lag_method %in% c("tsig", "fixed", "aic", "bic")) {
    stop("'lag_method' must be one of: 't-sig', 'fixed', 'aic', 'bic'.")
  }

  # Minimum recommended sample size

min_recommended <- ceiling((max_breaks + 1) / trimming)
  if (n < min_recommended) {
    warning(sprintf(
      "Sample size (%d) is small for %d breaks with trimming %.2f. ",
      n, max_breaks, trimming),
      sprintf("Recommended minimum: %d observations. ", min_recommended),
      "Results may be unreliable. Consider reducing max_breaks or trimming."
    )
  }

  # Combine data
  datap <- cbind(y, x)

  # Trimming in observations
  tb <- round(trimming * n)

  # Run the test
  result <- maki_test(datap, max_breaks, model, tb, max_lags, lag_method)

  # Get critical values
  cv <- get_critical_values(k, max_breaks, model)

  # Create return object
  out <- list(
    statistic = result$test_stat,
    critical_values = c("1%" = cv[1], "5%" = cv[2], "10%" = cv[3]),
    breakpoints = result$breakpoints[result$breakpoints > 0],
    break_fractions = result$breakpoints[result$breakpoints > 0] / n,
    nobs = n,
    model = model,
    model_name = c("Level Shift", "Level Shift with Trend",
                   "Regime Shift", "Regime Shift with Trend")[model + 1],
    max_breaks = max_breaks,
    n_regressors = k,
    lags_used = result$lags_used,
    trimming = trimming,
    lag_method = lag_method,
    reject = result$test_stat < cv[2]  # Reject at 5% level
  )

  class(out) <- "makicoint"
  return(out)
}

#' @export
print.makicoint <- function(x, ...) {
  cat("\n")
  cat("Maki (2012) Cointegration Test with Multiple Structural Breaks\n")
  cat(strrep("-", 65), "\n")

  cat(sprintf("Model: %s\n", x$model_name))
  cat(sprintf("Observations: %d | Regressors: %d | Max breaks: %d\n",
              x$nobs, x$n_regressors, x$max_breaks))
  cat(sprintf("Trimming: %.2f | Lags used: %d\n", x$trimming, x$lags_used))
  cat("\n")

  cat("H0: No cointegration\n")
  cat(sprintf("H1: Cointegration with up to %d break(s)\n", x$max_breaks))
  cat("\n")

  cat(sprintf("Test Statistic: %.4f\n", x$statistic))
  cat(sprintf("Critical Values:  1%%: %.3f  5%%: %.3f  10%%: %.3f\n",
              x$critical_values[1], x$critical_values[2], x$critical_values[3]))
  cat("\n")

  if (length(x$breakpoints) > 0) {
    cat("Estimated Break Points:\n")
    for (i in seq_along(x$breakpoints)) {
      cat(sprintf("  Break %d: Observation %d (fraction: %.4f)\n",
                  i, x$breakpoints[i], x$break_fractions[i]))
    }
    cat("\n")
  }

  if (x$statistic < x$critical_values[1]) {
    cat("Conclusion: Reject H0 at 1% significance level.\n")
    cat("            Evidence of cointegration with structural break(s).\n")
  } else if (x$statistic < x$critical_values[2]) {
    cat("Conclusion: Reject H0 at 5% significance level.\n")
    cat("            Evidence of cointegration with structural break(s).\n")
  } else if (x$statistic < x$critical_values[3]) {
    cat("Conclusion: Reject H0 at 10% significance level.\n")
    cat("            Evidence of cointegration with structural break(s).\n")
  } else {
    cat("Conclusion: Fail to reject H0.\n")
    cat("            No evidence of cointegration.\n")
  }
  cat(strrep("-", 65), "\n")

  invisible(x)
}

#' @export
summary.makicoint <- function(object, ...) {
  print(object, ...)
}


# ============================================================================
# Internal Functions
# ============================================================================

#' Main Maki test procedure
#' @noRd
maki_test <- function(datap, m, model, tb, max_lags, lag_method) {
  n <- nrow(datap)
  breakpoints <- rep(0, 5)
  lags_used <- 0

  if (m == 1) {
    result1 <- mbreak1(datap, n, model, tb, max_lags, lag_method)
    mintau <- result1$mintau
    breakpoints[1] <- result1$bp
    lags_used <- result1$lags_used
  }
  else if (m == 2) {
    result1 <- mbreak1(datap, n, model, tb, max_lags, lag_method)
    result2 <- mbreak2(datap, n, model, tb, result1$bp, max_lags, lag_method)

    alltau <- c(result1$mintau, result2$mintau)
    mintau <- min(alltau)

    bp <- sort(c(result1$bp, result2$bp))
    breakpoints[1:2] <- bp
    lags_used <- result2$lags_used
  }
  else if (m == 3) {
    result1 <- mbreak1(datap, n, model, tb, max_lags, lag_method)
    result2 <- mbreak2(datap, n, model, tb, result1$bp, max_lags, lag_method)

    bp <- sort(c(result1$bp, result2$bp))
    result3 <- mbreak3(datap, n, model, tb, bp, max_lags, lag_method)

    alltau <- c(result1$mintau, result2$mintau, result3$mintau)
    mintau <- min(alltau)

    bp <- sort(c(result1$bp, result2$bp, result3$bp))
    breakpoints[1:3] <- bp
    lags_used <- result3$lags_used
  }
  else if (m == 4) {
    result1 <- mbreak1(datap, n, model, tb, max_lags, lag_method)
    result2 <- mbreak2(datap, n, model, tb, result1$bp, max_lags, lag_method)

    bp <- sort(c(result1$bp, result2$bp))
    result3 <- mbreak3(datap, n, model, tb, bp, max_lags, lag_method)

    bp123 <- sort(c(result1$bp, result2$bp, result3$bp))
    result4 <- mbreak4(datap, n, model, tb, bp123, max_lags, lag_method)

    alltau <- c(result1$mintau, result2$mintau, result3$mintau, result4$mintau)
    mintau <- min(alltau)

    bp <- sort(c(result1$bp, result2$bp, result3$bp, result4$bp))
    breakpoints[1:4] <- bp
    lags_used <- result4$lags_used
  }
  else if (m == 5) {
    result1 <- mbreak1(datap, n, model, tb, max_lags, lag_method)
    result2 <- mbreak2(datap, n, model, tb, result1$bp, max_lags, lag_method)

    bp <- sort(c(result1$bp, result2$bp))
    result3 <- mbreak3(datap, n, model, tb, bp, max_lags, lag_method)

    bp123 <- sort(c(result1$bp, result2$bp, result3$bp))
    result4 <- mbreak4(datap, n, model, tb, bp123, max_lags, lag_method)

    bp1234 <- sort(c(result1$bp, result2$bp, result3$bp, result4$bp))
    result5 <- mbreak5(datap, n, model, tb, bp1234, max_lags, lag_method)

    alltau <- c(result1$mintau, result2$mintau, result3$mintau,
                result4$mintau, result5$mintau)
    mintau <- min(alltau)

    bp <- sort(c(result1$bp, result2$bp, result3$bp, result4$bp, result5$bp))
    breakpoints[1:5] <- bp
    lags_used <- result5$lags_used
  }

  list(test_stat = mintau, breakpoints = breakpoints, lags_used = lags_used)
}


#' Search for first break point
#' @noRd
mbreak1 <- function(datap, n, model, tb, max_lags, lag_method) {
  y <- datap[, 1]
  k <- ncol(datap)

  vectau <- rep(NA, n)
  vecssr <- rep(NA, n)
  lags_used <- 0

  for (i in (tb + 1):(n - tb)) {
    u <- rep(1, n)
    du <- c(rep(0, i), rep(1, n - i))

    X <- build_design_matrix(datap, model, n, k, list(i), u)

    e <- get_residuals(y, X)
    result <- adf_test(e, max_lags, lag_method)
    tau <- result$tau
    lags_used <- result$lags_used

    vectau[i] <- tau
    vecssr[i] <- sum(e^2)
  }

  mintau <- min(vectau[(tb + 1):(n - tb)], na.rm = TRUE)
  minidx <- which.min(vecssr[(tb + 1):(n - tb)]) + tb

  list(mintau = mintau, bp = minidx, lags_used = lags_used)
}


#' Search for second break point
#' @noRd
mbreak2 <- function(datap, n, model, tb, bp1_in, max_lags, lag_method) {
  if (bp1_in <= 0.1 * n) {
    result <- mbreak22(datap, n, model, tb, bp1_in, max_lags, lag_method)
  } else if (bp1_in >= 0.9 * n) {
    result <- mbreak21(datap, n, model, tb, bp1_in, max_lags, lag_method)
  } else {
    result1 <- mbreak21(datap, n, model, tb, bp1_in, max_lags, lag_method)
    result2 <- mbreak22(datap, n, model, tb, bp1_in, max_lags, lag_method)

    if (result1$mintau < result2$mintau) {
      result <- result1
    } else {
      result <- result2
    }
  }

  result
}


#' Search for second break point before first break
#' @noRd
mbreak21 <- function(datap, n, model, tb, bp1_in, max_lags, lag_method) {
  y <- datap[, 1]
  k <- ncol(datap)

  vectau <- rep(NA, n)
  vecssr <- rep(NA, n)
  lags_used <- 0

  if (bp1_in - tb < tb + 1) {
    return(list(mintau = Inf, bp = tb + 1, lags_used = 0))
  }

  for (i in (tb + 1):(bp1_in - tb)) {
    X <- build_design_matrix(datap, model, n, k, list(bp1_in, i), rep(1, n))

    e <- get_residuals(y, X)
    result <- adf_test(e, max_lags, lag_method)
    tau <- result$tau
    lags_used <- result$lags_used

    vectau[i] <- tau
    vecssr[i] <- sum(e^2)
  }

  valid_range <- (tb + 1):(bp1_in - tb)
  mintau <- min(vectau[valid_range], na.rm = TRUE)
  minidx <- valid_range[which.min(vecssr[valid_range])]

  list(mintau = mintau, bp = minidx, lags_used = lags_used)
}


#' Search for second break point after first break
#' @noRd
mbreak22 <- function(datap, n, model, tb, bp1_in, max_lags, lag_method) {
  y <- datap[, 1]
  k <- ncol(datap)

  vectau <- rep(NA, n)
  vecssr <- rep(NA, n)
  lags_used <- 0

  start_i <- bp1_in + tb + 1
  end_i <- n - tb

  if (start_i > end_i) {
    return(list(mintau = Inf, bp = start_i, lags_used = 0))
  }

  for (i in start_i:end_i) {
    X <- build_design_matrix(datap, model, n, k, list(bp1_in, i), rep(1, n))

    e <- get_residuals(y, X)
    result <- adf_test(e, max_lags, lag_method)
    tau <- result$tau
    lags_used <- result$lags_used

    vectau[i] <- tau
    vecssr[i] <- sum(e^2)
  }

  valid_range <- start_i:end_i
  mintau <- min(vectau[valid_range], na.rm = TRUE)
  minidx <- valid_range[which.min(vecssr[valid_range])]

  list(mintau = mintau, bp = minidx, lags_used = lags_used)
}


#' Search for third break point
#' @noRd
mbreak3 <- function(datap, n, model, tb, bp_in, max_lags, lag_method) {
  bp1 <- bp_in[1]
  bp2 <- bp_in[2]

  results <- list()

  # Region before bp1
  if (bp1 > 2 * tb) {
    results[[length(results) + 1]] <- mbreak_region(
      datap, n, model, tb, bp_in, tb + 1, bp1 - tb, max_lags, lag_method
    )
  }

  # Region between bp1 and bp2
  if (bp2 - bp1 > 2 * tb) {
    results[[length(results) + 1]] <- mbreak_region(
      datap, n, model, tb, bp_in, bp1 + tb + 1, bp2 - tb, max_lags, lag_method
    )
  }

  # Region after bp2
  if (n - bp2 > 2 * tb) {
    results[[length(results) + 1]] <- mbreak_region(
      datap, n, model, tb, bp_in, bp2 + tb + 1, n - tb, max_lags, lag_method
    )
  }

  if (length(results) == 0) {
    return(list(mintau = Inf, bp = round(n / 2), lags_used = 0))
  }

  alltau <- sapply(results, function(r) r$mintau)
  best_idx <- which.min(alltau)

  results[[best_idx]]
}


#' Search for fourth break point
#' @noRd
mbreak4 <- function(datap, n, model, tb, bp_in, max_lags, lag_method) {
  bp1 <- bp_in[1]
  bp2 <- bp_in[2]
  bp3 <- bp_in[3]

  results <- list()

  # Region before bp1
  if (bp1 > 2 * tb) {
    results[[length(results) + 1]] <- mbreak_region(
      datap, n, model, tb, bp_in, tb + 1, bp1 - tb, max_lags, lag_method
    )
  }

  # Region between bp1 and bp2
  if (bp2 - bp1 > 2 * tb) {
    results[[length(results) + 1]] <- mbreak_region(
      datap, n, model, tb, bp_in, bp1 + tb + 1, bp2 - tb, max_lags, lag_method
    )
  }

  # Region between bp2 and bp3
  if (bp3 - bp2 > 2 * tb) {
    results[[length(results) + 1]] <- mbreak_region(
      datap, n, model, tb, bp_in, bp2 + tb + 1, bp3 - tb, max_lags, lag_method
    )
  }

  # Region after bp3
  if (n - bp3 > 2 * tb) {
    results[[length(results) + 1]] <- mbreak_region(
      datap, n, model, tb, bp_in, bp3 + tb + 1, n - tb, max_lags, lag_method
    )
  }

  if (length(results) == 0) {
    return(list(mintau = Inf, bp = round(n / 2), lags_used = 0))
  }

  alltau <- sapply(results, function(r) r$mintau)
  best_idx <- which.min(alltau)

  results[[best_idx]]
}


#' Search for fifth break point
#' @noRd
mbreak5 <- function(datap, n, model, tb, bp_in, max_lags, lag_method) {
  bp1 <- bp_in[1]
  bp2 <- bp_in[2]
  bp3 <- bp_in[3]
  bp4 <- bp_in[4]

  results <- list()

  # Region before bp1
  if (bp1 > 2 * tb) {
    results[[length(results) + 1]] <- mbreak_region(
      datap, n, model, tb, bp_in, tb + 1, bp1 - tb, max_lags, lag_method
    )
  }

  # Region between bp1 and bp2
  if (bp2 - bp1 > 2 * tb) {
    results[[length(results) + 1]] <- mbreak_region(
      datap, n, model, tb, bp_in, bp1 + tb + 1, bp2 - tb, max_lags, lag_method
    )
  }

  # Region between bp2 and bp3
  if (bp3 - bp2 > 2 * tb) {
    results[[length(results) + 1]] <- mbreak_region(
      datap, n, model, tb, bp_in, bp2 + tb + 1, bp3 - tb, max_lags, lag_method
    )
  }

  # Region between bp3 and bp4
  if (bp4 - bp3 > 2 * tb) {
    results[[length(results) + 1]] <- mbreak_region(
      datap, n, model, tb, bp_in, bp3 + tb + 1, bp4 - tb, max_lags, lag_method
    )
  }

  # Region after bp4
  if (n - bp4 > 2 * tb) {
    results[[length(results) + 1]] <- mbreak_region(
      datap, n, model, tb, bp_in, bp4 + tb + 1, n - tb, max_lags, lag_method
    )
  }

  if (length(results) == 0) {
    return(list(mintau = Inf, bp = round(n / 2), lags_used = 0))
  }

  alltau <- sapply(results, function(r) r$mintau)
  best_idx <- which.min(alltau)

  results[[best_idx]]
}


#' Search for break point in a specific region
#' @noRd
mbreak_region <- function(datap, n, model, tb, bp_in, start_i, end_i,
                          max_lags, lag_method) {
  y <- datap[, 1]
  k <- ncol(datap)

  vectau <- rep(NA, n)
  vecssr <- rep(NA, n)
  lags_used <- 0

  if (start_i > end_i) {
    return(list(mintau = Inf, bp = start_i, lags_used = 0))
  }

  for (i in start_i:end_i) {
    all_breaks <- c(bp_in, i)
    X <- build_design_matrix(datap, model, n, k, as.list(all_breaks), rep(1, n))

    e <- get_residuals(y, X)
    result <- adf_test(e, max_lags, lag_method)
    tau <- result$tau
    lags_used <- result$lags_used

    vectau[i] <- tau
    vecssr[i] <- sum(e^2)
  }

  valid_range <- start_i:end_i
  mintau <- min(vectau[valid_range], na.rm = TRUE)
  minidx <- valid_range[which.min(vecssr[valid_range])]

  list(mintau = mintau, bp = minidx, lags_used = lags_used)
}


#' Build the design matrix for cointegrating regression
#' @noRd
build_design_matrix <- function(datap, model, n, k, breaks, u) {
  # Intercept
  X <- matrix(1, nrow = n, ncol = 1)

  # Add break dummies
  for (bp in breaks) {
    if (bp > 0 && bp < n) {
      du <- c(rep(0, bp), rep(1, n - bp))
      X <- cbind(X, du)
    }
  }

  # Model 0: Level shift - intercept + break dummies + regressors
  if (model == 0) {
    X <- cbind(X, datap[, 2:k, drop = FALSE])
  }
  # Model 1: Level shift with trend - add trend
  else if (model == 1) {
    tr <- 1:n
    X <- cbind(X, tr, datap[, 2:k, drop = FALSE])
  }
  # Model 2: Regime shift - add regressors and regime-shifted regressors
  else if (model == 2) {
    X <- cbind(X, datap[, 2:k, drop = FALSE])
    for (bp in breaks) {
      if (bp > 0 && bp < n) {
        dx <- rbind(matrix(0, nrow = bp, ncol = k - 1),
                    datap[(bp + 1):n, 2:k, drop = FALSE])
        X <- cbind(X, dx)
      }
    }
  }
  # Model 3: Regime shift with trend - add trend, trend breaks, and regime shifts
  else if (model == 3) {
    tr <- 1:n
    X <- cbind(X, tr)

    # Add trend breaks
    for (bp in breaks) {
      if (bp > 0 && bp < n) {
        dtr <- c(rep(0, bp), (bp + 1):n)
        X <- cbind(X, dtr)
      }
    }

    # Add regressors
    X <- cbind(X, datap[, 2:k, drop = FALSE])

    # Add regime-shifted regressors
    for (bp in breaks) {
      if (bp > 0 && bp < n) {
        dx <- rbind(matrix(0, nrow = bp, ncol = k - 1),
                    datap[(bp + 1):n, 2:k, drop = FALSE])
        X <- cbind(X, dx)
      }
    }
  }

  X
}


#' Get OLS residuals
#' @noRd
get_residuals <- function(y, X) {
  # Use solve with QR decomposition for numerical stability
  qr_X <- qr(X)
  b <- qr.solve(qr_X, y)
  e <- y - X %*% b
  as.vector(e)
}


#' Perform ADF test on residuals
#' @noRd
adf_test <- function(e, max_lags, lag_method) {
  n <- length(e)

  if (lag_method == "fixed") {
    lag <- max_lags
  } else if (lag_method == "tsig") {
    lag <- optimal_lag_tsig(e, max_lags)
  } else if (lag_method == "aic") {
    lag <- optimal_lag_ic(e, max_lags, "aic")
  } else if (lag_method == "bic") {
    lag <- optimal_lag_ic(e, max_lags, "bic")
  } else {
    lag <- optimal_lag_tsig(e, max_lags)
  }

  tau <- compute_adf_tau(e, lag)

  list(tau = tau, lags_used = lag)
}


#' Select optimal lag using t-sig method
#' @noRd
optimal_lag_tsig <- function(e, max_lags) {
  n <- length(e)
  dy <- diff(e)

  for (p in max_lags:1) {
    if (n - 1 - p < 3) next

    # Build regression matrix
    start_idx <- p + 1
    end_idx <- n - 1

    y_reg <- dy[start_idx:end_idx]
    X_reg <- e[start_idx:end_idx]

    if (p > 0) {
      for (j in 1:p) {
        X_reg <- cbind(X_reg, dy[(start_idx - j):(end_idx - j)])
      }
    }

    # Fit regression
    fit <- tryCatch({
      qr_X <- qr(X_reg)
      b <- qr.solve(qr_X, y_reg)
      resid <- y_reg - X_reg %*% b
      s2 <- sum(resid^2) / (length(y_reg) - ncol(X_reg))
      XtX_inv <- chol2inv(qr.R(qr_X))
      se <- sqrt(diag(s2 * XtX_inv))
      list(b = b, se = se)
    }, error = function(e) NULL)

    if (is.null(fit)) next

    # Check significance of last lag
    tstat <- fit$b[p + 1] / fit$se[p + 1]

    if (abs(tstat) > 1.654) {  # 10% significance level
      return(p)
    }
  }

  return(0)
}


#' Select optimal lag using information criterion
#' @noRd
optimal_lag_ic <- function(e, max_lags, criterion) {
  n <- length(e)
  dy <- diff(e)

  best_lag <- 0
  best_ic <- Inf

  for (p in 0:max_lags) {
    if (n - 1 - p < 3) next

    start_idx <- p + 1
    end_idx <- n - 1

    y_reg <- dy[start_idx:end_idx]
    X_reg <- e[start_idx:end_idx]

    if (p > 0) {
      for (j in 1:p) {
        X_reg <- cbind(X_reg, dy[(start_idx - j):(end_idx - j)])
      }
    }

    # Fit regression
    fit <- tryCatch({
      qr_X <- qr(X_reg)
      b <- qr.solve(qr_X, y_reg)
      resid <- y_reg - X_reg %*% b
      list(resid = resid)
    }, error = function(e) NULL)

    if (is.null(fit)) next

    nobs <- length(y_reg)
    npar <- ncol(as.matrix(X_reg))
    s2 <- sum(fit$resid^2) / nobs

    if (criterion == "aic") {
      ic <- log(s2) + 2 * npar / nobs
    } else {
      ic <- log(s2) + log(nobs) * npar / nobs
    }

    if (ic < best_ic) {
      best_ic <- ic
      best_lag <- p
    }
  }

  best_lag
}


#' Compute ADF t-statistic
#' @noRd
compute_adf_tau <- function(e, lag) {
  n <- length(e)
  dy <- diff(e)

  r <- 2 + lag
  start_idx <- r - 1
  end_idx <- n - 1

  if (start_idx > end_idx || start_idx < 1) {
    return(0)
  }

  y_reg <- dy[start_idx:end_idx]
  X_reg <- e[start_idx:end_idx]

  if (lag > 0) {
    for (q in 1:lag) {
      X_reg <- cbind(X_reg, dy[(start_idx - q):(end_idx - q)])
    }
  }

  # Fit regression and get t-statistic
  fit <- tryCatch({
    qr_X <- qr(X_reg)
    b <- qr.solve(qr_X, y_reg)
    resid <- y_reg - X_reg %*% b
    s2 <- sum(resid^2) / (length(y_reg) - ncol(as.matrix(X_reg)))
    XtX_inv <- chol2inv(qr.R(qr_X))
    se <- sqrt(diag(s2 * XtX_inv))
    list(b = b, se = se)
  }, error = function(e) list(b = 0, se = 1))

  tau <- fit$b[1] / fit$se[1]
  tau
}
