# Tests for makicoint package

test_that("makicoint returns correct structure", {
  set.seed(123)
  n <- 100
  x <- cumsum(rnorm(n))
  y <- 1 + 2*x + rnorm(n)

  result <- makicoint(y, x, model = 0, max_breaks = 1)

  expect_s3_class(result, "makicoint")
  expect_named(result, c("statistic", "critical_values", "breakpoints",
                         "break_fractions", "nobs", "model", "model_name",
                         "max_breaks", "n_regressors", "lags_used",
                         "trimming", "lag_method", "reject"))

  expect_equal(result$nobs, n)
  expect_equal(result$model, 0)
  expect_equal(result$max_breaks, 1)
  expect_equal(result$n_regressors, 1)
})


test_that("makicoint detects structural break", {
  set.seed(42)
  n <- 200
  x <- cumsum(rnorm(n))

  # Create clear structural break at t=100
  y <- 1 + 2*x + 5*(1:n > 100) + rnorm(n, sd = 0.5)

  result <- makicoint(y, x, model = 0, max_breaks = 1, trimming = 0.10)

  # Break should be detected near observation 100
  expect_true(length(result$breakpoints) >= 1)
  expect_true(abs(result$breakpoints[1] - 100) < 20)
})


test_that("makicoint handles multiple regressors", {
  set.seed(456)
  n <- 150
  x1 <- cumsum(rnorm(n))
  x2 <- cumsum(rnorm(n))
  y <- 1 + 2*x1 + 3*x2 + rnorm(n)

  result <- makicoint(y, cbind(x1, x2), model = 2, max_breaks = 2)

  expect_equal(result$n_regressors, 2)
  expect_equal(result$model, 2)
  expect_equal(result$max_breaks, 2)
})


test_that("makicoint validates inputs correctly", {
  set.seed(789)
  n <- 50
  x <- cumsum(rnorm(n))
  y <- 1 + x + rnorm(n)

  # Invalid model
  expect_error(makicoint(y, x, model = 5), "model")

  # Invalid max_breaks
  expect_error(makicoint(y, x, max_breaks = 0), "max_breaks")
  expect_error(makicoint(y, x, max_breaks = 6), "max_breaks")

  # Invalid trimming
  expect_error(makicoint(y, x, trimming = 0), "trimming")
  expect_error(makicoint(y, x, trimming = 0.6), "trimming")

  # Too many regressors
  x_many <- matrix(rnorm(n * 5), ncol = 5)
  expect_error(makicoint(y, x_many), "4 regressors")

  # Sample size too small
  expect_error(makicoint(y[1:20], x[1:20]), "30 observations")

  # Mismatched lengths
  expect_error(makicoint(y, x[1:40]), "same number of observations")
})


test_that("critical values are correct for model 0, k=1", {
  cv <- makicoint:::get_critical_values(k = 1, m = 1, model = 0)
  expect_equal(cv, c(-5.709, -4.602, -4.354), tolerance = 0.001)

  cv <- makicoint:::get_critical_values(k = 1, m = 3, model = 0)
  expect_equal(cv, c(-5.563, -5.083, -4.784), tolerance = 0.001)
})


test_that("critical values are correct for model 2, k=2", {
  cv <- makicoint:::get_critical_values(k = 2, m = 1, model = 2)
  expect_equal(cv, c(-6.020, -5.558, -5.287), tolerance = 0.001)

  cv <- makicoint:::get_critical_values(k = 2, m = 5, model = 2)
  expect_equal(cv, c(-7.839, -7.288, -6.976), tolerance = 0.001)
})


test_that("critical values are correct for model 3, k=4", {
  cv <- makicoint:::get_critical_values(k = 4, m = 1, model = 3)
  expect_equal(cv, c(-7.400, -6.911, -6.649), tolerance = 0.001)

  cv <- makicoint:::get_critical_values(k = 4, m = 5, model = 3)
  expect_equal(cv, c(-10.08, -9.482, -9.151), tolerance = 0.01)
})


test_that("all models run without error", {
  set.seed(101)
  n <- 100
  x <- cumsum(rnorm(n))
  y <- 1 + x + rnorm(n)

  for (m in 0:3) {
    result <- makicoint(y, x, model = m, max_breaks = 1)
    expect_s3_class(result, "makicoint")
    expect_equal(result$model, m)
  }
})


test_that("all lag methods work", {
  set.seed(202)
  n <- 100
  x <- cumsum(rnorm(n))
  y <- 1 + x + rnorm(n)

  for (method in c("t-sig", "fixed", "aic", "bic")) {
    result <- makicoint(y, x, model = 0, max_breaks = 1, lag_method = method)
    expect_s3_class(result, "makicoint")
  }
})


test_that("print method works", {
  set.seed(303)
  n <- 100
  x <- cumsum(rnorm(n))
  y <- 1 + x + rnorm(n)

  result <- makicoint(y, x, model = 0, max_breaks = 1)

  expect_output(print(result), "Maki")
  expect_output(print(result), "Test Statistic")
  expect_output(print(result), "Critical Values")
})


test_that("multiple breaks are detected", {
  set.seed(404)
  n <- 300
  x <- cumsum(rnorm(n))

  # Two structural breaks
  y <- 1 + 2*x + 3*(1:n > 100) + 4*(1:n > 200) + rnorm(n, sd = 0.5)

  result <- makicoint(y, x, model = 0, max_breaks = 2, trimming = 0.10)

  expect_equal(result$max_breaks, 2)
  expect_true(length(result$breakpoints) == 2)
})


test_that("data.frame input works", {
  set.seed(505)
  n <- 100
  df <- data.frame(x1 = cumsum(rnorm(n)), x2 = cumsum(rnorm(n)))
  y <- 1 + df$x1 + df$x2 + rnorm(n)

  result <- makicoint(y, df, model = 2, max_breaks = 1)

  expect_s3_class(result, "makicoint")
  expect_equal(result$n_regressors, 2)
})
