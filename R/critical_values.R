#' Get Critical Values for Maki (2012) Test
#'
#' @description Returns critical values from Maki (2012) Table 1 for the
#'   specified number of regressors, breaks, and model type.
#'
#' @param k Integer. Number of independent variables (1-4).
#' @param m Integer. Number of breaks (1-5).
#' @param model Integer. Model specification (0-3).
#'
#' @return Numeric vector of length 3 with 1%, 5%, and 10% critical values.
#'
#' @details
#' Critical values are from Table 1 of Maki (2012), which were obtained
#' via Monte Carlo simulation with 50,000 replications.
#'
#' @references
#' Maki, D. (2012). Tests for cointegration allowing for an unknown number of
#' breaks. \emph{Economic Modelling}, 29(5), 2011-2015.
#' \doi{10.1016/j.econmod.2012.04.022}
#'
#' @keywords internal
#' @noRd
get_critical_values <- function(k, m, model) {
  # Critical values from Maki (2012) Table 1
  # Dimensions: [k][model][m, cv_level]
  # cv_level: 1=1%, 2=5%, 3=10%

  # Model 0: Level shift
  cv_model0 <- list(
    # k=1
    matrix(c(
      -5.709, -4.602, -4.354,  # m=1
      -5.416, -4.892, -4.610,  # m=2
      -5.563, -5.083, -4.784,  # m=3
      -5.776, -5.230, -4.982,  # m=4
      -5.959, -5.426, -5.131   # m=5
    ), nrow = 5, byrow = TRUE),
    # k=2
    matrix(c(
      -5.541, -5.004, -4.733,
      -5.717, -5.211, -4.957,
      -5.943, -5.392, -5.125,
      -6.075, -5.550, -5.297,
      -6.296, -5.760, -5.491
    ), nrow = 5, byrow = TRUE),
    # k=3
    matrix(c(
      -5.820, -5.341, -5.101,
      -5.984, -5.517, -5.272,
      -6.229, -5.704, -5.427,
      -6.406, -5.871, -5.603,
      -6.555, -6.038, -5.773
    ), nrow = 5, byrow = TRUE),
    # k=4
    matrix(c(
      -6.139, -5.650, -5.386,
      -6.303, -5.839, -5.575,
      -6.501, -5.992, -5.714,
      -6.640, -6.132, -5.892,
      -6.856, -6.306, -6.039
    ), nrow = 5, byrow = TRUE)
  )

  # Model 1: Level shift with trend
  cv_model1 <- list(
    # k=1
    matrix(c(
      -5.524, -5.038, -4.784,
      -5.708, -5.196, -4.938,
      -5.833, -5.373, -5.106,
      -6.059, -5.508, -5.245,
      -6.193, -5.699, -5.449
    ), nrow = 5, byrow = TRUE),
    # k=2
    matrix(c(
      -5.840, -5.359, -5.117,
      -6.011, -5.518, -5.247,
      -6.169, -5.691, -5.408,
      -6.329, -5.831, -5.558,
      -6.530, -5.993, -5.722
    ), nrow = 5, byrow = TRUE),
    # k=3
    matrix(c(
      -6.144, -5.645, -5.398,
      -6.271, -5.796, -5.538,
      -6.472, -5.957, -5.682,
      -6.575, -6.086, -5.820,
      -6.784, -6.250, -5.976
    ), nrow = 5, byrow = TRUE),
    # k=4
    matrix(c(
      -6.361, -5.913, -5.686,
      -6.556, -6.055, -5.805,
      -6.741, -6.214, -5.974,
      -6.845, -6.373, -6.096,
      -7.053, -6.494, -6.220
    ), nrow = 5, byrow = TRUE)
  )

  # Model 2: Regime shift
  cv_model2 <- list(
    # k=1
    matrix(c(
      -5.457, -4.895, -4.626,
      -5.863, -5.363, -5.070,
      -6.251, -5.703, -5.402,
      -6.596, -6.011, -5.723,
      -6.915, -6.357, -6.057
    ), nrow = 5, byrow = TRUE),
    # k=2
    matrix(c(
      -6.020, -5.558, -5.287,
      -6.628, -6.093, -5.833,
      -7.031, -6.516, -6.210,
      -7.470, -6.872, -6.563,
      -7.839, -7.288, -6.976
    ), nrow = 5, byrow = TRUE),
    # k=3
    matrix(c(
      -6.565, -6.035, -5.773,
      -7.232, -6.702, -6.411,
      -7.767, -7.155, -6.868,
      -8.236, -7.625, -7.329,
      -8.673, -8.110, -7.796
    ), nrow = 5, byrow = TRUE),
    # k=4
    matrix(c(
      -7.021, -6.520, -6.242,
      -7.756, -7.244, -6.964,
      -8.336, -7.803, -7.481,
      -8.895, -8.292, -8.004,
      -9.441, -8.869, -8.541
    ), nrow = 5, byrow = TRUE)
  )

  # Model 3: Regime shift with trend
  cv_model3 <- list(
    # k=1
    matrix(c(
      -6.048, -5.541, -5.281,
      -6.620, -6.100, -5.845,
      -7.082, -6.524, -6.267,
      -7.553, -7.009, -6.712,
      -8.004, -7.414, -7.110
    ), nrow = 5, byrow = TRUE),
    # k=2
    matrix(c(
      -6.523, -6.055, -5.795,
      -7.153, -6.657, -6.397,
      -7.673, -7.145, -6.873,
      -8.217, -7.636, -7.341,
      -8.713, -8.129, -7.811
    ), nrow = 5, byrow = TRUE),
    # k=3
    matrix(c(
      -6.964, -6.464, -6.220,
      -7.737, -7.201, -6.926,
      -8.331, -7.743, -7.449,
      -8.851, -8.269, -7.960,
      -9.428, -8.800, -8.508
    ), nrow = 5, byrow = TRUE),
    # k=4
    matrix(c(
      -7.400, -6.911, -6.649,
      -8.167, -7.638, -7.381,
      -8.865, -8.254, -7.977,
      -9.433, -8.871, -8.574,
      -10.08, -9.482, -9.151
    ), nrow = 5, byrow = TRUE)
  )

  # Select appropriate critical values
  cv_all <- list(cv_model0, cv_model1, cv_model2, cv_model3)

  cv <- cv_all[[model + 1]][[k]][m, ]

  return(cv)
}
