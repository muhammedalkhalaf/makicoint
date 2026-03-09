\name{makicoint-package}
\alias{makicoint-package}
\docType{package}
\title{
\packageTitle{makicoint}
}
\description{
The makicoint package implements the Maki (2012) cointegration test
that allows for an unknown number of structural breaks in the
cointegrating relationship. This test extends the Gregory-Hansen (1996)
framework to handle multiple breaks, making it particularly useful for
analyzing long economic time series where multiple regime changes may
have occurred.
}
\details{
The package supports four model specifications:
\itemize{
  \item Model 0: Level shift - breaks only in the intercept
  \item Model 1: Level shift with trend - breaks in intercept with linear trend
  \item Model 2: Regime shift - breaks in both intercept and slope coefficients
  \item Model 3: Regime shift with trend - breaks in intercept, slope, and trend
}

The test procedure:
\enumerate{
  \item Estimates the cointegrating regression with break dummies
  \item Tests residuals for unit root using ADF test
  \item Searches over all possible break locations
  \item Returns minimum ADF statistic across all break configurations
}

The null hypothesis is no cointegration. The test is left-tailed; more
negative values of the test statistic provide stronger evidence against
the null (i.e., evidence for cointegration with structural breaks).
}
\author{
\packageAuthor{makicoint}

Maintainer: \packageMaintainer{makicoint}
}
\references{
Maki, D. (2012). Tests for cointegration allowing for an unknown number of
breaks. \emph{Economic Modelling}, 29(5), 2011-2015.
\doi{10.1016/j.econmod.2012.04.022}

Gregory, A. W., & Hansen, B. E. (1996). Residual-based tests for
cointegration in models with regime shifts. \emph{Journal of Econometrics},
70(1), 99-126. \doi{10.1016/0304-4076(69)41685-7}
}
\keyword{package}
\seealso{
\code{\link{makicoint}}
}
