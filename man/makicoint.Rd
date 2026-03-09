\name{makicoint}
\alias{makicoint}
\title{Maki Cointegration Test with Multiple Structural Breaks}
\usage{
makicoint(
  y,
  x,
  model = 2,
  max_breaks = 1,
  trimming = 0.1,
  max_lags = 12,
  lag_method = "t-sig"
)
}
\arguments{
\item{y}{Numeric vector. The dependent variable (must be I(1)).}

\item{x}{Numeric matrix or vector. The independent variable(s) (must be I(1)).
Maximum 4 regressors supported due to critical value availability.}

\item{model}{Integer (0-3). The model specification:
\itemize{
  \item 0: Level shift - breaks only in intercept
  \item 1: Level shift with trend - breaks in intercept, linear trend included
  \item 2: Regime shift - breaks in intercept and slope coefficients
  \item 3: Regime shift with trend - breaks in intercept, slope, and trend
}}

\item{max_breaks}{Integer (1-5). Maximum number of structural breaks to consider.}

\item{trimming}{Numeric (0-0.5). Trimming parameter defining the minimum
segment length as a fraction of sample size. Default is 0.10.}

\item{max_lags}{Integer. Maximum number of lags for ADF test. Default is 12.}

\item{lag_method}{Character. Lag selection method:
\itemize{
  \item "t-sig": General-to-specific based on t-statistic significance (default)
  \item "fixed": Use \code{max_lags} for all tests
  \item "aic": Minimize Akaike Information Criterion
  \item "bic": Minimize Bayesian Information Criterion
}}
}
\value{
An object of class "makicoint" containing:
\item{statistic}{The test statistic (minimum ADF t-statistic)}
\item{critical_values}{Named vector with 1\%, 5\%, and 10\% critical values}
\item{breakpoints}{Vector of estimated break point locations (observation indices)}
\item{break_fractions}{Vector of break point locations as fractions of sample}
\item{nobs}{Number of observations}
\item{model}{Model specification used}
\item{max_breaks}{Maximum breaks tested}
\item{lags_used}{Number of lags used in ADF test}
\item{trimming}{Trimming parameter}
\item{lag_method}{Lag selection method}
\item{reject}{Logical indicating whether null is rejected at 5\% level}
}
\description{
Performs the Maki (2012) cointegration test that allows for an
unknown number of structural breaks in the cointegrating relationship.
The test extends the Gregory-Hansen (1996) framework to multiple breaks.
}
\details{
The Maki (2012) test is based on residual-based cointegration testing with
structural breaks. The procedure:
\enumerate{
  \item Estimates the cointegrating regression with break dummies for each
    candidate break point
  \item Computes the ADF test statistic on the residuals
  \item Selects the break points that minimize the ADF statistic
  \item The minimum ADF statistic is compared to critical values
}

The null hypothesis is no cointegration. Rejection indicates evidence of
cointegration with structural break(s). The test is left-tailed; more
negative values provide stronger evidence against the null.
}
\examples{
# Generate cointegrated series with a structural break
set.seed(123)
n <- 200
x <- cumsum(rnorm(n))

# Create break in intercept at t=100
y <- 1 + 2*x + 3*(1:n > 100) + rnorm(n)

# Test for cointegration with breaks
result <- makicoint(y, x, model = 0, max_breaks = 2)
print(result)

}
\references{
Maki, D. (2012). Tests for cointegration allowing for an unknown number of
breaks. \emph{Economic Modelling}, 29(5), 2011-2015.
\doi{10.1016/j.econmod.2012.04.022}

Gregory, A. W., & Hansen, B. E. (1996). Residual-based tests for
cointegration in models with regime shifts. \emph{Journal of Econometrics},
70(1), 99-126. \doi{10.1016/0304-4076(69)41685-7}
}
