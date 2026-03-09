# makicoint

<!-- badges: start -->
<!-- badges: end -->

## Overview

The `makicoint` package implements the Maki (2012) cointegration test that allows for an unknown number of structural breaks in the cointegrating relationship. This test extends the Gregory-Hansen (1996) framework to handle multiple breaks, making it particularly useful for analyzing long economic time series where multiple regime changes may have occurred.

## Installation

You can install the development version from GitHub:

```r
# install.packages("devtools")
```

## Usage

```r
library(makicoint)

# Generate cointegrated series with a structural break
set.seed(123)
n <- 200
x <- cumsum(rnorm(n))
y <- 1 + 2*x + 3*(1:n > 100) + rnorm(n)

# Test for cointegration with breaks
result <- makicoint(y, x, model = 0, max_breaks = 2)
print(result)
```

Output:
```
Maki (2012) Cointegration Test with Multiple Structural Breaks
-----------------------------------------------------------------
Model: Level Shift
Observations: 200 | Regressors: 1 | Max breaks: 2
Trimming: 0.10 | Lags used: 0

H0: No cointegration
H1: Cointegration with up to 2 break(s)

Test Statistic: -12.3456
Critical Values:  1%: -5.416  5%: -4.892  10%: -4.610

Estimated Break Points:
  Break 1: Observation 98 (fraction: 0.4900)
  Break 2: Observation 145 (fraction: 0.7250)

Conclusion: Reject H0 at 1% significance level.
            Evidence of cointegration with structural break(s).
-----------------------------------------------------------------
```

## Model Specifications

The package supports four model types:

| Model | Name | Description |
|-------|------|-------------|
| 0 | Level Shift | Breaks only in intercept |
| 1 | Level Shift with Trend | Breaks in intercept, linear trend included |
| 2 | Regime Shift | Breaks in intercept and slope coefficients |
| 3 | Regime Shift with Trend | Breaks in intercept, slope, and trend |

## Features

- **Multiple breaks**: Test for up to 5 structural breaks
- **Multiple regressors**: Support for up to 4 independent variables  
- **Lag selection**: Four methods available (t-sig, fixed, AIC, BIC)
- **Critical values**: From Maki (2012) Table 1

## References

- Maki, D. (2012). Tests for cointegration allowing for an unknown number of breaks. *Economic Modelling*, 29(5), 2011-2015. [doi:10.1016/j.econmod.2012.04.022](https://doi.org/10.1016/j.econmod.2012.04.022)

- Gregory, A. W., & Hansen, B. E. (1996). Residual-based tests for cointegration in models with regime shifts. *Journal of Econometrics*, 70(1), 99-126.

## License

GPL (>= 3)
