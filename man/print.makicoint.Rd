\name{print.makicoint}
\alias{print.makicoint}
\alias{summary.makicoint}
\title{Print and Summary Methods for makicoint Objects}
\usage{
\method{print}{makicoint}(x, ...)

\method{summary}{makicoint}(object, ...)
}
\arguments{
\item{x, object}{An object of class "makicoint" returned by \code{\link{makicoint}}.}

\item{...}{Additional arguments (currently unused).}
}
\value{
The input object is returned invisibly.
}
\description{
Methods for printing and summarizing the results of the Maki cointegration test.
}
\examples{
set.seed(123)
n <- 100
x <- cumsum(rnorm(n))
y <- 1 + 2*x + rnorm(n)

result <- makicoint(y, x, model = 0, max_breaks = 1)
print(result)
summary(result)

}
\seealso{
\code{\link{makicoint}}
}
