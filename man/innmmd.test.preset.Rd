\name{innmmd.test.preset}
\alias{innmmd.test.preset}
\title{
A list of numeric vectors pre-calculated for \code{innmmd.test} function
}
\description{
When \code{innmmd} function has to be calculated repeatedly, we can save a lot of computation by pre-calculating some numeric vectors that \code{innmmd.test} function uses by using \code{innmmd.test.preset} function.
}
\usage{
innmmd.test.preset(n, d, parameters = NULL, B.parameters = 1000)
}
\arguments{
  \item{n}{
Sample size.
}
  \item{d}{
Number of random vectors.
}
  \item{parameters}{
If \code{NULL}, then \code{innmmd.test.preset} automatically calculates the parameters. Otherwise it calculates preset using the given set of parameters.
}
  \item{B.parameters}{
An integer value specifying the number of Monte-Carlo iterations to be used for automatic parameter generation procedure.
}
}
\details{
See reference.
}
\value{
A list of numeric vectors.
}
\references{
Roy et al. (2019)
}
\author{
Angshuman Roy
}
\examples{
#automated parameters
set.seed(0)
innmmd.test.preset(10,2,B.parameters=1000)

#manual parameters
innmmd.test.preset(10,2,parameters=c(0.3,0.5,1))
}
\keyword{nonparametric}
