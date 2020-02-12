\name{innmmd.test}
\alias{innmmd.test}
\title{
Test of independence based on ranks of nearest neighbours and MMD
}
\description{
This function runs a test for independence among several random vectors based on ranks of nearest neighbors and MMD.
}
\usage{
innmmd.test(Xlist = NULL , Dlist = NULL, alpha = 0.05, 
dist.type = "euclidean", lp = 2, parameters = NULL,
B.parameters = 1000, B = 100, test.preset = NULL)
}
\arguments{
  \item{Xlist}{
A list matrices where each matrix corresponds to each random vector. Each row of a particular matrix is an observation on the corresponding random vector.
}
  \item{Dlist}{
A list of "dist" objects where each distance object corresponds to each random vector. It is required if \code{Xlist} is \code{NULL}.
}
  \item{alpha}{
A numeric value between 0 and 1 specifying the level of the test.
}
  \item{dist.type}{
The distance measure to be used. This must be one of \code{"euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"}. This is similar to parameter \code{method} of \code{dist} function.
}
  \item{lp}{
The power of Minkowski distance. This is similar to parameter \code{p} of \code{dist} function.
}
  \item{parameters}{
If \code{NULL}, then \code{innmmd.test} automatically calculates the parameters. Otherwise it uses the given set of parameters.
}
  \item{B.parameters}{
An integer value specifying the number of Monte-Carlo iterations to be used for automatic parameter generation procedure.
}
  \item{B}{
An integer value specifying the number of Monte-Carlo iterations made in the permutation test.
}
  \item{test.preset}{
If the list of numeric vectors is precalculated by the \code{innmmd.test.preset} function, supply that list here.
}
}
\details{
See reference.
}
\value{
A list of following components:
  \item{Tsum.stat}{
The value of the test statistic Tsum.
}
  \item{Tsum.max}{
The value of the test statistic Tmax.
}
  \item{Tsum.cutoff}{
The critical value for the test statistic Tsum.
}
  \item{Tmax.cutoff}{
The critical value for the test statistic Tmax.
}
  \item{Tsum.pvalue}{
P-value of the hypothesis test based on Tsum.
}
  \item{Tmax.pvalue}{
P-value of the hypothesis test based on Tmax.
}
}
\references{
Roy et al. (2019)
}
\author{
Angshuman Roy
}
\examples{
set.seed(1)
x <- rnorm(10)
y <- rnorm(10)
z <- matrix(x+y+rnorm(10,sd=0.1),ncol=1)
xlist <- list(x,y,z)
set.seed(1)
innmmd.test(xlist,alpha=0.05)

##
set.seed(1)
x <- rnorm(10)
y <- rnorm(10)
z <- matrix(x+y+rnorm(10,sd=0.1),ncol=1)
xlist <- list(x,y,z)
dlist <- lapply(xlist,dist)
set.seed(1)
innmmd.test(Dlist=dlist,alpha=0.05)
}
\keyword{nonparametric}
\keyword{htest}
