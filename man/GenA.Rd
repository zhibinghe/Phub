% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Generation.R
\name{GenA}
\alias{GenA}
\title{Simulate a simple pattern 'Adjacent matrix'}
\usage{
GenA(n, n0, pf, pr, pi)
}
\arguments{
\item{n}{number of nodes}

\item{n0}{number of true nonzero components}

\item{pf}{probability of following}

\item{pr}{probality of rejection}

\item{pi}{a n-dimensional vector containg the probability of appearance in nonleader group}
}
\value{
a (n0+1)*n matrix; the first row is \code{pi}
}
\description{
Simulate a simple pattern 'Adjacent matrix'
}
\examples{
pi = rep(0.05,100)
A0 = GenA(100,5,0.4,0.1,pi)

}
