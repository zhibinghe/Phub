% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Mphub.R
\name{mphub}
\alias{mphub}
\title{\code{\link{phub}} with multiple initial points, both \code{A} and \code{rho} are randomly drawn from Unif(0,1)}
\usage{
mphub(
  G,
  M,
  lam,
  rep = 10,
  pen.type = c("plog", "plasso", "log"),
  iter.max = 1000,
  tol = 1e-04
)
}
\arguments{
\item{G}{observed group data}

\item{M}{user-specified number of components}

\item{lam}{tuning parameter for component selection, degenrate to a standard EM without penalty if lam=0}

\item{rep}{number of repetition}

\item{pen.type}{type of penalty, including 'log': penalization for logarithm of all components;
'plog': penalization for logarithm of partial components except the null component;
'plasso': penalization for lasso form of partial components}

\item{iter.max}{maximum iteration steps}

\item{tol}{threshold for shrinking rho to 0}
}
\value{
a list of components
\item{A}{M*n matrix containg correlation among nodes}
\item{rho}{M-dimensional vector containing estimated hub weights}
\item{l}{log-likelihood}
\item{iteration}{number of iterations to converge}
}
\description{
\code{\link{phub}} with multiple initial points, both \code{A} and \code{rho} are randomly drawn from Unif(0,1)
}
\examples{
n0 = 10; n=100; T=1000
A0 = GenA(100,n0,0.4,0.1,rep(0.05,n))
G0 = GenG(A0,T,c(0.2,rep(0.8/n0,n0)))
M = 15; mphub(G0,M,0.02)$rho
\donttest{M = 20; mphub(G0,M,0.035)}
\donttest{M = 50; mphub(G0,M,0.019)}

}
