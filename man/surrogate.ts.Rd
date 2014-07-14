\name{surrogate.ts}
\alias{surrogate.ts}
\title{
Create surrogate time series via Markov process
}
\description{
Create surrogate time series with the same short-term time correlation and overall
temporal pattern as the original time series using the Markov process described by
Cazelles and Stones (2003)
}
\usage{
surrogate.ts (ts, distr.ts = NULL, trans.ts = NULL, nbreaks = 10)
}
\arguments{
  \item{ts}{
time series in matrix format (\code{n} rows x 2 columns). The first column should contain the time 
steps and the second column should contain the values. If \code{ts} is a column vector 
instead of a matrix, then it will be automatically converted to a matrix with 
column 1 corresponding to a time index ranging from 1 to the length of \code{ts}}
  \item{distr.ts}{
binning of time series values. This parameter must be specified 
if \code{trans.ts} is not set to \code{NULL}. Default is \code{NULL}.}
  \item{trans.ts}{
transition matrix from bin \eqn{i} to bin \eqn{j}. Default is \code{NULL}.}
  \item{nbreaks}{
number of bins to use to group the time series values. Default is \code{10}.}
}
\details{
The values of the time series \eqn{x_n} are grouped into \code{nbreak} equally-sized bins. 
The transition matrix \eqn{M_{ij}} describing the probability of \eqn{x_{n+1}} belonging to 
bin \eqn{j} based on \eqn{x_n} belonging to bin \eqn{i} is defined using the relative
frequencies of the data such that: 
\eqn{M_{ij}=Pr(x_{n+1} \in b_{j} | x_{n} \in b_{i})}. The surrogate time series is then constructed
by randomly selecting a starting value and randomly selecting the next value from the proper bin
based on the transition matrix. This process is repeated until the surrogate time series has 
the same length as the original time series.
}

\value{
Returns a named list containing:

\item{surr.ts }{surrogate time series in matrix format}
\item{trans }{transition matrix \eqn{M_{ij}}}
\item{distr }{binning of time series values}
}
\references{
Cazelles, B., and L. Stone. 2003. Detection of imperfect population synchrony 
in an uncertain world. \emph{Journal of Animal Ecology} 72:953-968.
}
\author{
Tarik C. Gouhier (tarik.gouhier@gmail.com)
}

\seealso{
\code{\link{phase.sync}}

}
\examples{
t1=runif(100)
surr.t1=surrogate.ts(ts=t1)
plot(t1, t="l")
lines(surr.t1$surr.ts, col="red")
}
