\name{vario.func}
\alias{vario.func}
\title{
  vario.func
}
\description{
Compute the empirical variogram values for each bin
}
\usage{vario.func (x, y, glob.mean, glob.sd, glob.N, is.multivar = FALSE,
                        type = c("semivar", "cov", "pearson", 
                                 "spearman", "kendall", "moran", "geary"))
}

\arguments{
  \item{x}{First set of sites within bin/lag distance}
  \item{y}{Second set of sites within bin/lag distance}
  \item{glob.mean}{Global mean}
  \item{glob.sd}{Global standard deviation}
  \item{glob.N}{Global number of points}
  \item{is.multivar }{Is the data multivariate? Default is \code{FALSE}}
  \item{type}{Type of variogram to compute. Default is \code{semivar} for semivariance. 
  Other options include \code{cov} for covariance, \code{pearson} for Pearson
correlation, \code{spearman} for Spearman correlation, \code{kendall} for Kendall correlation,
\code{moran} for Moran's I, and \code{geary} for Geary's C  }
}
\value{
  Return the value.}
\author{
  Tarik C. Gouhier (tarik.gouhier@gmail.com)
}
\seealso{
  \code{\link{vario}}
}
\examples{
  # Internal function used by vario
}
\keyword{ variogram }
