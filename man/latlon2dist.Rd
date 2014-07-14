\name{latlon2dist}
\alias{latlon2dist}
\title{
latlon2dist
}
\description{
Calculate distance between a pair of coordinates
}
\usage{
latlon2dist (coords)
}
\arguments{
  \item{coords}{
4-element vector of coordinates with format: \code{(lat1, lon1, lat2, lon2)}
}
}
\value{
Returns the great circle distance distance between the pair of coordinates
}
\author{
Tarik C. Gouhier (tarik.gouhier@gmail.com)
}
\seealso{
\code{\link{coord2dist}}
}
\examples{
coords=c(32, -125, 43, -130)
# Compute great circle distance
latlon2dist(coords)
}
\keyword{ dist }
