\name{coord2dist}
\alias{coord2dist}
\title{
coord2dist
}
\description{
Calculate distance between all pairs of sites
}
\usage{
coord2dist (coords, is.latlon = TRUE, lower.tri = TRUE)
}
\arguments{
  \item{coords}{
\code{n} x 4 matrix of coordinates consisting of {lat or y, lon or x} pairs for each each site
}
  \item{is.latlon}{
are coordinates latitudes/longitudes? Default is \code{TRUE}
}
  \item{lower.tri}{
Return lower triangular part of the distance matrix? Default is \code{TRUE}
}
}
\value{
Returns the distance between all pairs of sites
}
\author{
Tarik C. Gouhier (tarik.gouhier@gmail.com)
}
\examples{
coords=rbind(c(32, -125), c(43, -130))
# Compute great circle distance
coord2dist(coords)
}
\keyword{ dist }
