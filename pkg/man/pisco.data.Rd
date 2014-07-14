\name{pisco.data}
\alias{pisco.data}
\docType{data}
\title{
PISCO multi-year and spatially-explicit mussel and environmental dataset
}
\description{
Contains the mean annual chl-a concentration, sea surface temperature, 
upwelling currents, and mussel abundance at 48 intertidal sites
along the West Coast of the United States from 2000-2003.
}
\usage{data(pisco.data)}
\format{
  A data frame with 192 observations on the following 7 variables.
  \describe{
    \item{\code{latitude}}{latitude (degrees North)}
    \item{\code{longitude}}{longitude (degrees West)}
    \item{\code{chl}}{mean annual remote sensed chlorophyll-a concentration}
    \item{\code{sst}}{mean annual remote sensed sea surface temperature}
    \item{\code{upwelling}}{mean annual remote sensed upwelling currents}
    \item{\code{mussel_abund}}{mean annual mussel cover (\emph{Mytilus californianus})}
    \item{\code{year}}{sampling year}
  }
}
\references{
Gouhier, T. C., F. Guichard, and B. A. Menge. 2010. Ecological processes can synchronize marine population dynamics over continental scales. \emph{Proceedings of the National Academy of Sciences} 107:8281-8286.
}
\examples{
data(pisco.data)
}
\keyword{datasets}
