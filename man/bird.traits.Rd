\name{bird.traits}
\alias{bird.traits}
\docType{data}
\title{
bird trait dataset
}
\description{
Contains the wing length, tail length, and bill length from 12 birds
}
\usage{data(bird.traits)}
\format{
  A data frame with 12 observations (birds) on the following 3 variables.
  \describe{
    \item{\code{wing.length}}{a numeric vector containing wing length in cm}
    \item{\code{tail.length}}{a numeric vector containing tail length in cm}
    \item{\code{bill.length}}{a numeric vector containing bill length in cm}
  }
}
\details{
Dataset from Zar (1999; page 444)
}
\source{
Zar, J. H. 1999. Biostatistical Analysis, Fourth edition. Prentice-Hall, Inc., Upper Saddle River, NJ.
}
\examples{
data(bird.traits)
(w=kendall.w(bird.traits))
}
\keyword{datasets}
