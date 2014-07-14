\name{plot.synchrony}
\alias{plot.synchrony}
\title{
  Plot \code{synchrony} objects
}
\description{
  Plot \code{synchrony} objects
}
\usage{
  \method{plot}{synchrony} (x, main = "", xlab = "Values from randomizations", 
                            ylab = "Frequency", line.col = "red", lty = 2, 
                            lwd = 1, col = "grey", \dots)
}
\arguments{
  \item{x}{
    \code{synchrony} object
  }
  \item{main}{
  main title of the figure
  }
  \item{xlab}{
    xlabel of the figure. Default is "Values from randomizations"
  }
  \item{ylab}{
    ylabel of the figure. Default is "Frequency"
  }
  \item{line.col}{
    color of the vertical line indicating the value observed in the data. Default is "red"
  }
  \item{lty}{
    line type. Default is 2 or dashed
  }
  \item{lwd}{
    line width. Default is 1
  }
  \item{col}{
    color of the bars. Default is grey  
  }
  \item{\dots}{
    other graphical parameters.
  }
}

\author{
  Tarik C. Gouhier (tarik.gouhier@gmail.com)
}

\examples{
comm.rand=matrix(runif(100), nrow=5, ncol=20)
comm.rand.sync=community.sync(comm.rand, nrands=20)
plot(comm.rand.sync)
}
