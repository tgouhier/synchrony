\name{peaks}
\alias{peaks}
\title{
Find the proportion of local minima/maxima common to both time series and compute 
its significance via Monte Carlo randomizations
}
\description{
Find the proportion of local minima/maxima common to both time series and compute 
its significance via Monte Carlo randomizations
}
\usage{
peaks (t1, t2, nrands = 0, type = 1, quiet = FALSE)
}
\arguments{
  \item{t1}{
time series 1 in matrix format (\code{n} rows x 2 columns). The first column 
should contain the time steps and the second column should contain the values. 
If \code{t1} is a column vector instead of a matrix, then it will be automatically 
converted to a matrix with column 1 corresponding to a time index ranging 
from 1 to the length of \code{t1}
}
  \item{t2}{
time series 2 in matrix format (\code{n} rows x 2 columns). The first column should contain the time 
steps and the second column should contain the values. If \code{t2} is a column
vector instead of a matrix, then it will be automatically converted to matrix 
with column 1 corresponding to a time index ranging from 1 to the length of \code{t2}.
}
  \item{nrands}{
number of randomizations. Default is \code{0}.
}
  \item{type}{
Randomization method. The \code{type=1} method randomly shuffles each time series, 
thus destroying both the autocorrelation structure of each time series and their 
cross-correlation. The \code{type=2} method shifts each time series by a random amount, 
thus preserving the autocorrelation structure but destroying the cross-correlation 
between the time series (Purves and Law 2002). Default is \code{type=1}
}

\item{quiet}{
Suppress progress bar when set to \code{TRUE}. Default is \code{FALSE}
}

}
\value{
Returns a named list containing:

\item{pval }{p-value computed by randomly shuffling both time series \code{nrands} times}
\item{rands }{proportion of local minima/maxima common to both time series for each randomization}
\item{obs }{proportion of local minima/maxima common to both time series in the observed dataset}
\item{index }{indices of local minima/maxima common to both time series in the observed dataset}
}
\references{
Buonaccorsi, J. P., J. S. Elkinton, S. R. Evans, and A. M. Liebhold. 2001. 
Measuring and testing for spatial synchrony. \emph{Ecology} 82:1668-1679.

Purves, D. W., and R. Law. 2002. Fine-scale spatial structure in a grassland community: quantifying the plant's eye view. 
\emph{Journal of Ecology} 90:121-129.
}
\author{
Tarik C. Gouhier (tarik.gouhier@gmail.com)
}

\examples{
t1=runif(100)
t2=runif(100)
(p=peaks(t1, t2))

}
