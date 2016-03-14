\name{community.sync}
\alias{community.sync}
\title{
Compute community-wide synchrony and its significance via Monte Carlo randomizations
}
\description{
Compute community-wide synchrony and its the significance via Monte Carlo randomizations. 
If all species fluctuate in perfect unison, the community-wide synchrony will be 1. 
If species undergo uncorrelated fluctuations, the community-wide synchrony will be 1/S.
The Monte Carlo randomizations are performed by shuffling the columns of the community
matrix independently. This function also returns the mean correlation between the columns
of the matrix.
}
\usage{
community.sync (data, nrands = 0, method = c("pearson", "kendall", "spearman"), 
                alternative = c("greater", "less"), type = 1, quiet = FALSE, ...)
}
\arguments{
 \item{data}{
community matrix in wide format where each row contains the abundance at each time step
and each column corresponds to a different species.
}
 \item{nrands}{
number of randomizations to perform (default is 0)
}
\item{method}{Method to compute mean correlation between columns? 
Options include \code{pearson}, \code{kendall}, and \code{spearman}. Default is \code{pearson}}

\item{alternative}{Alternative hypothesis. Options are
\code{less} and \code{greater}. Default is \code{greater}}

\item{type}{
Randomization method. The \code{type=1} method randomly shuffles each column of the data
matrix, thus destroying both the autocorrelation structure of each column and the cross-correlation
between columns. The \code{type=2} method shifts each column of the data matrix
by a random amount, thus preserving the autocorrelation structure of each column but destroying
the cross-correlation between columns (Purves and Law 2002). Default is \code{type=1}
}

\item{quiet}{
Suppress progress bar when set to \code{TRUE}. Default is \code{FALSE}
}

\item{\dots}{
Other parameters to \code{\link{cor}} function.
}
}

\details{
Loreau and de Mazancourt (2008) show that community-wide synchrony \eqn{\varphi} can be quantified
by computing the temporal variance \eqn{\sigma_{x_T}^2} of the community time series
\eqn{x_T(t)=\sum{x_i(t)}} and the sum of the temporal standard deviation of the time series
across all species \eqn{\left(\sum{\sigma_{x_i}}\right)^2} such that:
\eqn{\varphi=\frac{\sigma_{x_T}^2}{\left(\sum{\sigma_{x_i}}\right)^2}}
}

\value{
Returns a named list containing:

\item{obs }{the observed community synchrony}
\item{meancorr }{the mean correlation between the columns of the matrix}
\item{rands }{the community synchrony value the randomizations. 
This variable is only returned if \code{nrands > 0}.}
\item{pval }{p-value of observed community synchrony. 
This variable is only returned if \code{nrands > 0}.}
\item{alternative }{Alternative hypothesis. This variable is only returned if \code{nrands > 0}.}
}

\references{
Loreau, M., and C. de Mazancourt. 2008. Species synchrony and its drivers: 
Neutral and nonneutral community dynamics in fluctuating environments. 
\emph{The American Naturalist} 172:E48-E66.

Purves, D. W., and R. Law. 2002. Fine-scale spatial structure in a grassland community: quantifying the plant's eye view. 
\emph{Journal of Ecology} 90:121-129.
}

\author{
Tarik C. Gouhier (tarik.gouhier@gmail.com)
}

\examples{
# Community matrix for 20 species undergoing random fluctuations 
comm.rand=matrix(runif(100), nrow=5, ncol=20)
community.sync(comm.rand, nrands=20)$pval
# Community matrix for 20 species undergoing synchronized fluctuations 
comm.corr=matrix(rep(comm.rand[,1], 20), nrow=5, ncol=20)
community.sync(comm.corr, nrands=20)$pval
# On "real" data
data(bird.traits)
community.sync(bird.traits, nrands=20)$pval
}
