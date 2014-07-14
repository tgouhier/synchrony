\name{meancorr}
\alias{meancorr}
\title{
Compute mean column-wise correlation and determine its significance via Monte Carlo randomizations
}
\description{
Compute mean column-wise correlation and determine its significance via Monte Carlo randomizations. 
The Monte Carlo randomizations are performed by shuffling the columns of the community
matrix independently.
}
\usage{
meancorr (data, nrands = 0, alternative = c("two.tailed", "greater", "less"), 
                            method = c("pearson", "kendall", "spearman"), 
                            type = 1, quiet = FALSE, ...)
}
\arguments{
 \item{data}{
community matrix in wide format where each row contains the abundance at each time step
and each column corresponds to a different species.
}
 \item{nrands}{
number of randomizations to perform (default is 0)
}
\item{alternative}{Alternative hypothesis. Options include \code{greater} and 
\code{less} for the one-tailed test and \code{two.tailed}. Default is \code{two.tailed}.}
\item{method}{Method to compute correlation? Options include \code{pearson}, \code{kendall}, 
and \code{spearman}. Default is \code{pearson}}
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

\value{
Returns a named list containing:

\item{obs }{the observed mean correlation}
\item{rands }{the mean correlation for each randomization. 
This variable is only returned if \code{nrands > 0}.}
\item{pval }{p-value of observed mean correlation. 
This variable is only returned if \code{nrands > 0}.}
\item{alternative }{Alternative hypothesis.
This variable is only returned if \code{nrands > 0}.}
\item{method }{Method used to compute the mean correlation.}
}

\references{

Purves, D. W., and R. Law. 2002. Fine-scale spatial structure in a grassland community: quantifying the plant's eye view. 
\emph{Journal of Ecology} 90:121-129.
}

\author{
Tarik C. Gouhier (tarik.gouhier@gmail.com)
}

\examples{
# Community matrix for 20 species undergoing random fluctuations 
comm.rand=matrix(runif(100), nrow=5, ncol=20)
meancorr(comm.rand, nrands=20)$pval
# Community matrix for 20 species undergoing synchronized fluctuations 
comm.corr=matrix(rep(comm.rand[,1], 20), nrow=5, ncol=20)
meancorr(comm.corr, nrands=20)$pval
# On "real" data
data(bird.traits)
meancorr(bird.traits, nrands=20)$pval
}
