\name{kendall.w}
\alias{kendall.w}
\alias{concordance}
\title{
Kendall's W
}
\description{
Compute Kendall's coefficient of concordance (W)
}
\usage{
kendall.w (data, nrands = 0, type = 1, quiet = FALSE)
}
\arguments{
  \item{data}{
  matrix in wide format where each row represents a different sample and each column 
  represents a different variable.
}
  \item{nrands}{
  Number of randomizations to perform to determine significance. Default is 0.
}
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
}
\details{
Kendall's W is a non-parametric statistic that ranges from 0 to 1 and measures 
the level of agreement between multiple variables. When the number of observations \eqn{n>10}, 
its significance can be determined by using a \eqn{\chi^2} distribution with \eqn{df=n-1}.
Legendre (2005) shows that the \eqn{\chi^2} test is always too conservative (low power)
compared to the randomization test. Hence, both tests have been made available in
this function. The Monte Carlo randomizations are performed by shuffling the 
columns of the community matrix independently (Legendre 2005).
}
\value{
Returns a named list containing:

\item{w.uncorrected }{Kendall's W uncorrected for tied ranks}
\item{w.corrected }{Kendall's W corrected for tied ranks}
\item{pval }{p-value of Kendall's W}
\item{spearman.corr }{Spearman's ranked correlation}
\item{pval.rand }{p-value of Kendall's W based on randomization test. 
This variable is only returned if \code{nrands > 0}}
\item{rands }{randomizations. This variable is only returned if \code{nrands > 0}}
}
\references{
Buonaccorsi, J. P., J. S. Elkinton, S. R. Evans, and A. M. Liebhold. 2001. 
Measuring and testing for spatial synchrony. \emph{Ecology} 82:1668-1679.

Gouhier, T. C., and F. Guichard. 2007. Local disturbance cycles and the maintenance 
of spatial heterogeneity across scales in marine metapopulations. \emph{Ecology} 88:647-657.

Gouhier, T. C., F. Guichard, and A. Gonzalez. 2010. Synchrony and stability of 
food webs in metacommunities. \emph{The American Naturalist} 175:E16-E34.

Legendre, P. 2005. Species associations: the Kendall coefficient of concordance 
revisited. \emph{Journal of Agricultural, Biological, and Environmental Statistics} 10:226-245.

Purves, D. W., and R. Law. 2002. Fine-scale spatial structure in a grassland community: quantifying the plant's eye view. 
\emph{Journal of Ecology} 90:121-129.

Zar, J. H. 1999. Biostatistical Analysis, Fourth edition. Prentice-Hall, Inc., 
Upper Saddle River, NJ.
}
\author{
Tarik C. Gouhier (tarik.gouhier@gmail.com)
}

\examples{
data(bird.traits)
(w=kendall.w(bird.traits))
}
