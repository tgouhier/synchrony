\name{vario}
\alias{vario}
\title{
vario
}
\description{
Compute the empirical variogram and determine its significance via Monte Carlo randomizations
}
\usage{vario (n.bins = 20, size.bins = NULL, extent = 0.5, data, data2 = NULL, 
                  is.latlon = TRUE, is.centered = FALSE, nrands = 0, 
                  type = c("semivar", "cov", "pearson", 
                  "spearman", "kendall", "moran", "geary"),
                  alternative = c("one.tailed", "two.tailed"),
                  mult.test.corr = c("none", "holm", "hochberg", "bonferroni"),
                  quiet = FALSE)
}

\arguments{
  \item{n.bins}{Number of bins or lag distances. This argument is only used when 
  \code{size.bins=NULL}
  } 
  \item{size.bins}{Size of bins in units of distance (e.g., km). When specified, this argument overrides
  \code{n.bins}. Default is \code{NULL}}
  \item{extent}{Proportion of the spatial extent of the data over which to compute 
  the variogram. Default is 0.5 to limit potentially spurious results due to the 
  limited number of data points at large lag distances.
  }
  \item{data}{\code{n} x \code{m} matrix containing y-coordinates (or latitude), 
              x-coordinates (or longitude), and values. The values can either be
              a single column of observations at each site for univariate variograms
              or a matrix of observations at each site for multivariate variograms 
              (e.g., to compute spatial synchrony).
              }
  \item{data2}{\code{n} x \code{m} matrix containing y-coordinates (or latitude), 
              x-coordinates (or longitude), and values for second variable. The values 
              can either be a single column of observations at each site for univariate 
              variograms or a matrix of observations at each site for multivariate 
              variograms (e.g., to compute spatial synchrony).
              }
  \item{is.latlon}{Are coordinates latitudes/longitudes? Default is \code{TRUE}
}

\item{is.centered}{Should the variogram be centered by subtracting the regional mean
from each value? If so, the zero-line represents the regional mean. Default is \code{FALSE}
}
\item{nrands}{Number of randomizations to determine statistical significance 
  of variogram. Default is 0.
  } 
  \item{type}{Type of variogram to compute. Default is \code{semivar} for semivariance. 
  Other options include \code{cov} for covariance, \code{pearson} for Pearson
correlation, \code{spearman} for Spearman correlation, \code{kendall} for Kendall correlation,
\code{moran} for Moran's I, and \code{geary} for Geary's C
}
\item{alternative}{Conduct a one-tailed or a two-tailed test? Note that the statistical test
is to determine whether the local value within each lag distance is different from the regional
mean. If the variogram is centered, the null hypothesis is that the local values are equal to zero. 
If the variogram is not centered, the null hypothesis is that the local values are equal to the 
regional mean. Default is \code{one.tailed}}
\item{mult.test.corr}{Correct for multiple tests? Default is \code{"none"}. Other options include
\code{holm}, \code{hochberg} and \code{bonferroni}
}
\item{quiet}{
Suppress progress bar when set to \code{TRUE}. Default is \code{FALSE}
}
}
\details{
This function can be used to compute univariate correlograms using Moran's I, 
Geary's C, and the covariance function or variograms using the semivariance function. 
Multivariate (Mantel) correlograms can also be computed using the covariance function, 
Pearson's, Spearman's or Kendall's correlation coefficients. Cross-correlograms/variograms 
between \code{data1} and \code{data2} can be computed with the covariance function,
Pearson's, Spearman's or Kendall's correlation coefficients for multivariate
variograms and Moran's I, Geary's C, the covariance function, or semivariance 
for univariate variograms.
}

\value{
Returns a named list containing the following variables:

\item{bins }{Center of each lag/bin}
\item{mean.bin.dist }{Mean distance of each lag/bin}
\item{vario }{Variogram values in each lag/bin}
\item{npoints }{Number of pairs of points in each lag/bin}
\item{metric }{Type of variogram computed}
\item{is.centered }{Is the variogram centered?}
\item{regional.mean }{Regional mean value}
\item{pvals }{p-value for each lag/bin. 
This variable is only returned if \code{nrands > 0}.}
\item{rands }{\code{nrands} x \code{n.bins} matrix of randomizations.
This variable is only returned if \code{nrands > 0}.}
\item{alternative }{One-tailed or two-tailed test?
This variable is only returned if \code{nrands > 0}.}
\item{mult.test.corr }{Correct for multiple tests?
This variable is only returned if \code{nrands > 0}.}
\item{is.multivar }{Was the analysis performed on multivariate data?}
}
\references{
Bjornstad, O. N., and W. Falck. 2001. Nonparametric spatial covariance functions: 
Estimation and testing. \emph{Environmental and Ecological Statistics} 8:53-70.

Bjornstad, O. N., R. A. Ims, and X. Lambin. 1999. Spatial population dynamics: 
analyzing patterns and processes of population synchrony. \emph{Trends in Ecology & Evolution} 14:427-432.

Fortin, M. J., and M. R. T. Dale. 2005. Spatial Analysis: A Guide for Ecologists. Cambridge University Press.
}

\author{
Tarik C. Gouhier (tarik.gouhier@gmail.com)
}
\seealso{
\code{\link{vario.func}}
}
\examples{
data(pisco.data)
d=subset(pisco.data, subset=year==2000, select=c("latitude", "longitude", "sst"))
semiv=vario(data=d)
moran=vario(data=d, type="moran", nrands=100)
par(mfrow=c(2,1), mar=c(4.2, 4, 1, 1))
plot(semiv$mean.bin.dist, semiv$vario, xlab="Lag distance (km)", ylab="Semivariance")
plot(moran$mean.bin.dist, moran$vario, xlab="Lag distance (km)", ylab="Moran's I", t="l")
points(moran$mean.bin.dist[moran$pvals >= 0.05], moran$vario[moran$pvals >= 0.05], 
       bg="white", pch=21)
points(moran$mean.bin.dist[moran$pvals < 0.05], moran$vario[moran$pvals < 0.05], 
       bg="black", pch=21)
abline(h=0, lty=2)

# Compute spatial synchrony
d.upw=subset(pisco.data, select=c("latitude", "longitude", "year", "upwelling"))
d.cov=subset(pisco.data, select=c("latitude", "longitude", "year", "mussel_abund"))
# Reshape the data
d.upw.wide=reshape(data=d.upw, idvar=c("latitude", "longitude"), timevar=c("year"), 
                   direction="wide")
d.cov.wide=reshape(data=d.cov, idvar=c("latitude", "longitude"), timevar=c("year"), 
                   direction="wide")
# Generate variograms
v.upw=vario(n.bins=12, data=d.upw.wide, type="pearson", extent=1, nrands=999)
v.cov=vario(n.bins=12, data=d.cov.wide, type="pearson", extent=1, nrands=999)
## Fit variograms
v.cov.per=vario.fit(v.cov$vario, v.cov$mean.bin.dist, type="period", 
                    start.vals=list(a=1, b=3, c=0))
v.upw.lin=vario.fit(v.upw$vario, v.upw$mean.bin.dist, type="linear")

par(mfrow=c(2,1))
plot(v.cov, xlab="Lag distance (km)", bg.sig="red", col.nonsig="red", 
     main="Mussel cover", 
     rug=TRUE, ylim=c(-0.3, 0.3))
lines(v.cov$mean.bin.dist, v.cov.per$fit, col="red")
plot(v.upw, xlab="Lag distance (km)", bg.sig="blue", col.nonsig="blue", 
     main="Upwelling", rug=TRUE)
lines(v.upw$mean.bin.dist, v.upw.lin$fit, col="blue")
}
\keyword{ variogram }
