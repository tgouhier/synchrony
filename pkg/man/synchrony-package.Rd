\name{synchrony-package}
\alias{synchrony-package}
\alias{synchrony}
\docType{package}
\title{
Methods for computing spatial, temporal, and spatiotemporal statistics
}
\description{
Methods for computing spatial, temporal, and spatiotemporal statistics including:
empirical univariate, bivariate and multivariate variograms; 
fitting variogram models; phase locking and synchrony analysis; 
generating autocorrelated and cross-correlated matrices.
}
\details{
\tabular{ll}{
Package: \tab synchrony\cr
Type: \tab Package\cr
Version: \tab 0.2.2\cr
Date: \tab 2014-03-06\cr
License: \tab GPL (>=2)\cr
URL: \tab http://synchrony.r-forge.r-project.org\cr
LazyLoad: \tab yes\cr
}
}
\author{
Tarik C. Gouhier (tarik.gouhier@gmail.com)

Maintainer: Tarik C. Gouhier (tarik.gouhier@gmail.com)
}
\references{
Bjornstad, O. N., and W. Falck. 2001. Nonparametric spatial covariance functions: 
Estimation and testing. \emph{Environmental and Ecological Statistics} 8:53-70.

Bjornstad, O. N., R. A. Ims, and X. Lambin. 1999. Spatial population dynamics: 
analyzing patterns and processes of population synchrony. \emph{Trends in Ecology & Evolution} 14:427-432.

Buonaccorsi, J. P., J. S. Elkinton, S. R. Evans, and A. M. Liebhold. 2001. 
Measuring and testing for spatial synchrony. \emph{Ecology} 82:1668-1679.

Cazelles, B., and L. Stone. 2003. Detection of imperfect population synchrony in 
an uncertain world. \emph{Journal of Animal Ecology} 72:953-968.

Fortin, M. J., and M. R. T. Dale. 2005. Spatial Analysis: A Guide for Ecologists. Cambridge University Press.

Gouhier, T. C., and F. Guichard. 2007. Local disturbance cycles and the maintenance 
of spatial heterogeneity across scales in marine metapopulations. \emph{Ecology} 88:647-657.

Gouhier, T. C., F. Guichard, and A. Gonzalez. 2010. Synchrony and stability of 
food webs in metacommunities. \emph{The American Naturalist} 175:E16-E34.

Gouhier, T. C., F. Guichard, and B. A. Menge. 2010. Ecological processes can 
synchronize marine population dynamics over continental scales. 
\emph{Proceedings of the National Academy of Sciences} 107:8281-8286.

Loreau, M., and C. de Mazancourt. 2008. Species synchrony and its drivers: 
Neutral and nonneutral community dynamics in fluctuating environments. 
\emph{The American Naturalist} 172:E48-E66.

Purves, D. W., and R. Law. 2002. Fine-scale spatial structure in a grassland community: 
quantifying the plant's eye view. 
\emph{Journal of Ecology} 90:121-129.

Vasseur, D. A. 2007. Environmental colour intensifies the Moran effect when 
population dynamics are spatially heterogeneous. \emph{Oikos} 116:1726-1736.

Zar, J. H. 1999. Biostatistical Analysis, Fourth edition. Prentice-Hall, Inc., 
Upper Saddle River, NJ.
}
\keyword{ synchrony }
\keyword{ concordance }
\keyword{ phase locking }
\keyword{ spatial synchrony }
\keyword{ correlogram }
\keyword{ variogram }

\examples{
# Compute phase synchrony
t1=runif(100)
t2=runif(100)
sync=phase.sync(t1, t2)
# Distribution of phase difference
hist(sync$deltaphase$mod_phase_diff_2pi)

# Compute concordant peaks
p=peaks(t1, t2, nrands=100)
# Find proportion of time steps where both time series peak together
p$peaks
# Plot (null) distribution of proportion of time steps where both time
# series peak together
hist(p$rand)
# p-value of observed value
p$pval

# Compute Kendall's W 
data(bird.traits)
(w=kendall.w(bird.traits))

# Community matrix for 20 species undergoing random fluctuations 
comm.rand=matrix(runif(100), nrow=5, ncol=20)
community.sync(comm.rand, nrands=10)
# Community matrix for 20 species undergoing synchronized fluctuations 
comm.corr=matrix(rep(comm.rand[,1], 20), nrow=5, ncol=20)
community.sync(comm.corr, nrands=10)
}