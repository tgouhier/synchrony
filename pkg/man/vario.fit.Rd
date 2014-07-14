\name{vario.fit}
\alias{vario.fit}
\title{
vario.fit
}
\description{
Fit model to the empirical variogram
}
\usage{vario.fit (vario, bins, weights = rep(1, length(vario)),
                       type = c("spherical", "gaussian", "nugget", "linear", 
                       "exponential", "sill", "periodic", "hole"),
                       start.vals = list(c0 = 0, c1 = max(vario), 
                                          a = max(bins)/4, b=0.1, c=0.1),
                       control = list(maxit=10000)) 
}

\arguments{
  \item{vario}{Empirical variogram from \code{emp.vario} function}
  \item{bins}{Bins or lag distances from \code{emp.vario} function}
  \item{weights}{Vector of weights of the same length as \code{vario}. If \code{weights} is a vector 
  containing the number of points in each distance bin, the model will be fit via 
  weighted least squares with the weights corresponding to the proportion of points
  within each bin (i.e., weights sum to 1). Default is a vector of weights equal to 1}
  \item{type}{Type of variogram model to fit to the data. Default is \code{spherical}. 
  Other options are \code{gaussian}, \code{nugget}, \code{linear}, \code{exponential}, 
  \code{sill}, \code{periodic}, and \code{hole}
  }
\item{start.vals}{Named list containing the start values for the variogram model: 
\code{c0}: nugget, \code{c1}: sill, \code{a}: spatial range; \code{b}: slope; 
\code{c}: frequency}
  \item{control}{optional parameter for the \code{optim} function. 
  See \code{?optim} for details}
}
\value{
Return a named list containing the following variables:

\item{vario }{Empirical variogram values}
\item{bins }{Empirical variogram bins/lag distances}
\item{AIC }{AIC score of the model fit: \eqn{AIC=nlog\left(\frac{SSE}{n}\right)+2p}
 where \eqn{n} is the number of points in the variogram, \eqn{SSE=\sum{(\hat{x}_{i}-x_{i}})^2},
 and \eqn{p} is the number of parameters}
\item{RMSE }{Root Mean Square Error of the model fit: \eqn{\sqrt{\frac{SSE}{n}}}}
\item{params }{Named list containing the best model parameter estimates}
\item{fit }{Predicted variogram values from the model fit}
\item{nls.success }{did \code{nls} succeed?}
\item{convergence }{did \code{nls} or \code{optim} converge?}
}
\note{
Selecting proper initial values is critical for fitting a reasonable model to the 
empirical variogram. If these values are off, \code{nls} will fail and fall-back
functions will be used to determine the best parameter values that minimize the 
Root Mean Square Error (RMSE).
}

\author{
Tarik C. Gouhier (tarik.gouhier@gmail.com)
}
\seealso{
\code{\link{vario}, \link{vario.func}}
}
\examples{
# Load data
data(pisco.data)
# Environmental variogram
d=subset(pisco.data, subset=year==2000, select=c("latitude", "longitude", "upwelling"))
semiv=vario(data=d)
plot(semiv, xlab="Lag distance (km)")
mod.sph=vario.fit(semiv$vario, semiv$mean.bin.dist)
# Weighted least squares fit based on the number of points
mod.exp=vario.fit(semiv$vario, semiv$mean.bin.dist, 
                  weights=semiv$npoints/sum(semiv$npoints), 
                  type="expo")
mod.gau=vario.fit(semiv$vario, semiv$mean.bin.dist, type="gauss")
mod.lin=vario.fit(semiv$vario, semiv$mean.bin.dist, type="lin")
lines(semiv$mean.bin.dist, mod.sph$fit, col="red")
lines(semiv$mean.bin.dist, mod.exp$fit, col="black")
lines(semiv$mean.bin.dist, mod.gau$fit, col="blue")
lines(semiv$mean.bin.dist, mod.lin$fit, col="green")
legend(x="topleft", legend=paste(c("Spherical AIC:", "Exponential AIC:", 
                                   "Gaussian AIC:", "Linear AIC:"), 
                                   c(format(mod.sph$AIC, dig=2), 
                                   format(mod.exp$AIC, dig=2), 
                                   format(mod.gau$AIC, dig=2), 
       format(mod.lin$AIC, dig=2))), lty=1, col=c("red", "black", "blue", "green"), 
       bty="n")

# Correlogram
cover=subset(pisco.data, subset=year==2000, 
             select=c("latitude", "longitude", "mussel_abund"))
moran=vario(data=cover, type="moran")
mod.hol=vario.fit(moran$vario, moran$mean.bin.dist, 
                  type="hole", start.vals=list(c0=0.6, a=25, c1=0.01))
mod.per=vario.fit(moran$vario, moran$mean.bin.dist, type="period",
                  start.vals=list(a=1, b=3, c=0))
mod.lin=vario.fit(moran$vario, moran$mean.bin.dist, type="linear")
plot(moran, xlab="Lag distance (km)", ylim=c(-0.6, 0.8))
lines(moran$mean.bin.dist, mod.per$fit, col="red")
lines(moran$mean.bin.dist, mod.hol$fit, col="black")
lines(moran$mean.bin.dist, mod.lin$fit, col="blue")
legend(x="topleft", legend=paste(c("Periodic AIC:", "Hole AIC:", 
                                   "Linear AIC:"), 
                                   c(format(mod.per$AIC, dig=2), 
                                   format(mod.hol$AIC, dig=2), 
                                   format(mod.lin$AIC, dig=2))), 
                                   lty=1, col=c("red", "black", "blue"), bty="n")
}
\keyword{ variogram }
