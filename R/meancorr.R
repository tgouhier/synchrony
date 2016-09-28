#' Compute mean column-wise correlation and determine its significance via Monte
#' Carlo randomizations
#'
#' Compute mean column-wise correlation and determine its significance via Monte
#' Carlo randomizations. The Monte Carlo randomizations are performed by
#' shuffling the columns of the community matrix independently.
#'
#' @note n x m matrix with n=time step, m=species
#' @author Tarik C. Gouhier (tarik.gouhier@@gmail.com)
#'
#' @param data community matrix in wide format where each row contains the
#'   abundance at each time step and each column corresponds to a different
#'   species.
#' @param nrands number of randomizations to perform (default is 0)
#' @param alternative Alternative hypothesis. Options include \code{greater} and
#'   \code{less} for the one-tailed test and \code{two.tailed}. Default is
#'   \code{two.tailed}.
#' @param method Method to compute correlation? Options include \code{pearson},
#'   \code{kendall}, and \code{spearman}. Default is \code{pearson}
#' @param type Randomization method. The \code{type=1} method randomly shuffles
#'   each column of the data matrix, thus destroying both the autocorrelation
#'   structure of each column and the cross-correlation between columns. The
#'   \code{type=2} method shifts each column of the data matrix by a random
#'   amount, thus preserving the autocorrelation structure of each column but
#'   destroying the cross-correlation between columns (Purves and Law 2002).
#'   Default is \code{type=1}
#' @param quiet TRUE suppresses the progress bar. Default is FALSE.
#' @param ... arguments passed to \code{\link{cor}}
#'
#' @return Returns a named list containing:
#' \item{obs}{The observed mean correlation.}
#' \item{rands}{The mean correlation for each randomization.
#'              This variable is only returned if \code{nrands > 0}.}
#' \item{pval}{p-value of observed mean correlation.
#'             This variable is only returned if \code{nrands > 0}.}
#' \item{alternative}{Alternative hypothesis.
#'                    This variable is only returned if \code{nrands > 0}.}
#' \item{method}{Method used to compute the mean correlation.}
#'
#' @references
#' Purves, D. W., and R. Law. 2002. Fine-scale spatial structure in a grassland community: quantifying the plant's eye view.
#' \emph{Journal of Ecology} 90:121-129.
#'
#' @examples
#' # Community matrix for 20 species undergoing random fluctuations
#' comm.rand=matrix(runif(100), nrow=5, ncol=20)
#' meancorr(comm.rand, nrands=20)$pval
#'
#' # Community matrix for 20 species undergoing synchronized fluctuations
#' comm.corr=matrix(rep(comm.rand[,1], 20), nrow=5, ncol=20)
#' meancorr(comm.corr, nrands=20)$pval
#'
#' # On "real" data
#' data(bird.traits)
#' meancorr(bird.traits, nrands=20)$pval

#' @export
meancorr <- function (data, nrands = 0,
                            alternative=c("two.tailed", "greater", "less"),
                            method=c("pearson", "kendall", "spearman"),
                      type=1, quiet = FALSE, ...) {
  data=as.matrix(data)
  results=list()
  methods=c("pearson", "kendall", "spearman")
  method=match.arg(method, methods)

  alternatives=c("two.tailed", "greater", "less")
  alternative=match.arg(tolower(alternative), alternatives)

  results$obs=meancorr.aux (data, method=method, ...)

  if (nrands > 0) {
    nr=NROW(data)
    nc=NCOL(data)
    if (!quiet)
      prog.bar=txtProgressBar(min = 0, max = nrands, style = 3)
    results$rands=numeric(length=nrands+1)*NA
    for (i in 1:nrands) {
      if (type==1)
        rand.mat=apply(data, 2, sample)
      else {
        lags=sample(1:nr, size=nc, replace=TRUE)
        rand.mat=mlag(data, lags)
      }
      results$rands[i]=meancorr.aux(rand.mat, method=method, ...)
      if (!quiet)
        setTxtProgressBar(prog.bar, i)
    }
    results$rands[nrands+1]=results$obs

    if (alternative == "two.tailed") {
      pvals=sum(abs(results$rands) >= abs(results$obs))/(nrands+1)
    }
    else {
      if (alternative=="greater")
        pvals=sum(results$rands >= results$obs)/(nrands+1)
      else
        pvals=sum(results$rands <= results$obs)/(nrands+1)
    }

    results$pval=pvals
    results$alternative=alternative
  }
  results$method=method
  class(results)="synchrony"
  return (results)
}

meancorr.aux <- function (data, method=method, ...) {
  mean.corr=suppressWarnings(cor(data, method=method, ...))
  mean.corr=mean(mean.corr[lower.tri(mean.corr)], na.rm=TRUE)
  return (mean.corr)
}
