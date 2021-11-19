#' Detecting linear trend changes for univariate time series
#'
#' The main function of the package \code{\link{trendsegmentR}}. This function estimates the number and locations of change-points in linear trend of noisy data. The estimated change-points may contain point anomalies (segments including only one data point) if any. It also returns the estimated signal, the best linear fit for each segment between a pair of adjacent change-points. The algorithm includes three steps, Tail-Greedy Unbalanced Wavelet (TGUW) transform (\code{\link{TGUW}}), thresholding (\code{\link{thresholding}}) and inverse TGUW transform (\code{\link{invTGUW}}).
#'
#' The algorithm is described in H. Maeng and P. Fryzlewicz (2021), Detecting linear trend changes in data sequences, preprint.
#'
#' @param x A data vector to be examined for change-point detection.
#' @param th.const Thresholding parameter used in \code{\link{thresholding}}. The default is 1.3 and the exact magnitude of the threshold is \code{sigma * th.const * sqrt(2 * log(n))} where \code{n} is the length of data sequence \code{x} and \code{sigma} is estimated by Median Absolute Deviation (MAD) method under the Gaussian assumption for noise.
#' @param p Proportion of all possible remaining merges which specifies the number of merges allowed in a single pass over the data. This is used in \code{\link{TGUW}} and the default is 0.01.
#' @param bal The minimum ratio of the length of the shorter region to the length of the entire merging region especially when the merges of Type 2 (merging one initial and a paired smooth coefficient) or of Type 3 (merging two sets of (paired) smooth coefficients) are performed. If \code{bal < 1/n} and \code{minsegL = 1}, point anomalies (segments including only one data point) can possibly be detected, otherwise, the detection of point anomalies is not guaranteed. The default is set to 0.
#' @param minsegL The minimum segment length of estimated signal returned by \code{trendsegment}. The default is set to 5 for avoiding too frequent change-points.
#' @param continuous If continuous=TRUE, the estimated signal returned by \code{trendsegment} is continuous at change-points, otherwise discontinuous at change-points. The default is set to FALSE.
#' @param connected If connected=TRUE, the \code{trendsegment} algorithm puts the connected rule above the \code{minsegL}, otherwise it makes keeping the \code{minsegL} a priority. The default is set to FALSE.
#' @return A list with the following.
#' \item{x}{The original input vector \code{x}.}
#' \item{est}{The estimated piecewise-linear signal of \code{x}.}
#' \item{no.of.cpt}{The estimated number of change-points.}
#' \item{cpt}{The estimated locations of change-points.}
#' @author Hyeyoung Maeng \email{h.maeng4@@lancaster.ac.uk}, Piotr Fryzlewicz \email{p.fryzlewicz@@lse.ac.uk}
#' @seealso \code{\link{TGUW}}, \code{\link{thresholding}}, \code{\link{invTGUW}}
#' @examples
#' x <- c(rep(0,100), seq(0, 4, length.out = 100), rep(3, 100), seq(3, -1, length.out=99))
#' n <- length(x)
#' x <- x + rnorm(n)
#' tsfit <- trendsegment(x = x, bal = 0)
#' tsfit
#'
#' plot(x, type = "b", ylim = range(x, tsfit$est))
#' lines(tsfit$est, col=2, lwd=2)
#' @importFrom stats mad
#' @export

trendsegment <- function(x, th.const = 1.3, p = .01, bal = 0, minsegL = 5, continuous = FALSE, connected = FALSE){

  n <- length(x)

  # The estimator of the standard deviation of noise.
  # The default is obtained by Median Absolute Deviation (MAD) method under the Gaussian assumption for noise.
  sigma <- stats::mad(diff(diff(x)))/sqrt(6)
  lambda <- sigma * sqrt(2 * log(n)) * th.const

  if (n == 1) {

    est <- x
    no.of.cpt <- 0
    cpt <- integer(0)

  } else {

    x.d <- TGUW(x, p)
    x.d.d <- thresholding(x.d, lambda=lambda, bal=bal, minsegL=minsegL, connected=connected)

    if(connected == FALSE){

      cpt <- finding.cp(x.d.d)
      no.of.cpt <- length(cpt)

      if(continuous == FALSE){

        est <- rep(NA, n)
        sgmts <- cbind(c(1, cpt+1), c(cpt, n))
        for(i in 1:dim(sgmts)[1]){
          a <- c(sgmts[i, 1]:sgmts[i, 2])
          est[a] <- stats::lm(x[a]~a)$fitted.values
        }

      }else{

        splns <- splines::bs(1:n, knots = cpt, degree = 1, intercept = TRUE)
        est <- stats::lm.fit(splns, x)$fitted.values

      }

    }else{

      est <- invTGUW(x.d.d) # use inverse TGUW directly

      cpt <- finding.cp(est)
      no.of.cpt <- length(cpt)

      #x <- est$x

      if(continuous == FALSE){

        est <- est$ts.coeffs

      }else{

        splns <- splines::bs(1:n, knots = cpt, degree = 1, intercept = TRUE)
        est <- stats::lm.fit(splns, x)$fitted.values

      }

    }

  }

  list(x=x, est=est, no.of.cpt=no.of.cpt, cpt=cpt)
}
