#' Detecting linear trend changes and point anomalies for univariate time series
#'
#' The main function of the package \code{\link{trendsegmentR}}. This function estimates the number and locations of change-points in linear trend of noisy data. The estimated change-points may contain point anomalies if any. It also returns the estimated signal, the best linear fit for each segment between a pair of adjacent change-points. The algorithm includes three steps, Tail-Greedy Unbalanced Wavelet (TGUW) transform (\code{\link{TGUW}}), thresholding (\code{\link{thresholding}}) and inverse TGUW transform (\code{\link{invTGUW}}).
#'
#' The algorithm is described in H. Maeng and P. Fryzlewicz (2019), Detecting linear trend changes and point anomalies in data sequences, preprint.
#'
#' @param x A data vector to be examined for change-point detection.
#' @param th.const Thresholding parameter used in \code{\link{thresholding}}. The default is 1.3 and the exact magnitude of the threshold is \code{sigma * th.const * sqrt(2 * log(n))} where \code{n} is the length of data sequence \code{x} and \code{sigma} is estimated by Median Absolute Deviation (MAD) method under the Gaussian assumption for noise.
#' @param p Proportion of all possible remaining merges which specifies the number of merges allowed in a single pass over the data. This is used in \code{\link{TGUW}} and the default is 0.01.
#' @param bal The minimum ratio of the length of the shorter region to the length of the entire merging region especially when the merges of Type 2 (merging one initial and a paired smooth coefficient) or of Type 3 (merging two sets of (paired) smooth coefficients) are performed. If \code{bal < 1/n}, point anomalies can be detected, otherwise, the detection of point anomalies is not guaranteed. The default is set to 0 for detecting both point anomalies and linear trend changes.
#' @return A list with the following.
#' \item{x}{The original input vector \code{x}.}
#' \item{est}{The estimated piecewise-linear signal of \code{x}.}
#' \item{no.of.cpt}{The estimated number of change-points.}
#' \item{cpt}{The estimated locations of change-points.}
#' @author Hyeyoung Maeng, \email{h.maeng@@lse.ac.uk}
#' @seealso \code{\link{TGUW}}, \code{\link{thresholding}}, \code{\link{invTGUW}}
#' @examples
#' x <- c(rep(0,100), -4, seq(0, 4, length.out = 100), rep(3, 100), seq(3, -1, length.out=99))
#' n <- length(x)
#' x <- x + rnorm(n)
#' tsfit <- trendsegment(x = x, bal = 0)
#' tsfit
#'
#' plot(x, type = "b", ylim = range(x, tsfit$est))
#' lines(tsfit$est, col=2, lwd=2)
#' @importFrom stats mad
#' @export





trendsegment <- function(x, th.const = 1.3, p = .01, bal = 0){

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
    x.d.d <- thresholding(x.d, lambda=lambda, bal=bal)
    est <- invTGUW(x.d.d)

    cpt <- finding.cp(est)
    no.of.cpt <- length(cpt)

    x <- est$x
    est <- est$ts.coeffs
  }

  list(x=x, est=est, no.of.cpt=no.of.cpt, cpt=cpt)
}
