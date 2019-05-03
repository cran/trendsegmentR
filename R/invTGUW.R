#' Inverse Tail-Greedy Unbalanced Wavelet (TGUW) transformation
#'
#' This function performs the inverse TGUW transformation by undoing the orthonormal transformation of \code{TGUW} in reverse order.
#' Details of the inverse TGUW transformation can be found in H. Maeng and P. Fryzlewicz (2019), Detecting linear trend changes and point anomalies in data sequences, preprint.
#'
#' @param ts.obj A list returned by \code{thresholding}.
#' @return
#' \item{ts.obj}{The modified ts.obj is the result of the inverse TGUW transformation and \code{ts.coeffs} now presents the estimated piecewise-linear signal of the data.}
#' @author Hyeyoung Maeng, \email{h.maeng@@lse.ac.uk}
#' @seealso \code{\link{trendsegment}}, \code{\link{TGUW}}, \code{\link{thresholding}}
#' @examples
#' x <- c(1:10, 0, rep(5,9))
#' n <- length(x)
#' x <- x + rnorm(n)
#' tguwfit <- TGUW(x)
#' th.const <- 1.3
#' lambda <- (stats::mad(diff(diff(x)))/sqrt(6)) * sqrt(2 * log(n)) * th.const
#' thrfit <- thresholding(ts.obj = tguwfit, lambda = lambda)
#' invfit <- invTGUW(ts.obj = thrfit)
#' invfit
#' @export





invTGUW <- function(ts.obj) {

  n <- ts.obj$n

  for (i in (n-2):1) {

    ######################################################################################################
    ##### STEP 1: Choose the corresponding high-pass filters (h1, h2, h3) and use it to obtain M.inv #####
    ######################################################################################################
    ### dim(ts.obj$high.pass) = (3)x(n-2)
    ### Choose the column of high-pass filter from the final one, and apply the "sym.orthmatrix" function to have M matrix
    M.inv <- t(orthmatrix(ts.obj$merging.hist[2,,i]))


    #########################################################################################################################################
    ##### STEP 2: Select the index and the corresponding (min(detail), x) and recover 3 values of x by using M.inv obtained in STEP 1.  #####
    #########################################################################################################################################
    ### Choose the index from the final one
    ind <- ts.obj$merging.hist[1,,i]

    ### Choose the corresponding combination of detail and x value
    ts.obj$merging.hist[3,2,i] <- ts.obj$ts.coeffs[ind[1]]
    ts.obj$merging.hist[3,3,i] <- ts.obj$ts.coeffs[ind[2]]

    tmp <- c(ts.obj$merging.hist[3,1,i], ts.obj$merging.hist[3,2,i], ts.obj$merging.hist[3,3,i])

    ### Obtain "3" values of x from the combination of detail coefficient and "2" values of x, tmp=c( min(details), x ),
    ### by multiplying M.inv (inverse of transformation M)
    rcstr.tmp <- M.inv %*% tmp

    ########################################################################################################################################
    ##### STEP 3: Update res with old x from previous res (which has indices unrelated to this step) and newly obtained 3 values of X  #####
    ########################################################################################################################################
    ### The old "res" part with index [(ind+2):(n-1)], which is unrelated part of this step, is used as a replacement of new res with index [(ind+3):n]

    ### newly obtained three values (rcstr.tmp) will be used in new res
    ts.obj$ts.coeffs[ind[1]] <- rcstr.tmp[1,1]
    ts.obj$ts.coeffs[ind[2]] <- rcstr.tmp[2,1]
    ts.obj$ts.coeffs[ind[3]] <- rcstr.tmp[3,1]

  }

  return(ts.obj)

}
