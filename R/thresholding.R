#' Noise reduction from the sequence of detail coefficients returned by the Tail-Greedy Unbalanced Wavelet (TGUW) transformation
#'
#' This function is used inside \code{\link{trendsegment}} and performs the thresholding of the detail coefficients returned by the Tail-Greedy Unbalanced Wavelet (TGUW) transformation.
#' The denoising is achieved by a prespecified threshold in a "connected" way in that it prunes the branches if and only if
#' the detail coefficient itself and all of its children coefficients are below some thresholds in its size. Also, the
#' "two together" rule is applied to any paired detail coefficients returned by Type 3
#' merging (merging two sets of paired smooth coefficients) in the sense that both detail coefficients should be survived if at least one of their size is over threshold. For details, see H. Maeng and P. Fryzlewicz (2021), Detecting linear trend changes in data sequences, preprint.
#'
#' @param ts.obj A list returned by \code{TGUW}.
#' @param lambda The magnitude of the threshold. It has a form of \code{sigma * th.const * sqrt(2 * log(n))} where \code{n} is the length of input data \code{x}, the default of \code{th.const} is 1.3 and the \code{sigma} can be estimated by Median Absolute Deviation (MAD) method under the Gaussian assumption for noise.
#' @param minsegL The minimum segment length of estimated signal returned by \code{trendsegment}.
#' @param bal The minimum ratio of the length of the shorter region to the length of the entire merging region especially when the merges of Type 2 (merging one initial and a paired smooth coefficient) or of Type 3 (merging two sets of (paired) smooth coefficients) are performed. Only triplets which satisfy this balancedness condition survives in denoising. Point anomalies can be detected only if \code{bal < 1/n} and \code{minsegL = 1}. The default is set to 0.
#' @param connected If connected=TRUE, the thresholding puts the connected rule above the \code{minsegL}, otherwise it makes keeping the \code{minsegL} a priority.
#' @return
#' \item{ts.obj}{The modified ts.obj containing zero detail coefficients in the \code{merging.hist} if not survived from thresholding.}
#' @author Hyeyoung Maeng \email{h.maeng4@@lancaster.ac.uk}, Piotr Fryzlewicz \email{p.fryzlewicz@@lse.ac.uk}
#' @seealso \code{\link{trendsegment}}, \code{\link{TGUW}}, \code{\link{invTGUW}}
#' @examples
#' x <- c(1:10, rep(5,9))
#' n <- length(x)
#' x <- x + rnorm(n)
#' tguwfit <- TGUW(x)
#' th.const <- 1.3
#' lambda <- (stats::mad(diff(diff(x)))/sqrt(6)) * sqrt(2 * log(n)) * th.const
#' thrfit <- thresholding(ts.obj = tguwfit, lambda = lambda, minsegL = 5, bal = 0, connected = FALSE)
#' thrfit
#' @export

thresholding <- function(ts.obj, lambda, minsegL, bal = 0, connected = FALSE) {

  details <- ts.obj$merging.hist[3,1,]
  twotogether <- ts.obj$twotogether
  n <- ts.obj$n

  ##### 1) make the tree shape connected
  protected <- rep(0, n)

  for (i in 1:(n-2)) {

    ### !protected[ts.obj$merging.hist[1,1:2,i]] == TRUE
    ### <-> protected[the indices of the regions merged] == 0
    ### <-> Those indices have not been survived and their children detail coefficients are estimated by zero in previous stages
    ### <-> if the children detail coefficients are not zero then we skip their parents detail coefficients and just survive it
    if (!protected[ts.obj$merging.hist[1,1,i]] & !protected[ts.obj$merging.hist[1,2,i]] &
        !protected[ts.obj$merging.hist[1,3,i]])

      ts.obj$merging.hist[3,1,i] <- ts.obj$merging.hist[3,1,i] * (
        ######################### TRUE==1, FALSE==0 #########################
        ### condition 1 : whether |detail coef| is greather than lambda
        (abs(ts.obj$merging.hist[3,1,i]) > lambda) &
          (ts.obj$merging.hist[4,1,i] > bal) &
          (ts.obj$merging.hist[4,2,i] > bal) &
          (ifelse(ts.obj$merging.hist[4,3,i] < 1, (1 >= minsegL), (min(ts.obj$merging.hist[4,1:2,i]*ts.obj$merging.hist[4,3,i]) >= minsegL)))
      )

    if (abs(ts.obj$merging.hist[3,1,i]) > 0) protected[ts.obj$merging.hist[1,1,i]] <- protected[ts.obj$merging.hist[1,2,i]] <- 1
  }

  ##### 2) make the pairs both survive if at least one is survived
  paired <- matrix(which(ts.obj$twotogether!=0), nrow=2)

  if(length(paired)>0){
    for(i in 1:ncol(paired)){
      overzero <- is.element(paired[,i], which(abs(ts.obj$merging.hist[3,1,])>0))
      zero <- is.element(paired[,i], which(abs(ts.obj$merging.hist[3,1,])==0))
      if(sum(overzero)==1 & sum(zero)==1){
        ts.obj$merging.hist[3,1,paired[which(overzero==F),i]] <- details[paired[which(overzero==F),i]]
      }
    }
  }

  if(connected == FALSE){

    svved <- which(ts.obj$merging.hist[3,1,]!=0) # survived indices
    should.remove <- c()
    for(i in svved){
      if(any(ts.obj$merging.hist[4, 1:2, i]*ts.obj$merging.hist[4, 3, i] < minsegL)){
        should.remove <- c(should.remove, i)
      }
    }
    if(length(should.remove)>0){
      ts.obj$merging.hist[3,1,should.remove] <- 0
    }

  }

  ts.obj$lambda <- lambda

  return(ts.obj)

}

