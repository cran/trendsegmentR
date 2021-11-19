#' Tail-Greedy Unbalanced Wavelet (TGUW) transformation of a vector
#'
#' Performs the bottom-up unbalanced wavelet decomposition. This function is used inside \code{\link{trendsegment}}.
#' Details of the TGUW transformation can be found in H. Maeng and P. Fryzlewicz (2021), Detecting linear trend changes in data sequences, preprint.
#'
#' @param x An input vector to be decomposed.
#' @param p Proportion of all possible remaining merges which specifies the number of merges allowed in a single pass over the data. The default is 0.01.
#' @return A list with the followings:
#' \item{x}{The original input vector \code{x}.}
#' \item{n}{The length of \code{x}.}
#' \item{twotogether}{A vector indicating locations of the detail coefficients returned by Type 3 merges (merging two sets of paired smooth coefficients). This is used in \code{\link{thresholding}} to apply the "two together" rule which makes both detail coefficients (paired by a Type 3 merge) survived if at least one of their size is over threshold.}
#' \item{merging.hist}{An array of dimension 4 by 3 by \code{n}-2 which has the full record of the \code{n}-2 merges in the TGUW transformation. Each matrix contains the information of each merge. The first row shows the indices of merged smooth coefficients in increasing order and the second row gives the value of detail filter coefficients which is the weight vector for computing the corresponding detail coefficient. The third row shows the (detail coefficient, first smooth coefficient, second smooth coefficient) obtained by an orthonormal transform. The fourth row gives the balancedness of merging. If it is Type 1 merging (three initial smooth coefficients) then the fourth row is always (1/3, 1/3, 1/3). In Type 2 and Type 3 mergings, the values depend on the ratio of the length of the left and right wings to the entire merged region and only first two components of the fourth row are filled with the corresponding ratios (sum to 1) but the third one is left as NA.}
#' \item{ts.coeffs}{The transformed \code{x} by the TGUW transformation.}
#' @author Hyeyoung Maeng \email{h.maeng4@@lancaster.ac.uk}, Piotr Fryzlewicz \email{p.fryzlewicz@@lse.ac.uk}
#' @seealso \code{\link{trendsegment}}, \code{\link{thresholding}}, \code{\link{invTGUW}}
#' @examples
#' x <- c(1:10, rep(5,9))
#' n <- length(x)
#' x <- x + rnorm(n)
#' tguwfit <- TGUW(x)
#' tguwfit
#' @export


TGUW <- function(x, p = .01) {

  n <- length(x)
  noe <- n-2 # to be updated for each scale j

  weights.const <- rep(1, n)
  weights.lin <- 1:n
  idx <- 1:n
  paired <- c()

  edges <- matrix(0, noe, 3) # to be updated for each scale j

  edges[,1] <- 1:(n-2)
  edges[,2] <- 2:(n-1)
  edges[,3] <- 3:n

  edges <- cbind(edges, 0)

  merging.hist <- array(0, dim=c(4, 3, n-2))

  ts.coeffs <- as.numeric(x)

  steps.left <- n-1 # to be updated for each scale j
  current.step <- 0 # to be updated for each scale j

  twotogether <- c()

  while (dim(edges)[1]) {


    max.current.steps <- ceiling(p * steps.left)
    removable.nodes <- rep(1, max(idx))

    ###############################################################################
    ##### STEP 1: obtain a vector of high-pass filter (h) and compute details #####
    ###############################################################################

    if(any(edges[, 4]>0)){

      ### 1) 2 sets of pairs
      pr <- matrix(which(edges[,4]!=0), nrow=2)

      ## i) first edges
      cd1 <- computeDET(edges=edges, edgerow=pr[1,], weights.const=weights.const, weights.lin=weights.lin, ts.coeffs=ts.coeffs)
      detcoef <- cd1$detcoef
      p1.details <- cd1$det

      ## ii) second edges
      M0 <- apply(detcoef, 2, orthmatrix)
      upd.wc <- cbind(rowSums(cd1$wc*t(M0[c(2,5,8),])), rowSums(cd1$wc*t(M0[c(3,6,9),])), weights.const[edges[pr[2,],3]])
      upd.wl <- cbind(rowSums(cd1$wl*t(M0[c(2,5,8),])), rowSums(cd1$wl*t(M0[c(3,6,9),])), weights.lin[edges[pr[2,],3]])
      upd.ts.coeffs <- cbind(rowSums(cd1$tc*t(M0[c(2,5,8),])), rowSums(cd1$tc*t(M0[c(3,6,9),])), ts.coeffs[edges[pr[2,],3]])
      p2.details <- colSums(apply(cbind(upd.wc, upd.wl), 1, filter)*t(upd.ts.coeffs))

      p.detail <- rep(apply(cbind(abs(p1.details), abs(p2.details)), 1, max), each=2)

      if(length(pr)!=dim(edges)[1]){
        ### 2) non-paired
        cd3 <- computeDET(edges=edges, edgerow=c(1:dim(edges)[1])[-c(pr)], weights.const=weights.const, weights.lin=weights.lin, ts.coeffs=ts.coeffs)
        np.details <- cd3$det
        details <- c(p.detail, np.details)
      } else{
        details <- p.detail
      }

    } else{
      ### compute details for every combination of (x_i, x_{i+1}, x_{i+2}) with high-pass filter0
      cd <- computeDET(edges=edges, edgerow=c(1:dim(edges)[1]), weights.const=weights.const, weights.lin=weights.lin, ts.coeffs=ts.coeffs)
      details <- cd$det

    }

    ######################################################################################
    ##### STEP 2: sort detail coefs and indices / find the index of X to be removed  #####
    ######################################################################################

    ord.det <- order(abs(details)) # it gives the index of "abs(details)" with its ascending order
    cand <- rbind(ord.det, edges[ord.det,4])

    ### The first run
    eitr <- 1 # edge.indices.2b.removed
    tei <- 1 # traverse.edges.index
    if(cand[2,1] > 0){
      removable.nodes[edges[ord.det[1:2],1]] <- removable.nodes[edges[ord.det[1:2],2]] <- removable.nodes[edges[ord.det[1:2],3]] <- 0
      tei <- tei + 1
      eitr <- c(eitr, tei)
    } else{
      removable.nodes[edges[ord.det[1],1]] <- removable.nodes[edges[ord.det[1],2]] <- removable.nodes[edges[ord.det[1],3]] <- 0
    }

    ### More
    while  ( ( length(eitr) < max.current.steps ) & ( tei < noe ) ) {

      tei <- tei + 1

      if(cand[2, tei] > 0){
        if(sum(removable.nodes[edges[ord.det[tei:(tei+1)],1]])==2 & sum(removable.nodes[edges[ord.det[tei:(tei+1)],2]])==2 & sum(removable.nodes[edges[ord.det[tei:(tei+1)],3]])==2){
          removable.nodes[edges[ord.det[tei:(tei+1)],1]] <- removable.nodes[edges[ord.det[tei:(tei+1)],2]] <- removable.nodes[edges[ord.det[tei:(tei+1)],3]] <- 0
          eitr <- c(eitr, tei, tei+1)
        }
        tei <- tei + 1
      } else{
        if(removable.nodes[edges[ord.det[tei],1]] & removable.nodes[edges[ord.det[tei],2]] & removable.nodes[edges[ord.det[tei],3]]){
          eitr <- c(eitr, tei)
          removable.nodes[edges[ord.det[tei],1]] <- removable.nodes[edges[ord.det[tei],2]] <- removable.nodes[edges[ord.det[tei],3]] <- 0
        }
      }
    }

    # get the first indices of two consecutive points which will be combined
    details.min.ind <- ord.det[eitr]

    # number of pairs to be combined
    no.of.current.steps <- length(eitr)

    ##########################################################################################
    ##### STEP 3: obtain the low-pass filters and smoothed (X, weight.const, weight.lin) #####
    ##########################################################################################
    ##### STEP 4: UPDATING the corresponding old (X, weight.const, weight.lin)           #####
    ##### with smoothed (X, weight.const, weight.lin)                                    #####
    ##########################################################################################
    ### create low-pass filter (l1, l2) by using high-pass filter h
    ### t(M)%*%M is an identity matrix (which means M is a orthonormal matrix)

    ee <- matrix(edges[details.min.ind,], no.of.current.steps, 4)
    #twotogether <- c(twotogether, c(ee[,4]))
    idx0 <- idx # keep original idx

    if(sum(ee[,4]>0)==0){

      twotogether <- c(twotogether, c(ee[,4]))

      ee <- ee[, -4, drop=F]

      udt <- updating(ee=ee, weights.const=weights.const, weights.lin=weights.lin, ts.coeffs=ts.coeffs, idx=idx)
      weights.const <- udt$weights.const
      weights.lin <- udt$weights.lin
      ts.coeffs <- udt$ts.coeffs
      idx <- udt$idx

      merging.hist[1,,(current.step+1):(current.step+no.of.current.steps)] <- t(ee)
      merging.hist[2,,(current.step+1):(current.step+no.of.current.steps)] <- udt$h
      merging.hist[3,,(current.step+1):(current.step+no.of.current.steps)] <- t(udt$tc1)
      merging.hist[4,,(current.step+1):(current.step+no.of.current.steps)] <- balance.np(paired=paired, ee=ee, idx=idx0, no.of.current.steps=no.of.current.steps, n=n)

    } else{

      twotogether <- c(twotogether, c(ee[which(ee[,4]!=0), 4]))

      ### 1) 2 sets of pairs
      pr <- matrix(which(ee[,4]!=0), nrow=2)
      ee[pr[2,], 1:2] <- ee[pr[1,], 1:2]

      #ee <- ee[, -4, drop=F]


      ee.p1 <- ee[pr[1,], -4, drop=F]
      ee.p2 <- ee[pr[2,], -4, drop=F]

      ######################################## i) first paired edges
      udt <- updating(ee=ee.p1, weights.const=weights.const, weights.lin=weights.lin, ts.coeffs=ts.coeffs, idx=idx)
      weights.const <- udt$weights.const
      weights.lin <- udt$weights.lin
      ts.coeffs <- udt$ts.coeffs
      idx <- udt$idx

      merge1 <- seq((current.step+1), (current.step+length(pr)), by=2)
      merge2 <- merge1 + 1

      #merging.hist[1,,(current.step+1):(current.step+dim(ee.p1)[1])] <- t(ee.p1)
      #merging.hist[2,,(current.step+1):(current.step+dim(ee.p1)[1])] <- udt$h
      #merging.hist[3,,(current.step+1):(current.step+dim(ee.p1)[1])] <- t(udt$tc1)
      #merging.hist[4,,(current.step+1):(current.step+dim(ee.p1)[1])] <- balance.p(pr=pr, ee.p1=ee.p1, idx=idx0, ee.p2=ee.p2, n=n)
      merging.hist[1,,merge1] <- t(ee.p1)
      merging.hist[2,,merge1] <- udt$h
      merging.hist[3,,merge1] <- t(udt$tc1)
      merging.hist[4,,merge1] <- balance.p(pr=pr, ee.p1=ee.p1, idx=idx0, ee.p2=ee.p2, n=n)

      ######################################## ii) second paired edges

      udt <- updating(ee=ee.p2, weights.const=weights.const, weights.lin=weights.lin, ts.coeffs=ts.coeffs, idx=idx)
      weights.const <- udt$weights.const
      weights.lin <- udt$weights.lin
      ts.coeffs <- udt$ts.coeffs
      idx <- udt$idx

      #merging.hist[1,,(current.step+dim(ee.p1)[1]+1):(current.step+length(pr))] <- t(ee.p2)
      #merging.hist[2,,(current.step+dim(ee.p1)[1]+1):(current.step+length(pr))] <- udt$h
      #merging.hist[3,,(current.step+dim(ee.p1)[1]+1):(current.step+length(pr))] <- t(udt$tc1)
      #merging.hist[4,,(current.step+dim(ee.p1)[1]+1):(current.step+length(pr))] <- balance.p(pr=pr, ee.p1=ee.p1, idx=idx0, ee.p2=ee.p2, n=n)
      merging.hist[1,,merge2] <- t(ee.p2)
      merging.hist[2,,merge2] <- udt$h
      merging.hist[3,,merge2] <- t(udt$tc1)
      merging.hist[4,,merge2] <- balance.p(pr=pr, ee.p1=ee.p1, idx=idx0, ee.p2=ee.p2, n=n)

      if(length(pr)!=dim(ee)[1]){
        ######################################## iii) non paired
        twotogether <- c(twotogether, c(ee[which(ee[,4]==0), 4]))
        ee.np <- ee[-c(pr), -4, drop=F]

        udt <- updating(ee=ee.np, weights.const=weights.const, weights.lin=weights.lin, ts.coeffs=ts.coeffs, idx=idx)
        weights.const <- udt$weights.const
        weights.lin <- udt$weights.lin
        ts.coeffs <- udt$ts.coeffs
        idx <- udt$idx

        merging.hist[1,,(current.step+length(pr)+1):(current.step+no.of.current.steps)] <- t(ee.np)
        merging.hist[2,,(current.step+length(pr)+1):(current.step+no.of.current.steps)] <- udt$h
        merging.hist[3,,(current.step+length(pr)+1):(current.step+no.of.current.steps)] <- t(udt$tc1)
        merging.hist[4,,(current.step+length(pr)+1):(current.step+no.of.current.steps)] <- balance.np(paired=paired, ee=ee.np, idx=idx0, no.of.current.steps=no.of.current.steps-length(pr), n=n)

      }

      ee <- ee[, -4, drop=F]

    }

    ############################################
    ##### STEP 5: Updating other variables #####
    ############################################

    paired <- sort(unique(c(paired, c(ee[,1:2]))))
    if(length(which(!is.na(match(paired, c(ee[,3])))))>0){
      paired <- sort(paired[-c(which(!is.na(match(paired, c(ee[,3])))))])
    }

    edges <- matrix(0, length(idx)-2, 3) # to be updated for each scale j
    edges[,1] <- idx[1:(length(idx)-2)]
    edges[,2] <- idx[2:(length(idx)-1)]
    edges[,3] <- idx[3:(length(idx))]

    ### the edges rows matched with edges.0 are survived
    matchpair <- ifelse(is.na(matrix(match(edges, paired), ncol=3)), NA, 1)
    rs <- which(rowSums(matchpair, na.rm=T)==3)
    if(length(rs)>0){
      edges <- as.matrix(rbind(edges[rs,], edges[-rs,]))
      matchpair <- rbind(matchpair[rs,], matchpair[-rs,])
      edges <- as.matrix(cbind(edges, c(rep(1:(length(rs)/2), each=2), rep(0, dim(edges)[1]-length(rs)))))
    } else{
      edges <- as.matrix(cbind(edges, rep(0, dim(edges)[1])))
    }

    ### finding the other survived edges rows
    removed <- c(which(rowSums(matchpair, na.rm=T)==1), which(matchpair[,1]==1 & is.na(matchpair[,2]) & matchpair[,3]==1))
    if(length(removed)>0){
      edges <- edges[-removed,,drop=F]
    }

    # define noe again as the possible combinations of edges
    noe <- dim(edges)[1]
    steps.left <- steps.left - no.of.current.steps
    current.step <- current.step + no.of.current.steps

    ### give zero to the index of paring for last remained edge(1,2,x)
    if(noe==1 & dim(edges)[1]==1){
      edges[,4] <- 0
    }

  }

  return(list(x=x, n = n, twotogether=twotogether, merging.hist=merging.hist, ts.coeffs=ts.coeffs))

}
