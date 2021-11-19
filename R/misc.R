# Finding a detail filter vector for computing a detail coefficient
#
# This function is used inside \code{\link{TGUW}} function but is typically not called directly by the user. This returns a vector of detail filter, a weight vector of the triplet of smooth coefficients where the corresponding detail coefficient is obtained as a weighted sum of three smooth coefficients.
#
# The detail filter is obtained in such a way to produce zero detail coefficient only when the corresponding (raw) observations in merging regions have a perfect linear trend, as the detail coefficient itself represents the extent of non-linearity in the corresponding region of data. This implies that the smaller the size of the detail coefficient, the closer the alignment of the corresponding data section with linearity.
#
# @param a A vector with length 6. First three elements contain a triplet selected from the weight vector of constancy and the last three elements correspond to a selected triplet of the weight vector of linearity.
# @return
# \item{df}{The detail filter vector of length 3 which is used as a weight vector for the corresponding triplet of smooth coefficients.}
# @author Hyeyoung Maeng \email{h.maeng4@@lancaster.ac.uk}, Piotr Fryzlewicz \email{p.fryzlewicz@@lse.ac.uk}
# @seealso \code{\link{TGUW}}
# @examples
# x <- c(rep(1, 3), 1:3)
# filter(x)
#
# y <- c(rep(1, 3), 4:6)
# filter(y) # same as filter(x)
# @export

filter <- function(a) {

  #	work out a 3-tap high-pass filter which annihilates vectors a and b

  w <- -sqrt( (a[5]*a[1] - a[4]*a[2])^2 / ((a[2]*a[6] - a[3]*a[5])^2 + (a[4]*a[3] - a[1]*a[6])^2 + (a[5]*a[1] - a[4]*a[2])^2))

  u <- w * (a[2]*a[6] - a[3]*a[5])/(a[5]*a[1] - a[4]*a[2])

  v <- w * (a[3]*a[4] - a[1]*a[6])/(a[5]*a[1] - a[4]*a[2])

  df <- c(u, v, w)

  if (any(is.na(df))) {
    z <- filter(a[6:1])
    df <- z[3:1]
  }
  return(df)

}





# Finding an orthonormal matrix based on a detail filter vector
#
# This function is used inside \code{\link{TGUW}} function but is typically not called directly by the user. This gives an orthonormal matrix with dimension 3 by 3 by computing two low filter vectors based on a given detail filter vector. The returned orthonormal matrix is firstly used in the orthonormal transformation of \code{\link{TGUW}} when updating three neighbouring smooth coefficients into one detail and two smooth coefficients, and its inverse matrix is used in the inverse TGUW transformation (\code{\link{invTGUW}}).
#
# @param d A detail filter returned by \code{\link{filter}} which has a form of a vector with length 3.
# @return
# \item{M}{The orthonormal matrix with dimension 3 by 3 which is used in \code{\link{TGUW}} and \code{\link{invTGUW}}.}
# @author Hyeyoung Maeng \email{h.maeng4@@lancaster.ac.uk}, Piotr Fryzlewicz \email{p.fryzlewicz@@lse.ac.uk}
# @seealso \code{\link{TGUW}}, \code{\link{invTGUW}}, \code{\link{filter}}
# @examples
# x <- c(rep(1, 3), 1:3)
# df <- filter(x) # detail filter
# orthmatrix(df)
# @export

orthmatrix <- function(d) {

  M <- matrix(0, 3, 3)

  ##### STEP 1 : Normalisaion of high-pass filter row #####
  M[1,] <- d
  M[1,] <- M[1,] / sqrt(sum(M[1,]^2))
  u <- M[1, 1]
  v <- M[1, 2]
  w <- M[1, 3]

  ##### STEP 2 : Choose the vector which gives zero for the inner product of the first one #####
  M[2,] <- c(1-u^2, -u*v, -u*w)
  M[3,] <- c(0, -w, v)

  ##### STEP 3 : Normalisation #####
  M[2,] <- M[2,] / sqrt(sum(M[2,]^2))
  M[3,] <- M[3,] / sqrt(sum(M[3,]^2))

  return(M)

}





# Finding the detected change-points from the object returned by the inverse TGUW transformation
# INTERNAL function in \code{\link{trendsegment}} to obtain the estimated change-points by trendsegment algorithm. This function is typically not called directly by the user.
#  @param ts.obj An object returned by \code{invTGUW}.
#  @return
#  \item{cp}{The estimated change-points}

finding.cp <- function(ts.obj){

  if(length(ts.obj$twotogether)==1){
    all.edges <- matrix(c(ts.obj$merging.hist[1,,], ts.obj$twotogether), ncol=1) #twotogether shows pairs
  } else{
    all.edges <- rbind(ts.obj$merging.hist[1,,], ts.obj$twotogether) #twotogether shows pairs
  }
  survived.edges <- all.edges[ ,which(abs(ts.obj$merging.hist[3,1,])>1e-10), drop=F]

  if(length(survived.edges)==0){
    cp <- c()
  } else if(length(survived.edges)>0 & dim(survived.edges)[2]>1){
    i <- 1
    cp <- c()
    while(i<dim(survived.edges)[2]){
      part.obj <- ts.obj$merging.hist[1,c(1,2),] # for case 3): to differentiate it from case 4)
      matched <- which(!is.na(match(data.frame(part.obj), data.frame(matrix(survived.edges[c(2:3),i], ncol=1))))) # for case 3): to differentiate it from case 4
      ### 1) (xx from a chunk, xx from another chunk)
      if((survived.edges[4,i]!=0) & diff(survived.edges[-4,i])[1]==1 & diff(survived.edges[-4,i+1])[1]==1){
        cp <- c(cp, survived.edges[c(1, 3), i])
        i <- i+2
      } else if((survived.edges[4,i]!=0) & diff(survived.edges[-4,i])[2]==1 & diff(survived.edges[-4,i+1])[1]==1){
        cp <- c(cp, survived.edges[c(1, 3), i+1])
        i <- i+2
      } else if((survived.edges[4,i]==0) & diff(survived.edges[-4,i])[1]==1 & diff(survived.edges[-4,i])[2]!=1){ ### 2) (xx from a chunk, x)
        cp <-c(cp, survived.edges[c(1,3), i])
        i <- i+1
      } else if((survived.edges[4,i]==0) & diff(survived.edges[-4,i])[1]==1 & diff(survived.edges[-4,i])[2]==1 & length(matched)>0){ ### 3) (x, xx from a chunk)
        cp <-c(cp, survived.edges[c(1,2), i])
        i <- i+1
      } else {
        cp <-c(cp, survived.edges[c(1,2,3), i]) ### 4) (x, x, x)
        i <- i+1
      }
    }
    cp <- unique(sort(cp))
    cp <- c(cp, ts.obj$n+1)
  } else if(length(survived.edges)>0 & dim(survived.edges)[2]==1 & survived.edges[4,1]==0 & diff(survived.edges[-4,1])[1]==1 & diff(survived.edges[-4,1])[2]!=1){
    cp <- c()
    cp <- c(cp, survived.edges[c(1,3), 1])
  } else if(length(survived.edges)>0 & dim(survived.edges)[2]==1 & survived.edges[4,1]==0 & diff(survived.edges[-4,1])[1]==1 & diff(survived.edges[-4,1])[2]==1){ ### 3) (x, xx from a chunk)
    cp <- c()
    cp <- c(cp, survived.edges[c(1,2), 1])
  } else {
    cp <- c()
  }

  ### last adjustment
  # if the number of observation is three and the detail is survived then c.p=c(2,3)
  if(ts.obj$n==3 & length(cp)>0 & dim(survived.edges)[2]==1){
    cp <- cp[-1]
    cp <- c(cp, ts.obj$n)
  } else {
    cp <- cp[which(cp<=ts.obj$n & cp>1)]
  }

  ### for comparing with other methods
  if(length(cp)>0){
    cp <- cp-1
  }
  return(cp)
}





# Computing the detail coefficient for all triplets at each scale
# INTERNAL function in \code{\link{TGUW}} to obtain the sequence of detail coefficients for all candidate triplets of smooth coefficients. This function is typically not called directly by the user.
# @param edges A matrix with 4 columns. Each row corresponds to indices of each triplet and the row length is the number of allowed merges.
# @param edgerow A vector containing the selected row indices of \code{edges}. Detail coefficients are computed only for these chosen rows.
# @param weights.const The up-to-date weight vector of constancy.
# @param weights.lin The up-to-date weight vector of linearity.
# @param ts.coeffs The up-to-date vector of transformed \code{x}.
# @return A list with the following:
# \item{detcoef}{The matrix containing each detail filter vector in each row.}
# \item{det}{The vector of detail coefficients for selected triplets in \code{edgerow}.}
# \item{wc}{The matrix containing selected triplets of the weight vector of constancy.}
# \item{wl}{The matrix containing selected triplets of the weight vector of linearity.}
# \item{tc}{The matrix containing selected triplets of transformed \code{x}.}

computeDET <- function(edges = edges, edgerow = edgerow, weights.const = weights.const, weights.lin = weights.lin, ts.coeffs = ts.coeffs){

  sub.wc <- cbind(weights.const[edges[edgerow,1]], weights.const[edges[edgerow,2]], weights.const[edges[edgerow,3]])
  sub.wl <- cbind(weights.lin[edges[edgerow,1]], weights.lin[edges[edgerow,2]], weights.lin[edges[edgerow,3]])
  sub.tc <- cbind(ts.coeffs[edges[edgerow,1]], ts.coeffs[edges[edgerow,2]], ts.coeffs[edges[edgerow,3]])
  detcoef <- apply(cbind(sub.wc, sub.wl), 1, filter)
  details <- colSums(detcoef*t(sub.tc))
  return(list(detcoef=detcoef, det=details, wc=sub.wc, wl=sub.wl, tc=sub.tc))
}








# Updating two weight vectors (of constancy and of linearity) and data sequence for all chosen triplets by orthonormal transforms
# INTERNAL function in \code{\link{TGUW}} to update weight vectors and data sequence for all selected triplets by orthonormal transform. This function is typically not called directly by the user.
# @param ee A matrix with 4 columns. Each row corresponds to indices of selected triplet and the row length is the number of merges achieved at a certain scale.
# @param weights.const The up-to-date weight vector of constancy.
# @param weights.lin The up-to-date weight vector of linearity.
# @param ts.coeffs The up-to-date transformed input \code{x}.
# @param idx The survived indices of \code{x} which correspond to the indices of smooth coefficients in \code{x}.
# @return A list with the following:
# \item{weights.const}{The updated weight vector of constancy through the orthonormal transformations performed for all selected triplets.}
# \item{weights.lin}{The updated weight vector of linearity through the orthonormal transformations performed for all selected triplets.}
# \item{ts.coeffs}{The updated data sequence \code{x} through the orthonormal transformations performed for all selected triplets.}
# \item{idx}{The updated survived indices of \code{x} through the orthonormal transformations performed for all selected triplets.}
# \item{h}{The matrix containing the detail filter vector for all selected triplets. This is recorded in the second row of \code{merging.hist}.}
# \item{tc1}{The vector containing the transformed triplets of smooth coefficients at a certain scale. This is recorded in the third row of \code{merging.hist}.}

updating <- function(ee = ee, weights.const = weights.const, weights.lin = weights.lin, ts.coeffs = ts.coeffs, idx = idx){

  wc0 <- cbind(weights.const[ee[,1]], weights.const[ee[,2]], weights.const[ee[,3]])
  wl0 <- cbind(weights.lin[ee[,1]], weights.lin[ee[,2]], weights.lin[ee[,3]])
  tc0 <- cbind(ts.coeffs[ee[,1]], ts.coeffs[ee[,2]], ts.coeffs[ee[,3]])

  h <- apply(cbind(wc0,wl0), 1, filter)
  M0 <- apply(h, 2, orthmatrix)

  wc1 <- cbind(rowSums(wc0*t(M0[c(1,4,7),])), rowSums(wc0*t(M0[c(2,5,8),])), rowSums(wc0*t(M0[c(3,6,9),])))
  wl1 <- cbind(rowSums(wl0*t(M0[c(1,4,7),])), rowSums(wl0*t(M0[c(2,5,8),])), rowSums(wl0*t(M0[c(3,6,9),])))
  tc1 <- cbind(rowSums(tc0*t(M0[c(1,4,7),])), rowSums(tc0*t(M0[c(2,5,8),])), rowSums(tc0*t(M0[c(3,6,9),])))

  eating.up0 <- ee[,1]
  eating.up1 <- ee[,2]
  eaten.up <- ee[,3]
  idx <- idx[-which(is.na(match(idx, c(eaten.up)))==F)]

  ### 1) updating X
  ts.coeffs[eating.up0] <- c(tc1[,2])
  ts.coeffs[eating.up1] <- c(tc1[,3])
  ts.coeffs[eaten.up] <- c(tc1[,1])

  ### 2) updating weight.const
  weights.const[eaten.up] <- c(wc1[,1])
  weights.const[eating.up0] <- c(wc1[,2])
  weights.const[eating.up1] <- c(wc1[,3])

  ### 3) updating weight.lin
  weights.lin[eaten.up] <- c(wl1[,1])
  weights.lin[eating.up0] <- c(wl1[,2])
  weights.lin[eating.up1] <- c(wl1[,3])

  return(list(weights.const=weights.const, weights.lin=weights.lin, ts.coeffs=ts.coeffs, idx=idx, h=h, tc1=tc1))
}







# Computing the balancedness for Type 1 or Type 2 merges
# INTERNAL function in \code{\link{TGUW}} to compute the balancedness of each merge when the merge is either Type 1 (merging three initial smooth coefficients) or Type 2 (merging one initial and a paired smooth coefficient). This function is typically not called directly by the user.
# @param paired The vector containing all pairs of indices which are under "two together" rule.
# @param ee A matrix with 4 columns. Each row corresponds to indices of selected triplet and the row length is the number of merges done at a certain scale.
# @param idx The survived indices of \code{x}.
# @param no.of.current.steps The number of merges performed at a certain scale.
# @param n The length of input data \code{x}.
# @return
# \item{blnc}{The matrix containing balancedness of each merge in each column.}

balance.np <- function(paired=paired, ee=ee, idx, no.of.current.steps=no.of.current.steps, n=n){
  prd <- !is.na(matrix(match(ee, paired), ncol=3))
  blnc <- matrix(NA, nrow=3, ncol=no.of.current.steps)
  firsttwo <- which(rowSums(prd[,1:2, drop=F])==2)
  lasttwo <- which(rowSums(prd[,2:3, drop=F])==2)

  if(length(c(firsttwo, lasttwo))>0){
    nopair <- c(1:dim(ee)[1])[-c(firsttwo, lasttwo)]
  } else{
    nopair <- c(1:dim(ee)[1])
  }

  if(length(firsttwo)>0){
    prtn <- ee[firsttwo,3] - ee[firsttwo,1]
    #blnc[1:2, firsttwo] <- rbind(prtn/(prtn+1), 1/(prtn+1))
    blnc[1:3, firsttwo] <- rbind(prtn/(prtn+1), 1/(prtn+1), (prtn+1))
  }
  if(length(lasttwo)>0){
    prtn <- idx[match(ee[lasttwo, 3], idx)+1] - ee[lasttwo, 2]
    if(sum(is.na(prtn))>0){
      prtn[which(is.na(prtn))] <- n - ee[which(is.na(prtn)), 2] + 1
    }
    #blnc[1:2, lasttwo] <- rbind(1/(prtn+1), prtn/(prtn+1))
    blnc[1:3, lasttwo] <- rbind(1/(prtn+1), prtn/(prtn+1), (prtn+1))
  }
  if(length(nopair)>0){
    blnc[1:3, nopair] <- 1/3
  }
  return(blnc)
}







# Computing the balancedness for Type 3 merges
# INTERNAL function in \code{\link{TGUW}} to compute the balancedness of each merge when the merge is categorised into Type 3 (merging two sets of (paired) smooth coefficients). This function is typically not called directly by the user.
# @param pr The matrix indicating indices of Type 3 merges.
# @param ee.p1 The matrix containing the indices of triplets of first merges from all pairs.
# @param idx The survived indices of \code{x}.
# @param ee.p2 The matrix containing the indices of triplets of second merges from all pairs.
# @param n The length of input data \code{x}.
# @return
# \item{blnc}{The matrix containing balancedness of each merge in each column.}

balance.p <- function(pr=pr, ee.p1=ee.p1, idx, ee.p2=ee.p2, n=n){

  blnc <- matrix(NA, nrow=3, ncol=dim(pr)[2])

  c1 <- ee.p1[,3] - ee.p1[,1]
  c2 <- idx[match(ee.p2[, 3], idx)+1] - ee.p1[, 3]
  if(sum(is.na(c2))>0){
    c2[which(is.na(c2))] <- n - ee.p1[which(is.na(c2)), 3] + 1
  }

  #blnc[1:2,] <-  t(matrix(c(c1/(c1+c2), c2/(c1+c2)), nrow=2))
  #blnc[1:2,] <-  t(cbind(c1/(c1+c2), c2/(c1+c2)))
  blnc[1:3,] <-  rbind(c1/(c1+c2), c2/(c1+c2), (c1+c2))
  return(blnc)
}

