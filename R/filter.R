#' Finding a detail filter vector for computing a detail coefficient
#'
#' This function is used inside \code{\link{TGUW}} function but is typically not called directly by the user. This returns a vector of detail filter, a weight vector of the triplet of smooth coefficients where the corresponding detail coefficient is obtained as a weighted sum of three smooth coefficients.
#'
#' The detail filter is obtained in such a way to produce zero detail coefficient only when the corresponding (raw) observations in merging regions have a perfect linear trend, as the detail coefficient itself represents the extent of non-linearity in the corresponding region of data. This implies that the smaller the size of the detail coefficient, the closer the alignment of the corresponding data section with linearity. For details, see H. Maeng and P. Fryzlewicz (2019), Detecting linear trend changes and point anomalies in data sequences, preprint.
#'
#' @param a A vector with length 6. First three elements contain a triplet selected from the weight vector of constancy and the last three elements correspond to a selected triplet of the weight vector of linearity.
#' @return
#' \item{df}{The detail filter vector of length 3 which is used as a weight vector for the corresponding triplet of smooth coefficients.}
#' @author Hyeyoung Maeng, \email{h.maeng@@lse.ac.uk}
#' @seealso \code{\link{TGUW}}
#' @examples
#' x <- c(rep(1, 3), 1:3)
#' filter(x)
#'
#' y <- c(rep(1, 3), 4:6)
#' filter(y) # same as filter(x)
#' @export





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
