#' Finding an orthonormal matrix based on a detail filter vector
#'
#' This function is used inside \code{\link{TGUW}} function but is typically not called directly by the user. This gives an orthonormal matrix with dimension 3 by 3 by computing two low filter vectors based on a given detail filter vector. The returned orthonormal matrix is firstly used in the orthonormal transformation of \code{\link{TGUW}} when updating three neighbouring smooth coefficients into one detail and two smooth coefficients, and its inverse matrix is used in the inverse TGUW transformation (\code{\link{invTGUW}}).
#'
#' For details, see H. Maeng and P. Fryzlewicz (2019), Detecting linear trend changes and point anomalies in data sequences, preprint.
#'
#' @param d A detail filter returned by \code{\link{filter}} which has a form of a vector with length 3.
#' @return
#' \item{M}{The orthonormal matrix with dimension 3 by 3 which is used in \code{\link{TGUW}} and \code{\link{invTGUW}}.}
#' @author Hyeyoung Maeng, \email{h.maeng@@lse.ac.uk}
#' @seealso \code{\link{TGUW}}, \code{\link{invTGUW}}, \code{\link{filter}}
#' @examples
#' x <- c(rep(1, 3), 1:3)
#' df <- filter(x) # detail filter
#' orthmatrix(df)
#' @export



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
