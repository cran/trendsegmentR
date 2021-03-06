% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/orthmatrix.R
\name{orthmatrix}
\alias{orthmatrix}
\title{Finding an orthonormal matrix based on a detail filter vector}
\usage{
orthmatrix(d)
}
\arguments{
\item{d}{A detail filter returned by \code{\link{filter}} which has a form of a vector with length 3.}
}
\value{
\item{M}{The orthonormal matrix with dimension 3 by 3 which is used in \code{\link{TGUW}} and \code{\link{invTGUW}}.}
}
\description{
This function is used inside \code{\link{TGUW}} function but is typically not called directly by the user. This gives an orthonormal matrix with dimension 3 by 3 by computing two low filter vectors based on a given detail filter vector. The returned orthonormal matrix is firstly used in the orthonormal transformation of \code{\link{TGUW}} when updating three neighbouring smooth coefficients into one detail and two smooth coefficients, and its inverse matrix is used in the inverse TGUW transformation (\code{\link{invTGUW}}).
}
\details{
For details, see H. Maeng and P. Fryzlewicz (2019), Detecting linear trend changes and point anomalies in data sequences, preprint.
}
\examples{
x <- c(rep(1, 3), 1:3)
df <- filter(x) # detail filter
orthmatrix(df)
}
\seealso{
\code{\link{TGUW}}, \code{\link{invTGUW}}, \code{\link{filter}}
}
\author{
Hyeyoung Maeng, \email{h.maeng@lse.ac.uk}
}
