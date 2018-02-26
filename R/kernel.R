#' Epanechnikov Kernel
#'
#' @param datavec  a numeric vector
#' @param Bndwdth  a bandwidth
#' @export
#' @return kernel function result
#' @examples
#' Kh.Ep(2:7,5)
Kh.Ep <- function(datavec, Bndwdth) {
  KH <- numeric(length(datavec))
  xh<- datavec/Bndwdth
  Ind<- ( abs(xh) >=1 )
  KH[Ind] <- 0
  KH[!Ind] <- 0.75*(1-xh[!Ind]^2)/Bndwdth
  KH
}


#' Biweight kernel
#'
#' @param datavec   a numeric vector
#' @param Bndwdth   a bandwidth of the kernel
#' @export
#' @return  kernel function result
#' @examples
#' # same usage as Kh.Ep
Kh.Bw <- function(datavec, Bndwdth) {
  KH <- numeric(length(datavec))
  xh<- datavec/Bndwdth
  Ind<- ( abs(xh) >=1 )
  KH[Ind] <- 0
  KH[!Ind] <- (15/16)*(1-xh[!Ind]^2)^2/Bndwdth
  KH
}



#' Normal kernel
#'
#' @param datavec   a numeric vector
#' @param Bndwdth   a bandwidth of the kernel
#'
#' @importFrom stats dnorm
#' @export
#' @return  kernel function result
#' @examples
#' Kh.Nm(2:7,5)
Kh.Nm  <- function(datavec, Bndwdth) {
  KH <- numeric(length(datavec))
  xh<- datavec/Bndwdth
  KH <- dnorm(xh, 0, 1)/Bndwdth
  KH
}




#' Title  Nadaraya-Watson Kernel estimator at x0
#'
#' @param Xvec,Yvec  numeric vectors of data values, Xvec and Yvec must have the same length.
#' @param X0      a given value
#' @param Kernel  a character string indicating which kernel function is to be used. Use of "Ep", "Bw", or "Nm" for Epanechnikov, Biweight or Normal kernel function.
#' @param Bndwdth a bandwidth of the kernel
#' @param Wt      a weight vector or a constant. For longitudinal data, Wt=1/N corresponds to measurement uniform weight and Wt=1/(nni) corresponds subject uniform weight.
#'
#' @return   The kernel estimator at x0
#' @export
#' @examples
#' X <- seq(0, 1, len=100)
#' Y <- (X- 0.5)^3 - 2*(X-0.5)^2+ rnorm(100, 0, 0.1)
#' NW.WtKernel(X, Y,  X0=0.5, Kernel="Ep", Bndwdth=0.3, Wt=1 )
#' NW.WtKernel(X, Y,  X0=0.5, Kernel="Nm", Bndwdth=0.3, Wt=1 )
#' @references{ \enumerate{
#'  \item Fan, J. and Gijbels, I. Local Polynomial Modeling and Its Applications.
#'      Chapman & Hall, London, United Kingdom, 1996.
#'  \item  Wu, C.O. and Tian, X.   Nonparametric Models for Longitudinal Data. Chapman & Hall/CRC. To appear.}}
NW.WtKernel <- function(Xvec, Yvec, X0, Kernel="Ep", Bndwdth,  Wt=1)
{
  Xi0 <- Xvec- X0

  if (Kernel=="Ep")      { Wi <- Kh.Ep(Xi0,Bndwdth )*Wt }
  else if (Kernel=="Bw") { Wi <- Kh.Bw(Xi0,Bndwdth )*Wt }
  else if (Kernel=="Nm") { Wi <- Kh.Nm(Xi0,Bndwdth)*Wt  }

  Est0 <- sum(Wi*Yvec)/sum(Wi)
  Est0
}




#'  Nadaraya-Watson Kernel estimator
#' @param Xint  a vector of x interval to generate the local linear fit
#' @param Xvec,Yvec  numeric vectors of data values, Xvec and Yvec must have the same length.
#' @param bw   a bandwidth of the kernel
#' @param Wt   a weight vector
#' @param Kernel  a character string indicating which kernel function is to be used. Use of "Ep", "Bw", or "Nm" for Epanechnikov, Biweight or Normal kernel function.
#'
#' @export
#' @examples
#' X <- seq(0, 1, len=100)
#' Y <- (X- 0.5)^3 - 2*(X-0.5)^2+ rnorm(100, 0, 0.1)
#' kernel.fit(seq(0,1,0.1), X, Y, Kernel="Ep", bw=0.1, Wt=1   )
#' @references{ \enumerate{
#'  \item Fan, J. and Gijbels, I. Local Polynomial Modeling and Its Applications.
#'      Chapman & Hall, London, United Kingdom, 1996.
#'  \item  Wu, C.O. and Tian, X.   Nonparametric Models for Longitudinal Data. Chapman & Hall/CRC. To appear.}}
kernel.fit <- function(Xint , Xvec, Yvec, bw,  Kernel="Ep",  Wt=1)
{
  nID <- length(Xint)
  Yfit <- numeric(nID)
  for (i in 1:nID)
  {
    X1  <-  Xint[i]
    Yfit[i]<- NW.WtKernel(Xvec,  Yvec, X0=X1, Kernel, Bndwdth=bw,  Wt)
  }
  Yfit
}






