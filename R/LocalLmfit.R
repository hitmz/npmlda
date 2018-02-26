## Local linear fit at X0 ##

#' Local linear fit at X0 with Epanechnikov kernel
#'
#' @param Xvec,Yvec  numeric vectors of data values, Xvec and Yvec must have the same length.
#' @param X0        a given value
#' @param Bndwdth   a bandwidth of the kernel
#' @param Wt        a weight vector or a constant. For longitudinal data, Wt=1/N corresponds to measurement uniform weight and Wt=1/(nni) corresponds subject uniform weight.
#' @export
#'
#' @examples  # see usage of LocalLm
LocalLm.X0 <- function(Xvec, Yvec, X0, Bndwdth , Wt=1)
{
  Xi0 <- Xvec- X0
  Wi  <- Kh.Ep(Xi0,Bndwdth ) * Wt
  Sn2 <- sum(Wi*Xi0^2)
  Sn1 <- sum(Wi*Xi0 )
  Est0 <- (Sn2*sum(Wi*Yvec)-Sn1*sum(Wi*Xi0*Yvec))/(Sn2*sum(Wi)-Sn1^2)
  Est0
}


## Local linear fit for a x-interval##

#' Local linear fit with Epanechnikov kernel
#'
#' @param Xint  a vector of x interval to generate the local linear fit
#' @param Xvec,Yvec  numeric vectors of data values, Xvec and Yvec must have the same length.
#' @param bw   a bandwidth of the kernel
#' @param Wt   a weight vector
#'
#' @export
#'
#' @examples
#' data(BMACS)
#' Time.int<- seq(0.1,5.9,  by=0.1)
#' LocalFit.Y <- with(BMACS, LocalLm(Time.int, Time, CD4, bw=0.9, Wt=1))
LocalLm <- function(Xint, Xvec, Yvec, bw , Wt=1)
{
  nID <- length(Xint)
  Yfit <- numeric(nID)
  for (i in 1:nID)
  {
    X1  <-  Xint[i]
    Yfit[i] <- LocalLm.X0(Xvec,  Yvec, X1, Bndwdth=bw , Wt)
  }
  Yfit
}

## Local linear CV function ##
#' Leave one-subject cross-validation score for local linear fit
#'
#' @param Xvec,Yvec  numeric vectors of data values, Xvec and Yvec must have the same length.
#' @param bw  a bandwidth of the Epanechnikov kernel
#' @param ID  subject ID of the data value
#' @param Wt  a weight vector, may be subject-specific. a weight vector or a constant. For longitudinal data, Wt=1/N corresponds to measurement uniform weight and Wt=1/(nni) corresponds subject uniform weight.
#'
#' @export
#'
CVlm <- function( Xvec, Yvec, bw, ID, Wt)
{
  NN <-length( Yvec)
  Yest <- numeric(NN)
  nID <- length(unique(ID))
  for (i in 1:nID)
  {
    Xsub <- Xvec[ ID != i ]
    Ysub <- Yvec[ ID != i ]
    X.IDi <- Xvec[ ID == i ]
    Weighti<- Wt[ ID != i ]
    Yest[ID == i] <-  LocalLm(Xint=X.IDi, Xsub, Ysub , bw , Weighti)
  }
  sum( Wt*(Yvec- Yest)^2)
}


