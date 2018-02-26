#' Polynomial-spline fit with equally-spaced knots
#'
#' @param Xint    a vector of x interval to generate the local linear fit
#' @param Xvec,Yvec  numeric vectors of data values, Xvec and Yvec must have the same length.
#' @param nKnots    number of equally-spaced knots
#' @param Degree    degree of polynomial splines
#' @param Wt        a weight vector or a constant. For longitudinal data, Wt=1/N corresponds to measurement uniform weight and Wt=1/(nni) corresponds subject uniform weight.
#' @importFrom splines bs
#' @importFrom stats lm
#' @importFrom stats coef
#' @export
#' @references{ Wu, C.O. and Tian, X.   Nonparametric Models for Longitudinal Data. Chapman & Hall/CRC.
#'             To appear.}

Spline.fit <- function(Xint,  Xvec, Yvec , nKnots=2, Degree =3, Wt=1)
{
  DR<-  range(Xvec)
  Knots<- seq(DR[1], DR[2], length=nKnots+2)[c(-1, -(nKnots+2))]
  bs.Xsub <- bs(Xvec, knots =Knots, degree=Degree )
  if (length(Wt)==1 ) {Lm.fitW <- lm(Yvec ~ bs.Xsub) }
  else   { Lm.fitW <-lm(Yvec ~ bs.Xsub, weights=Wt)  }

  New.X<- bs(x= Xint, knots=Knots, degree=Degree )
  Spline.EstW <- as.numeric(cbind(1, New.X) %*% coef(Lm.fitW))
  Spline.EstW
}


#' Leave one-subject cross-validation score for spline fit
#'
#' @param Xvec,Yvec  numeric vectors of data values, Xvec and Yvec must have the same length.
#' @param ID    subject ID of the data value
#' @param nKnots    number of equally-spaced knots
#' @param Degree    degree of polynomial splines
#' @param Wt  a weight vector. For longitudinal data, Wt=1/N corresponds to measurement uniform weight and Wt=1/(nni) corresponds subject uniform weight.
#' @importFrom splines bs
#' @importFrom stats lm
#' @export
#' @references{ Wu, C.O. and Tian, X.   Nonparametric Models for Longitudinal Data. Chapman & Hall/CRC.
#'             To appear.}
CVspline <- function( Xvec, Yvec,  ID,  nKnots, Degree, Wt )
{
  NN <-length( Yvec)
  Yest <- numeric(NN)
  nID <- length(unique(ID))
  for (i in 1:nID)
  {
    #print(i)
    Xsub <- Xvec[ ID != i ]
    Ysub <- Yvec[ ID != i ]
    X.IDi <- Xvec[ ID == i ]
    Weighti<- Wt[ ID != i ]
    Yest[ID == i] <-  Spline.fit(X.IDi, Xsub, Ysub , nKnots, Degree , Weighti)
  }
  sum( Wt*(Yvec- Yest)^2)
}
