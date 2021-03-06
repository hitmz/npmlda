##' Least square local linear fit at t0
#'
#' @param t0       a given time point
#' @param Tvec,Yvec   numeric vectors of time and outcome values, Tvec and Yvec must have the same length.
#' @param X1,X2,X3   three covariate vectors
#' @param Weight    the weight vector
#' @param Bndwdth   a bandwidth of the Epanechnikov kernel
#' @importFrom stats lm
#' @export
#' @references{ Wu, C.O. and Tian, X.  Nonparametric Models for Longitudinal Data: With Implementation in R. Chapman & Hall/CRC.
#'             2018.}
#' @examples  # see usage of LocalLm.Beta
#'
#'
LocalLm.Beta.t0 <- function(t0, Tvec, X1, X2, X3, Yvec, Bndwdth, Weight )
{
  Tij0 <-  Tvec- t0
  Wij  <- Kh.Ep(Tij0 ,Bndwdth )
  Base <- data.frame(1, Tij0)
  DX1 <- Base *X1
  DX2 <- Base *X2
  DX3 <- Base *X3
  DesignX<- as.matrix(data.frame(Base , DX1, DX2, DX3))
  Est0 <- lm(Yvec~ DesignX-1, weights=Weight * Wij )$coef   # 8 values, only need intercept bi0
  Beta.t0 <-  Est0[c(1,3,5,7)]
  Beta.t0
}


#' Least square local linear fit
#'
#' @param Tint   a time interval
#' @param Tvec,Yvec  numeric vectors of time and outcome values, Tvec and Yvec must have the same length.
#' @param X1,X2,X3   three covariate vectors
#' @param Weight    the weight vector
#' @param Bndwdth   a bandwidth of the Epanechnikov kernel
#' @importFrom stats lm
#' @export
#' @references{ Wu, C.O. and Tian, X.  Nonparametric Models for Longitudinal Data: With Implementation in R. Chapman & Hall/CRC.
#'             2018.}
#' @examples
#' data(NGHS)
#' NGHS$Black <- (NGHS$RACE==2)*1
#' NGHS<- NGHS[!is.na(NGHS$SBP) & !is.na(NGHS$BMIPCT) & !is.na(NGHS$HTPCT ),]
#' Ct <-   data.frame(table(NGHS$ID))
#' names(Ct)<- c('ID', 'ni')
#' NGHS<- merge(NGHS, Ct, by= 'ID')
#' nID<- dim(Ct)[1]
#' Age.grid <- seq(9, 19, by=0.5) #21
#' NGHS$HTPCTc<- NGHS$HTPCT-50
#' NGHS$BMIPCTc<- NGHS$BMIPCT-50
#' Beta <- with(NGHS, LocalLm.Beta(Age.grid, AGE, X1=Black, X2=HTPCTc, X3=BMIPCTc, SBP, Bndwdth=3.5, Weight=1/ni))
LocalLm.Beta<- function(Tint, Tvec, X1, X2, X3, Yvec, Bndwdth, Weight )
{
  nt<- length(Tint)
  Beta.Vec<- matrix(NA, ncol=4, nrow=nt)
  for (i in 1:nt)
  {
    t1  <- Tint[i]
    Beta.Vec[i,] <-   LocalLm.Beta.t0( t1, Tvec, X1, X2, X3, Yvec, Bndwdth, Weight)
  }
  Beta.Vec
}



