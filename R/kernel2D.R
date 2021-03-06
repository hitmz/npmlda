#' Multiplicative Epanechnikov Kernel (2-dim)
#'
#' @param datavec1, datavec2   two numeric vectors of same length
#' @param Bndwdth1, Bndwdth2   two bandwidths for two vectors
#'
#' @export
#' @return 2-dim kernel function result
#' @examples
#' Kh2D(2:7, 2:7, 5, 5)
Kh2D  <- function(datavec1,datavec2, Bndwdth1, Bndwdth2) {
  KH1 <-  numeric(length(datavec1))
  KH2 <-  numeric(length(datavec2))

  u1 <- datavec1/Bndwdth1
  u2 <- datavec2/Bndwdth2

  # For each row of (u1i, u2i) is the vector input for bandwith function

  Ind1<- ( abs(u1 ) >=1 )
  Ind2<- ( abs(u2 ) >=1 )

  KH1[Ind1] <- 0
  KH1[!Ind1] <- 0.75*(1-u1[!Ind1]^2)/Bndwdth1

  KH2[Ind2] <- 0
  KH2[!Ind2] <- 0.75*(1-u2[!Ind2]^2)/Bndwdth2

  KH<- KH1*KH2
  KH
}



#' 2-dim Kernel function for longitudinal data
#'
#' @param IDls        the vector of subject ID in a longitudinal sample
#' @param Xvec, Yvec  numeric vectors of data values, Xvec and Yvec must have the same length
#' @param X01, X02    two given values of Xvec
#' @param Bndwdth1,Bndwdth2  two given bandwidths
#'
#' @return 2-dim kernel fit result
#' @export
#' @references{ Wu, C.O. and Tian, X.  Nonparametric Models for Longitudinal Data: With Implementation in R. Chapman & Hall/CRC.
#'             2018.}
Kernel2D <- function(IDls, Xvec, Yvec, X01,X02,  Bndwdth1, Bndwdth2 )
{
  IDD1 <- as.numeric(names(table(IDls)))
  nID <- length(IDD1)
  IndS1 <- KerS1<- numeric(nID )

  for ( i in 1:nID )
  {
    #if (i%%200==0) print(i)
    Xi <- Xvec[IDls==IDD1[i]]
    Yi <- Yvec[IDls==IDD1[i]]

    Tvec1 <- Xi - X01
    Tvec2 <- Xi - X02
    Ind1 <-  (1: length(Xi))[abs(Tvec1)<=Bndwdth1 ]
    Ind2 <-  (1: length(Xi))[abs(Tvec2)<=Bndwdth2 ]

    nn1 <- length( Ind1)
    nn2 <- length( Ind2)

    Tvec1.ker <- rep( Tvec1[Ind1 ] , rep(nn2,nn1 ) )
    Tvec2.ker <- rep( Tvec2[Ind2 ] , nn1 )

    KerU.i <- Kh2D(Tvec1.ker, Tvec2.ker, Bndwdth1, Bndwdth2)

    Indic.s1 <-  Yi[Ind1]
    Indic.s2 <-  Yi[Ind2]

    IndS1.i <- as.vector(Indic.s2%*%t(Indic.s1) )

    IndS1[i] <- sum(   KerU.i*IndS1.i  )
    KerS1[i] <- sum(KerU.i)
  }
  Est <-  sum(IndS1)/sum(KerS1)
  Est
}


