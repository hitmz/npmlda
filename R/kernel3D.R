#' Multiplicative Epanechnikov Kernel (3-dim)
#'
#' @param datavec1, datavec2, datavec3  three numeric vectors of same length
#' @param Bndwdth1, Bndwdth2, Bndwdth3  three bandwidths for three vectors
#'
#' @return 3-dim kernel function result
#' @export
#'
Kh3D  <- function(datavec1,datavec2, datavec3,Bndwdth1, Bndwdth2, Bndwdth3) {

  KH1 <-  numeric(length(datavec1))
  KH2 <-  numeric(length(datavec2))
  KH3 <-  numeric(length(datavec3))

  u1 <- datavec1/Bndwdth1
  u2 <- datavec2/Bndwdth2
  u3 <- datavec3/Bndwdth3

  # For each row of (u1i, u2i, u3i) is the vector input for bandwith function

  Ind1<- ( abs(u1 ) >=1 )
  Ind2<- ( abs(u2 ) >=1 )
  Ind3<- ( abs(u3 ) >=1 )

  KH1[Ind1] <- 0
  KH1[!Ind1] <- 0.75*(1-u1[!Ind1]^2)/Bndwdth1

  KH2[Ind2] <- 0
  KH2[!Ind2] <- 0.75*(1-u2[!Ind2]^2)/Bndwdth2

  KH3[Ind3] <- 0
  KH3[!Ind3] <- 0.75*(1-u3[!Ind3]^2)/Bndwdth3


  KH<- KH1*KH2*KH3
  KH
}



#' 3-dim Kernel function for longitudinal data to get Pr(y1(t1),y2(t2)|x(t1))
#'
#' @param IDls       the vector of subject ID in a longitudinal sample
#' @param Y,X,Time   numeric vectors of outcome, covariate and time of the the same length
#' @param T1,T2      twp given time points
#' @param X0          a given covariate value
#' @param Bndwdth1,Bndwdth2,Bndwdth3    three bandwidths around two time and one covariate value
#'
#' @return  3-dim Kernel function results
#' @export
#'
Kernel3D <- function(IDls=ID, Y, Time, X,  T1, T2, X0, Bndwdth1, Bndwdth2, Bndwdth3 )
{

  IDD1 <- as.numeric(names(table(IDls)))
  nID <- length(IDD1)
  IndS1 <- KerS1<- numeric(nID )

  for ( i in 1:nID )
  {
    #if (i%%200==0) print(i)
    Xi <- Time[IDls==IDD1[i]]
    Yi <- Y[IDls==IDD1[i]]
    Hti <- X[IDls==IDD1[i]]

    Tvec1 <- Xi - T1 # two time points S1 and S2
    Tvec2 <- Xi - T2

    Hvec1 <- Hti  - X0

    Ind1 <-  (1: length(Xi))[abs(Tvec1)<=Bndwdth1  & abs(Hvec1)<= Bndwdth3 ]
    Ind2 <-  (1: length(Xi))[abs(Tvec2)<=Bndwdth2 ]

    nn1 <- length( Ind1) # data points around S1
    nn2 <- length( Ind2) # data points around S2

    Hvec.ker  <- rep(Hvec1[Ind1 ],rep(nn1*nn2, nn1))      # length=nn1*(nn1*nn2),repeat each h (nn1*nn2) time
    Tvec1.ker <- rep(rep(Tvec1[Ind1 ], rep(nn2,nn1)),nn1)
    Tvec2.ker <- rep(rep(Tvec2[Ind2 ], nn1) ,nn1)

    KerU.i <- Kh3D(Tvec1.ker, Tvec2.ker, Hvec.ker,Bndwdth1, Bndwdth2, Bndwdth3)

    Indic.s1 <-  Yi[Ind1]
    Indic.s2 <-  Yi[Ind2]

    IndS1.i <- rep(as.vector(Indic.s2%*%t(Indic.s1) ), nn1)


    IndS1[i] <- sum(   KerU.i*IndS1.i  )
    KerS1[i] <- sum(KerU.i)
  }
  Est <-  sum(IndS1)/sum(KerS1)

  Est

}


#' 3-dim Kernel function for longitudinal data to get  Pr(y2(t2)|x(t1))
#'
#' @param IDls   the vector of subject ID in a longitudinal sample
#' @param Y,X,Time   numeric vectors of outcome, covariate and time of the the same length
#' @param T1,T2      twp given time points
#' @param X0          a given covariate value
#' @param Bndwdth1,Bndwdth2,Bndwdth3    three bandwidths around two time and one covariate value
#'
#' @return 3-dim Kernel function results
#' @export
#'
Kernel3D.S2 <- function(IDls=ID,  Y, Time, X,  T1, T2, X0, Bndwdth1, Bndwdth2, Bndwdth3 )
{
  IDD1 <- as.numeric(names(table(IDls)))
  nID <- length(IDD1)
  IndS1 <- KerS1<- numeric(nID )

  for ( i in 1:nID )
  {
    #if (i%%200==0) print(i)
    Xi <- Time[IDls==IDD1[i]]
    Yi <- Y[IDls==IDD1[i]]
    Hti <- X[IDls==IDD1[i]]

    Tvec1 <- Xi - T1
    Tvec2 <- Xi - T2

    Hvec1 <- Hti  -X0

    Ind1 <-  (1: length(Xi))[abs(Tvec1)<=Bndwdth1  & abs(Hvec1)<= Bndwdth3 ]
    Ind2 <-  (1: length(Xi))[abs(Tvec2)<=Bndwdth2 ]

    # need to get u1, u2
    nn1 <- length( Ind1)
    nn2 <- length( Ind2)


    Hvec.ker<-   rep(Hvec1[Ind1 ],rep(nn1*nn2, nn1))
    Tvec1.ker <- rep(rep(Tvec1[Ind1 ], rep(nn2,nn1)),nn1)
    Tvec2.ker <- rep(rep(Tvec2[Ind2 ], nn1) ,nn1)

    KerU.i <- Kh3D(Tvec1.ker, Tvec2.ker, Hvec.ker,Bndwdth1, Bndwdth2, Bndwdth3)


    IndS1.i <-  rep(rep(Yi[Ind2 ], nn1) ,nn1)

    IndS1[i] <- sum(   KerU.i*IndS1.i  )
    KerS1[i] <- sum(KerU.i)
  }
  Est <-  sum(IndS1)/sum(KerS1)

  Est

}







