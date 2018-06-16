#' Function Xi(s)
#'
#' @param s  a number or a vector
#'
#' @export
#' @return   value of the function with give s
#' @examples
#' Xi(0)
#' Xi(c(-1000, -10,-5, 0, 5,10, 1000 ))
Xi<- function(s)
{
  s<- as.vector(s)
  Xi <- numeric(length(s))
  Xi[s == 0]  <-  0.5
  Xi[s > 200] <-  0
  Xi[s < -200] <-  1
  Ind<- ( s !=0 & abs(s)<=200)
  si <- s[Ind]
  Xi[Ind] <-  ((si -1)*exp(si)+1)/(exp(si)-1)^2
  Xi
}



#' Derivative of the function Xi(s)
#'
#' @param s   a number or a vector
#'
#' @return  value of the function DXi with give s
#' @export
#'
#' @examples
#' DXi(c(-1000, -10,-5, 0, 5,10, 1000 ))
DXi<- function(s)
{
  s<- as.vector(s)
  DXi <- numeric(length(s))
  DXi[s == 0]  <-  1/6
  DXi[abs(s) > 200 ] <-   0
  Ind<- ( s != 0 & abs(s)<=200)
  si <- s[Ind]
  DXi[Ind] <-  (exp(si)*(si-2)+si+2)*exp(si)/(exp(si)-1)^3
  DXi
}




#' An equation solver with Newton's method with 2 variables
#'
#' @param Zij  2-dim covariate vector
#' @param b0,Ub   inital values
#' @param Indicator  Indicator of Yi1> Yi2
#' @param difflmt   limit to stop the interations
#' @param MaxIter  maximum no. of interations
#'
#' @return the root of the equation
#' @export
Newton2var<- function(Zij, b0 , Ub, Indicator, difflmt=1e-14, MaxIter=100)
{
  i<- 0

  #Initial the 3*3 Jacobian matrix.
  DU<- matrix(0, 2,2)

  while( (t(Ub)%*%Ub)  >  difflmt)
  {
    diag(DU) <- c( sum(Zij[,1]^2* DXi(Zij%*%b0)), sum(Zij[,2]^2* DXi(Zij%*%b0)) )

    DXi.Z    <- DXi(Zij%*%b0) ; #nT*nT*1
    DU[1,2]  <-  sum(Zij[,1]*Zij[,2]* DXi.Z )
    DU[2,1]  <-  sum(Zij[,1]*Zij[,2]* DXi.Z )

    Delta <- solve(DU, Ub)
    # newton update
    b0 <-   b0 - Delta

    i<- i+1

    if(i > MaxIter)
    {
      print('Stopped after maximum iteration')
      return (b0 )
      break
    }

    Diff<-  as.vector(Indicator)*1-  Xi(Zij%*%b0); #nT*nT*1

    U1  <-  sum(Zij[,1]* Diff)
    U2  <-  sum(Zij[,2]* Diff)

    Ub <- c(U1, U2)
  }
  #------------------------------------------------------------
  #output results
  b0
}



#' An equation solver with Newton's method with 1 variable
#'
#' @param Z12vec
#' @param h0
#' @param Vh
#' @param HZB
#' @param Ind.Y
#' @param Diff
#' @param ORR
#' @param MaxIter
#'
#' @return the root of the equation
#' @export
Newton1var<- function(Z12vec, h0 , Vh ,HZB, Ind.Y , Diff=1e-8,  ORR, MaxIter=100)
{
  hint <- h0
  j<- 0

  while( abs(Vh)  >  Diff)
  {

    DV <- sum(exp(HZB)/((1+ exp(HZB))^2)) #  fixed the square, outside of (1+e^x)^2
    #correction!!!: note for square 7/16/2015

    #print(c("j",j))
    #print(c("h0", h0))
    #print(c("Vh", Vh))
    #print(c("DV", DV))

    if (abs(DV)< 1e-6)
    {
      return (h0 )
      break
    }

    Delta.h <- Vh/DV

    #print(c("Delta", Delta.h))

    # newton update
    h0 <-   h0 - Delta.h

    #print(c("h0", h0))

    j<- j+1

    if(j > MaxIter)
    {
      print('Stopped after maximum iteration')
      return (h0 )
      break
    }

    HZB <- h0 + Z12vec%*%ORR
    Vh <- Ind.Y - sum(1/(1+exp(HZB)))

    if (abs(Vh-Ind.Y)< 1e-8)
    {# print('Stopped not converge')
      if (abs(h0)<30) {return (h0)}  else  {return (hint)}
      break
    }

  }

  #------------------------------------------------------------
  #output results
  h0
}














