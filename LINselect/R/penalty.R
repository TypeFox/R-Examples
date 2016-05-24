###########################################################################
#
# calcul de la penalite
#
###########################################################################

EDkhi1 <- function(x,D,N,q) {
  # ---------------------------------------------------------------
  # Case : L not too  small
  # called by calcEDkhi
  # INPUT
  # x scalar or vector
  # D, N : positive integer
  # q : scalar
  # OUTPUT
  # values in x : scalar or vector
  # CALL:
  # This function is called by 'calcEDkhi'
  # ---------------------------------------------------------------
  fct1 <- pf(q=x/(D+2), df1=D+2, df2=N, lower.tail = FALSE)
  fct2 <- (x/D)*pf(q=(N+2)*x/(N*D), df1=D, df2=N+2, lower.tail = FALSE)
  return(q+fct2-fct1)
} #  fin EDkhi1

EDkhi2  <- function(x,D,N,logq) {
  # ---------------------------------------------------------------
  # Case : L small
  # called by calcEDkhi
  # INPUT
  # x scalar or vector
  # D, N : positive integer
  # logq : scalar
  # OUTPUT
  # values in x : scalar or vector
  # CALL:
  # This function is called by 'calcEDkhi'
 # ---------------------------------------------------------------
  t1 <- log(2*(2*x+N*D)/(N*(N+2)*x))-lbeta(1+D/2,N/2)
  t3 <- (N/2)*log(N/(N+x))+(D/2)*log(x/(N+x))
  return(logq-t1-t3)
} #  fin EDkhi2

penalty <- function
  ### Calculate the penalty function for estimators selection.
  ##reference<<  See Baraud  et al. 2010 
  ## \url{http://hal.archives-ouvertes.fr/hal-00502156/fr/} \cr

(Delta, ##<<vector with \code{Dmax}+1 components : weights in the penalty function.
 n, ##<<integer : number of observatons.
 p, ##<<integer : number of variables.
 K ##<<scalar : constant in the penalty function.
 ) {
  # SH 19/02/2013 : penalite et penalty different legerement dans les
  # entrees. Il faudra corriger ce doublon
  # called by VSelect
   # VOIR: AB:08/03/2011 a rajouté l'argument p
  # SH 19/02/2013 ok
  # ---------------------------------------------------------------
  # INPUT
  # Delta vector of length Dmax+1
  # n,p: integer
  # K : scalar 
  # OUTPUT
  # vector of length dmax+1

  # CALL:
  # This function calls 'calcEDkhi'
  # ---------------------------------------------------------------
  Dmax <- length(Delta)-1
  Dm <- 0:Dmax
  pen <- Dm*0
  Nm <- n-Dm
  Lm <- Delta
  EDkhi <- calcEDkhi(p, Dm+1,Nm-1,Lm)
  if ( EDkhi$err != 0)
    {
      # error or warning cases
      mess <- paste(EDkhi$mess,
                    "\n err=", EDkhi$err,
                    "n n=", n,
                    "Dm=", Dm)
      if (EDkhi$err ==2) {
        ##note<< The values of the penalty function greater than
        ##1e+08 are set to 1e+08.
        warning(mess)
      }
      else
        stop(mess)
      ##note<<  If for some \code{Delta(d)} the
        ##equation \eqn{\phi}\code{(x) = exp(-Delta(d)/(d+1))}  has no
        ##solution, then the execution is stopped.
    } # fin err
    EDkhiR <- EDkhi$EDkhi2


  pen <- K*(Nm/(Nm-1))*EDkhiR
  names(pen) <- Dm
  
  return(pen)
  ### A vector with the same length as Delta: for each \code{d}=0, ..., \code{Dmax}, let
  ### \code{N}=\code{n}-\code{d}, \code{D}=\code{d+1} and \cr
  ### \code{pen(d) = x K N/(N-1)} where  \code{x} satisfies
  ###\cr
  ###\cr
  ### \eqn{\phi}\code{(x) = exp(-Delta(d))}, when \code{Delta(d)<50}, 
  ###    \cr
  ###     where  \eqn{\phi}\code{(x)=pf(q=x/(D+2),df1=D+2,df2=N-1,lower.tail=F)-(x/D)pf(q=(N+1)x/D(N-1),df1=D,df2=N+1,lower.tail=F)}
  ###\cr
  ###\cr
  ### \eqn{\psi}\code{(x) = Delta(d)}, when
  ###   \code{Delta(d)}\eqn{\ge}\code{50}, \cr
  ###    where  \eqn{\psi}\code{(x)=lbeta(1+D/2,(N-1)/2)-log(2(2x+(N-1)D)/((N-1)(N+2)x))-(N-1)/2log((N-1)/(N-1+x))-(D/2)log(x/(N-1+x))
 ### }
} #  fin penalty

calcEDkhi <- function(p, D,N,L)
  # VOIR: AB:08/03/2011 a rajouté l'argument p
# EDkhi <- calcEDkhi(Dm+1,Nm-1,LmC)
  # called by penalite and penalty
# D<-Dm+1; N=Nm-1; L=Lm
# SH : Nouvelle version
  # ---------------------------------------------------------------
  # INPUT
  # p : integer
  # D, N vectors of integers
  # K : vector
  # OUTPUT
  # array : list with components
  #   EDkhi2 value of EDkhi
  #   err : integer. Error or Warning code.
  #   mess : error message if err not equal 0
  # CALL:
  # This function calls 'EDkhi1', 'EDkhi2'.
  # It is called by 'penalty'
  # ---------------------------------------------------------------
   {
     Bsup <- 1e+08
     err <- 0
     mess <- NULL
     xq <- 0*D
     for (i in 2:length(D)) {
       if (xq[i-1] > Bsup) {
         xq[i:length(D)] <- Inf
         err <- 2
         mess <- paste("the values of the penalty function greater than", Bsup, "are set to Inf")
        break
      }
      if (i==2) {
        xInf <- 0
        if (N[2]>5) {
          Delta <- (L[2]+log(5)+1/N[2])/(1-5/N[2])
          U <- sqrt((1+2*D[2]/(N[2]+2))*2*Delta/D[2])
          xSup <- D[2]*(1+exp(2*Delta/(N[2]+2))*U)**2
        }
        if (N[2]<=5) {
          xSup <- seq(10,100,by=10)*p 
          fxSup <- EDkhi1(xSup,D[i],N[i],exp(-L[i]))
          if (max(fxSup) <0 ) {
            err <- 4
            mess <- paste("the values of the penalty function cannot be calculated\n", xq)
          break
          }
          if (max(fxSup) >=0 ) {
            xSup <- min(xSup[fxSup>=0])
          }
        }
      }
      else {
#       (i>2)
        xInf <- xq[i-1]
        xSup <- xq[i-1]*seq(10,100,by=10)
        fxSup <- EDkhi1(xSup,D[i],N[i],exp(-L[i]))
        if (max(fxSup) <0 ) {
          err <- 3
          mess <- paste("the values of the penalty function cannot be calculated\n", xq)
          break
        }
        xSup <- min(xSup[fxSup>=0])
      }
      if (L[i] < 50) {
        xx <-
          try(uniroot(EDkhi1,lower=xInf,upper=xSup,D=D[i],
                      N=N[i],q=exp(-L[i])))
        if (!is.list(xx)) {
          err <- 1
          mess <- xx
          break
        }
        xq[i] <- xx$root
      }
      else {
# L[i] >= 50
        fxSup <- EDkhi2(xSup,D[i],N[i],-L[i])
        if (fxSup <0) xSup<- 2*xSup
        z <-
          try(uniroot(EDkhi2,lower=xInf,upper=xSup,D=D[i],
                      N=N[i],logq=-L[i]))
        if (!is.list(xx)) {
          err <- 1
          mess <- xx
          break
        }
        xq[i] <- z$root
      }
    }
    return(list(EDkhi2=xq,err=err, mess=mess))
  } #  fin calcEDkhi

