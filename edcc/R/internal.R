  .echXbar1 <- function(x,delta=2,lambda=.05,P0=NULL,P1=NULL,C0=NULL,C1=NULL,Cr=25,Cf=50,T0=0.0167,Tc=1,Tf=0,Tr=0,a=1,b=.1,d1=1,d2=1, sided = "two"){
    h <- x[1]; L <- x[2]; n <- x[3]
    delta.std <- delta*sqrt(n)
    if(sided == "one"){
      ARL1 <- 1/pnorm(-L)
      ARL2 <- 1/pnorm(-L + abs(delta.std))
    }
    if(sided == "two"){
      alpha <- 2*pnorm(-L)
      beta <- pnorm(L - delta.std)-pnorm(- L- delta.std)
      ARL1 <- 1/alpha
      ARL2 <- 1/(1-beta)
    }
    tau <- (1-(1+lambda*h)*exp(-lambda*h))/(lambda*(1-exp(-lambda*h)))
    s <- 1/(exp(lambda*h)-1)
    if(!is.null(P0)&!is.null(P1)){
      ECT <- 1/lambda+(1-d1)*s*Tf/ARL1 - tau + n*T0 + h*ARL2 + Tc + Tr
      ECP <- P0/lambda + P1*(-tau+n*T0+h*ARL2+d1*Tc+d2*Tr) - s*Cf/ARL1 - Cr - (a+b*n)*(1/lambda-tau+n*T0+h*ARL2+d1*Tc+d2*Tr)/h
      ECH <- P0 - ECP/ECT
    }else
    if(!is.null(C0)&!is.null(C1)){
      ECT <- 1/lambda+(1-d1)*s*Tf/ARL1 - tau + n*T0 + h*ARL2 + Tc + Tr
      ECC <- C0/lambda + C1*(-tau+n*T0+h*ARL2+d1*Tc+d2*Tr) + s*Cf/ARL1+Cr+(a+b*n)*(1/lambda-tau+n*T0+h*ARL2+d1*Tc+d2*Tr)/h
      ECH <- ECC/ECT
    }else
    stop("You should at least give a pair of value to P0,P1 or C0,C1")
    return(ECH)
  }
  .echXbar2 <- function(x,n=5,delta=2,lambda=.05,P0=NULL,P1=NULL,C0=NULL,C1=NULL,Cr=25,Cf=50,T0=0.0167,Tc=1,Tf=0,Tr=0,a=1,b=.1,d1=1,d2=1, sided = "two"){
    h <- x[1]; L <- x[2]
    delta.std <- delta*sqrt(n)
    if(sided == "one"){
      ARL1 <- 1/pnorm(-L)
      ARL2 <- 1/pnorm(-L + abs(delta.std))
    }
    if(sided == "two"){
      alpha <- 2*pnorm(-L)
      beta <- pnorm(L - delta.std)-pnorm(- L- delta.std)
      ARL1 <- 1/alpha
      ARL2 <- 1/(1-beta)
    }
    tau <- (1-(1+lambda*h)*exp(-lambda*h))/(lambda*(1-exp(-lambda*h)))
    s <- 1/(exp(lambda*h)-1)
    if(!is.null(P0)&!is.null(P1)){
      ECT <- 1/lambda+(1-d1)*s*Tf/ARL1 - tau + n*T0 + h*ARL2 + Tc + Tr
      ECP <- P0/lambda + P1*(-tau+n*T0+h*ARL2+d1*Tc+d2*Tr) - s*Cf/ARL1 - Cr - (a+b*n)*(1/lambda-tau+n*T0+h*ARL2+d1*Tc+d2*Tr)/h
      ECH <- P0 - ECP/ECT
    }else
    if(!is.null(C0)&!is.null(C1)){
      ECT <- 1/lambda+(1-d1)*s*Tf/ARL1 - tau + n*T0 + h*ARL2 + Tc + Tr
      ECC <- C0/lambda + C1*(-tau+n*T0+h*ARL2+d1*Tc+d2*Tr) + s*Cf/ARL1+Cr+(a+b*n)*(1/lambda-tau+n*T0+h*ARL2+d1*Tc+d2*Tr)/h
      ECH <- ECC/ECT
    }else
    stop("You should at least give a pair of value to P0,P1 or C0,C1")
    return(ECH)
  }
.echCusum1 <- function(x,delta = 2,lambda = .01, P0 = NULL, P1 = NULL,C0 = NULL,C1 = NULL, Cr = 20, Cf = 10,T0 = 0, Tc = .1,Tf = .1,Tr = 0.2, a = .5, b = .1, d1 = 1, d2 = 1,sided = "one"){
  h <- x[1]; H <- x[2]; n <- x[3]
    delta.std <- sqrt(n)*delta          #standardization for delta
    k <- delta.std/2
    ARL1 <- as.numeric(xcusum.arl(k,H,0,sided=sided))
    ARL2 <- as.numeric(xcusum.arl(k,H,delta.std,sided=sided))
    tau <- (1-(1+lambda*h)*exp(-lambda*h))/(lambda*(1-exp(-lambda*h)))
    s <- 1/(exp(lambda*h)-1)
    if(!is.null(P0)&!is.null(P1)){
      ECT <- 1/lambda+(1-d1)*s*Tf/ARL1 - tau + n*T0 + h*ARL2 + Tc + Tr
      ECP <- P0/lambda + P1*(-tau+n*T0+h*ARL2+d1*Tc+d2*Tr) - s*Cf/ARL1 - Cr - (a+b*n)*(1/lambda-tau+n*T0+h*ARL2+d1*Tc+d2*Tr)/h
      ECH <- P0 - ECP/ECT
    }else
    if(!is.null(C0)&!is.null(C1)){
      ECT <- 1/lambda+(1-d1)*s*Tf/ARL1 - tau + n*T0 + h*ARL2 + Tc + Tr
      ECC <- C0/lambda + C1*(-tau+n*T0+h*ARL2+d1*Tc+d2*Tr) + s*Cf/ARL1+Cr+(a+b*n)*(1/lambda-tau+n*T0+h*ARL2+d1*Tc+d2*Tr)/h
      ECH <- ECC/ECT
    }else
    stop("You should at least give a pair of value to P0,P1 or C0,C1")
    return(ECH)
  }

.echCusum2 <- function(x,n=5,delta = 2,lambda = .01, P0 = NULL, P1 = NULL,C0 = NULL,C1 = NULL, Cr = 20, Cf = 10,T0 = 0, Tc = .1,Tf = .1,Tr = 0.2, a = .5, b = .1, d1 = 1, d2 = 1,sided = "one"){
  h <- x[1]; H <- x[2]
    delta.std <- sqrt(n)*delta          #standardization for delta
    k <- delta.std/2
    ARL1 <- as.numeric(xcusum.arl(k,H,0,sided=sided))
    ARL2 <- as.numeric(xcusum.arl(k,H,delta.std,sided=sided))
    tau <- (1-(1+lambda*h)*exp(-lambda*h))/(lambda*(1-exp(-lambda*h)))
    s <- 1/(exp(lambda*h)-1)
    if(!is.null(P0)&!is.null(P1)){
      ECT <- 1/lambda+(1-d1)*s*Tf/ARL1 - tau + n*T0 + h*ARL2 + Tc + Tr
      ECP <- P0/lambda + P1*(-tau+n*T0+h*ARL2+d1*Tc+d2*Tr) - s*Cf/ARL1 - Cr - (a+b*n)*(1/lambda-tau+n*T0+h*ARL2+d1*Tc+d2*Tr)/h
      ECH <- P0 - ECP/ECT
    }else
    if(!is.null(C0)&!is.null(C1)){
      ECT <- 1/lambda+(1-d1)*s*Tf/ARL1 - tau + n*T0 + h*ARL2 + Tc + Tr
      ECC <- C0/lambda + C1*(-tau+n*T0+h*ARL2+d1*Tc+d2*Tr) + s*Cf/ARL1+Cr+(a+b*n)*(1/lambda-tau+n*T0+h*ARL2+d1*Tc+d2*Tr)/h
      ECH <- ECC/ECT
    }else
    stop("You should at least give a pair of value to P0,P1 or C0,C1")
    return(ECH)
  }

.echEwma1 <- function(x,w=0.5,delta = 2,lambda = .05, P0 = NULL, P1 = NULL,C0 = NULL,C1 = NULL, Cr = 25, Cf = 10,T0 = 0.0167,Tc = 1, Tf = 0, Tr = 0, a = 1, b = .1,d1=1,d2=1,sided="two"){
   h <- x[1]; k <- x[2]; n <- x[3]
     delta.std <- sqrt(n)*delta #standardization fordelta
     ARL1 <- as.numeric(xewma.arl(w,k,0,sided=sided))
     ARL2 <- as.numeric(xewma.arl(w,k,delta.std,sided=sided))
     tau <- (1-(1+lambda*h)*exp(-lambda*h))/(lambda*(1-exp(-lambda*h)))
     s <- 1/(exp(lambda*h)-1)
     if(!is.null(P0)&!is.null(P1)){
       ECT <- 1/lambda+(1-d1)*s*Tf/ARL1 - tau + n*T0 + h*ARL2 + Tc + Tr
       ECP <- P0/lambda + P1*(-tau+n*T0+h*ARL2+d1*Tc+d2*Tr) - s*Cf/ARL1 - Cr - (a+b*n)*(1/lambda-tau+n*T0+h*ARL2+d1*Tc+d2*Tr)/h
       ECH <- P0 - ECP/ECT
     }else
     if(!is.null(C0)&!is.null(C1)){
       ECT <- 1/lambda+(1-d1)*s*Tf/ARL1 - tau + n*T0 + h*ARL2 + Tc + Tr
       ECC <- C0/lambda + C1*(-tau+n*T0+h*ARL2+d1*Tc+d2*Tr) + s*Cf/ARL1+Cr+(a+b*n)*(1/lambda-tau+n*T0+h*ARL2+d1*Tc+d2*Tr)/h
       ECH <- ECC/ECT
     }else
     stop("You should at least give a pair of value to P0,P1 or C0,C1")
     return(ECH)
   }
 .echEwma2 <- function(x,w=0.5,n=5,delta = 2,lambda = .05, P0 = NULL, P1 = NULL,C0 = NULL,C1 = NULL, Cr = 25, Cf = 10,T0 = 0.0167,Tc = 1, Tf = 0, Tr = 0, a = 1, b = .1,d1=1,d2=1,sided="two"){
   h <- x[1];  k <- x[2]
     delta.std <- sqrt(n)*delta #standardization fordelta
     ARL1 <- as.numeric(xewma.arl(w,k,0,sided=sided))
     ARL2 <- as.numeric(xewma.arl(w,k,delta.std,sided=sided))
     tau <- (1-(1+lambda*h)*exp(-lambda*h))/(lambda*(1-exp(-lambda*h)))
     s <- 1/(exp(lambda*h)-1)
     if(!is.null(P0)&!is.null(P1)){
       ECT <- 1/lambda+(1-d1)*s*Tf/ARL1 - tau + n*T0 + h*ARL2 + Tc + Tr
       ECP <- P0/lambda + P1*(-tau+n*T0+h*ARL2+d1*Tc+d2*Tr) - s*Cf/ARL1 - Cr - (a+b*n)*(1/lambda-tau+n*T0+h*ARL2+d1*Tc+d2*Tr)/h
       ECH <- P0 - ECP/ECT
     }else
     if(!is.null(C0)&!is.null(C1)){
       ECT <- 1/lambda+(1-d1)*s*Tf/ARL1 - tau + n*T0 + h*ARL2 + Tc + Tr
       ECC <- C0/lambda + C1*(-tau+n*T0+h*ARL2+d1*Tc+d2*Tr) + s*Cf/ARL1+Cr+(a+b*n)*(1/lambda-tau+n*T0+h*ARL2+d1*Tc+d2*Tr)/h
       ECH <- ECC/ECT
     }else
     stop("You should at least give a pair of value to P0,P1 or C0,C1")
     return(ECH)
   }
