##' Calculate the optimum parameters, n(sample size), h(sampling
##' interval) and L(number of s.d. from control limits to center line)
##' for Economic Design of the X-bar control chart .
##'
##' When parameter \code{par} is specified, optimization algorithms
##' would be used as default. \code{par} can be specified as:
##' \code{par = c(h, L)} where \code{h} and \code{L} are the initial
##' values of smapling interval and control limit when \code{n} is
##' specified; or \code{par = c(h, L, n)}. Good inital values may lead
##' to good optimum results.
##' 
##' When parameters \code{h}, \code{L}, \code{n} are all undefined,
##' \code{ecoXbar} function will try to find the global optimum point
##' to minimize the ECH (Expected Cost per Hour) using optimization
##' algorithms (\code{optim} function), but in this case \code{n}
##' would not be integer. It is usually helpful for the experimenter
##' to find the region where the optimum point may exist quickly. When
##' \code{h} and \code{L} are undefined but n is given as an integer
##' vector, \code{ecoXbar} function will try to find the optimum point
##' for each \code{n} value using optimization algorithms. When
##' \code{h}, \code{L} and \code{n} are all given, \code{ecoXbar}
##' function will use the "grid method" to calculate the optimum point,
##' that is ECH for all the combinations of the parameters will be
##' calculated. The "grid method" way is much slower than using
##' optimization algorithms, but it would be a good choice when
##' optimization algorithms fail to work well.
##'
##' For cost parameters either {P0, P1} or {C0, C1} is needed.  If P0
##' and P1 are given, they will be used first, else C0 and C1 will be
##' used.  For economic design of the X-bar chart, when \code{d1} and
##' \code{d2} are both 1, only if the difference between P0 and P1
##' keeps the same, the results are identical. If the difference
##' between C0 and C1 keeps the same, the optimum parameters are
##' almost the same but the ECH(Expected Cost per Hour) values will
##' change.
##'
##' \code{echXbar} is used to calculate the ECH (Expected Cost per
##' Hour) for one given design point.
##' 
##' @title Economic design for the X-bar control chart
##' @param h sampling interval. It can be a numeric vector or left
##' undefined. See 'Details'
##' @param L number of standard deviations from control limits to
##' center line. It can be a numeric vector or left undefined. See
##' 'Details'
##' @param n sample size. It can be an integer vector or left
##' undefined. See 'Details'
##' @param lambda we assume the in-control time follows an exponential
##' distribution with mean 1/lambda. Default value is 0.05.
##' @param delta shift in process mean in standard deviation units
##' when assignable cause occurs (delta = (mu1 - mu0)/sigma), where
##' sigma is the standard deviation of observations; mu0 is the
##' in-control process mean; mu1 is the out-of-control process mean.
##' Default value is 2.
##' @param P0 profit per hour earned by the process operating in
##' control. See 'Details'.
##' @param P1 profit per hour earned by the process operating out of
##' control(P0 > P1).
##' @param C0 cost per hour due to nonconformities produced while the
##' process is in control.
##' @param C1 cost per hour due to nonconformities produced while the
##' process is out of control.(C1 > C0)
##' @param Cr cost for searching and repairing the assignable cause,
##' including any downtime.
##' @param Cf cost per false alarm, including the cost of searching
##' for the cause and the cost of downtime if production ceases during
##' search.
##' @param T0 time to sample and chart one item.
##' @param Tc expected time to discover the assignable cause.
##' @param Tf expected search time when false alarm occurs.
##' @param Tr expected time to repair the process.
##' @param a fixed cost per sample.
##' @param b cost per unit sampled.
##' @param d1 flag for whether production continues during searches
##' (1-yes, 0-no). Default value is 1.
##' @param d2 flag for whether production continues during repairs
##' (1-yes, 0-no). Default value is 1.
##' @param nlevels number of contour levels desired. Default value is
##' 30. It works only when \code{contour.plot} is TRUE.
##' @param sided distinguish between one- and two-sided X-bar chart by
##' choosing ``one'' or ``two'' respectively.  When \code{sided =
##' ``one''}, \code{delta > 0} means the control chart for detecting a
##' positive shift, and vice versa. Default is ``two''.
##' @param par initial values for the parameters to be optimized
##' over. It can be a vector of length 2 or 3. See 'Details'
##' @param contour.plot a logical value indicating whether a contour
##' plot should be drawn. Default is FALSE. Only works when the
##' parameters h, L and n are all specified.
##' @param call.print a logical value indicating whether the "call"
##' should be printed on the contour plot. Default is TRUE
##' @param ... other arguments to be passed to \code{optim} function.
##' @return The \code{ecoXbar} function returns an object of class
##' "edcc", which is a list of elements \code{optimum},
##' \code{cost.frame}, \code{FAR} and \code{ATS}. \code{optimum} is a
##' vector with the optimum parameters and the corresponding ECH
##' value; \code{cost.frame} is a dataframe with the optimum
##' parameters and the corresponding ECH values for all given
##' \code{n}(if \code{n} is not specified, \code{cost.frame} won't be
##' returned); \code{FAR} indicates the false alarm rate during the
##' in-control time, which is calculated as lambda*(average number of
##' false alarm); \code{ATS} indicates the average time to signal
##' after the occurrence of an assignable cause, calculated as h*ARL2
##' - tau, where tau is the expected time of occurrence of the
##' assignable cause given it occurs between the i-th and (i+1)st
##' samples.
##' The \code{echXbar} function returns the calculated ECH value only.
##' @seealso \code{\link{ecoCusum}}, \code{\link{ecoEwma}},
##' \code{\link{contour}}, \code{\link[stats]{optim}},
##' \code{\link{update.edcc}}
##' @references Weicheng Zhu, Changsoon Park (2013), {edcc: An R
##' Package for the Economic Design of the Control Chart}. \emph{Journal
##' of Statistical Software}, 52(9),
##' 1-24. \url{http://www.jstatsoft.org/v52/i09/}
##' 
##' Douglas (2009). \emph{Statistical quality control: a modern
##' introduction}, sixth edition,  463-471.
##' 
##' Lorenzen and Vance (1986). The economic design of control charts:
##' a unified approach, \emph{Technometrics}, 28. 3-10.
##' @examples 
##' # Douglas (2009). Statistical quality control: a modern introduction, sixth edition, p470.
##' # In control profit per hour is 110, out of control profit per hour is 10
##' ecoXbar(P0=110,P1=10)
##' # In control profit per hour is 150, out of control profit per hour
##' # is 50, the result is identical with the previous one, because the
##' #difference between P0 and P1 are the same
##' ecoXbar(P0=150,P1=50)
##' # In control cost per hour is 0, out of control cost per hour is 100.
##' # The result is the same with the previous one
##' ecoXbar(C0=0,C1=100)
##' # The optimum parameters are the same with the previous one,
##' # but Cost values are different. See 'details'
##' ecoXbar(C0=10,C1=110)
##' # Based on the global optimum above, we specify the range of the
##' # parameters like this
##' x <- ecoXbar(h=seq(0.7,0.9,by=.01),L=seq(2.8,3.3,by=.01),n=4:6,P0=110,
##' P1=10,nlevels=50,contour.plot=TRUE)
##' x
##' # Modify the contour plot 
##' contour(x,nlevels=20,lty=2,col=2,call.print=FALSE)
##' # update the parameters
##' update(x,P0=NULL,P1=NULL,C0=10,C1=110)
##' @export
ecoXbar <- function(h, L, n, delta = 2, lambda = .05, P0 = NULL, P1 = NULL, C0 = NULL, C1 = NULL, Cr = 25, Cf = 50, T0 = 0.0167, Tc = 1, Tf = 0, Tr = 0, a = 1, b = .1, d1 = 1, d2 = 1, nlevels = 30, sided = "two", par = NULL, contour.plot = FALSE, call.print = TRUE, ...){
  if(!is.null(par)){
    if(!is.vector(par)) stop("par should be a vector of length 2 or 3!") else
    if(!(length(par)==3 || length(par)==2)) stop("par should be a vector of length 2 or 3!") else
    if(length(par)==3){
      opt <- optim(par,.echXbar1,lambda=lambda,delta=delta,P0=P0,P1=P1,C0=C0,C1=C1,Cr=Cr,Cf=Cf,T0=T0,Tc=Tc,Tf=Tf,Tr=Tr,a=a,b=b,d1=d1,d2=d2,sided=sided,...)
      optimum <- as.data.frame(t(c(opt$par,opt$value)))
      gn <- round(as.numeric(optimum[3])) # global optimum of n
      cost.frame <- NULL
      for(k in c(gn-1, gn, gn+1)){
        opt <- optim(par[1:2], .echXbar2,n=k,lambda=lambda,delta=delta,P0=P0,P1=P1,C0=C0,C1=C1,Cr=Cr,Cf=Cf,T0=T0,Tc=Tc,Tf=Tf,Tr=Tr,a=a,b=b,d1=d1,d2=d2,sided=sided,...)
        cost.frame <- rbind(cost.frame,c(opt$par,k,opt$value))
      }
      optimum <- cost.frame[which(cost.frame[,4]==min(cost.frame[,4])),]
      names(optimum) <- c("Optimum h","Optimum L","Optimum n","ECH")
      optXbar <- list(optimum=optimum)
    } else{
      cost.frame <- NULL
      for(k in n){
        opt <- optim(par,.echXbar2,n=k,lambda=lambda,delta=delta,P0=P0,P1=P1,C0=C0,C1=C1,Cr=Cr,Cf=Cf,T0=T0,Tc=Tc,Tf=Tf,Tr=Tr,a=a,b=b,d1=d1,d2=d2,sided=sided,...)
        cost.frame <- rbind(cost.frame,c(opt$par,k,opt$value))
      }
      optimum <- cost.frame[which(cost.frame[,4]==min(cost.frame[,4])),]
      names(optimum) <- c("Optimum h","Optimum L","Optimum n","ECH")
      colnames(cost.frame) <- c("Optimum h","Optimum L","Optimum n","ECH")
      rownames(cost.frame) <- rep("",length(n))
      optXbar <- list(optimum=optimum,cost.frame=cost.frame)
    }
  }else
  if(missing(h) && missing(L) && missing(n)){
    opt <- optim(c(1,3,5),.echXbar1,lambda=lambda,delta=delta,P0=P0,P1=P1,C0=C0,C1=C1,Cr=Cr,Cf=Cf,T0=T0,Tc=Tc,Tf=Tf,Tr=Tr,a=a,b=b,d1=d1,d2=d2,sided=sided,...)
#    cat('The optimum parameters maybe around h=', opt$par[1],' L=', opt$par[2],' n=', opt$par[3],' and the minimum ECH is about', opt$value,'\n')
    optimum <- as.data.frame(t(c(opt$par,opt$value)))
    gn <- round(as.numeric(optimum[3])) # global optimum of n
    cost.frame <- NULL
    for(k in c(gn-1, gn, gn+1)){
      opt <- optim(c(1,3),.echXbar2,n=k,delta=delta,lambda=lambda,P0=P0,P1=P1,C0=C0,C1=C1,Cr=Cr,Cf=Cf,T0=T0,Tc=Tc,Tf=Tf,Tr=Tr,a=a,b=b,d1=d1,d2=d2,sided=sided,...)
      cost.frame <- rbind(cost.frame,c(opt$par,n=k,opt$value))
    }
    optimum <- cost.frame[which(cost.frame[,4]==min(cost.frame[,4])),]
    names(optimum) <- c("Optimum h","Optimum L","Optimum n","ECH")
    optXbar <- list(optimum=optimum)
  }else
  if(missing(h) && missing(L)){
    cost.frame <- NULL
    for(k in n){
      opt <- optim(c(1,3),.echXbar2,n=k,delta=delta,lambda=lambda,P0=P0,P1=P1,C0=C0,C1=C1,Cr=Cr,Cf=Cf,T0=T0,Tc=Tc,Tf=Tf,Tr=Tr,a=a,b=b,d1=d1,d2=d2,sided=sided,...)
      cost.frame <- rbind(cost.frame,c(opt$par,n=k,opt$value))
    }
    optimum <- cost.frame[which(cost.frame[,4]==min(cost.frame[,4])),]
    names(optimum) <- c("Optimum h","Optimum L","Optimum n","ECH")
    colnames(cost.frame) <- c("Optimum h","Optimum L","Optimum n","ECH")
    rownames(cost.frame) <- rep("",length(n))
    optXbar <- list(optimum=optimum,cost.frame=cost.frame)
  }else{
    cost.frame <- NULL
    opt.mat=outer(h,L,FUN=echXbar,n=n[1],delta=delta,lambda=lambda,P0=P0,P1=P1,C0=C0,C1=C1,Cr=Cr,Cf=Cf,T0=T0,Tc=Tc,Tf=Tf,Tr=Tr,a=a,b=b,d1=d1,d2=d2,sided=sided)
    aa <- which(opt.mat==min(opt.mat),arr.ind=TRUE)
    opt.ech <- min(opt.mat)
    cost.frame <- rbind(cost.frame,c(h[aa[1,1]],L[aa[1,2]],n[1],opt.ech))
    for(k in n[-1]){
      mat=outer(h,L,FUN=echXbar,n=k,delta=delta,lambda=lambda,P0=P0,P1=P1,C0=C0,C1=C1,Cr=Cr,Cf=Cf,T0=T0,Tc=Tc,Tf=Tf,Tr=Tr,a=a,b=b,d1=d1,d2=d2,sided=sided)
      aa <- which(mat==min(mat),arr.ind=TRUE)
      cost.frame <- rbind(cost.frame,c(h[aa[1,1]],L[aa[1,2]],k, ech <- min(mat)))
      if(ech < opt.ech){
        opt.ech <- ech; opt.mat <- mat
      }
    }
    optimum <- cost.frame[which(cost.frame[,4]==min(cost.frame[,4])),]
    names(optimum) <- c("Optimum h","Optimum L","Optimum n","ECH")
    colnames(cost.frame) <- c("Optimum h","Optimum L","Optimum n","ECH")
    rownames(cost.frame) <- rep("",length(n))
    ## FIXME
    if(contour.plot){
      if(call.print){
        ca <- match.call()
        ca[[1]] <- NULL
        par(mar=c(6.1,4.1,4.1,2.1))
        contour(h,L,opt.mat,main=strwrap(paste("ecoXbar(",paste(names(ca),"=",unlist(ca),sep="",collapse=", "),")")),xlab="h",ylab="L",nlevels=nlevels)   
    } else{
      par(mar=c(6.1,4.1,2.1,2.1))
      contour(h,L,opt.mat,xlab="h",ylab="L",nlevels=nlevels)
      }
      mtext(sprintf('Opt h=%s   Opt L=%s   n=%s   ECH=%s',optimum[1], optimum[2], optimum[3],round(optimum[4],digits=4)),side=1,line=4.5)
      points(optimum[1],optimum[2],pch=3)
    }
    optXbar <- list(optimum=optimum,cost.frame=cost.frame, opt.mat=opt.mat)
  }
  ## FIXME: alpha~sided
  opth <- as.numeric(optXbar$optimum[1])
  optL <- as.numeric(optXbar$optimum[2])
  optn <- as.numeric(optXbar$optimum[3])

  if(sided == "one"){
      ARL1 <- 1/pnorm(-optL)
      ARL2 <- 1/pnorm(-optL + abs(delta*sqrt(optn)))
    }
  if(sided == "two"){
      alpha <- 2*pnorm(-optL)
      beta <- pnorm(optL - delta*sqrt(optn))-pnorm(-optL- delta*sqrt(optn))
      ARL1 <- 1/alpha
      ARL2 <- 1/(1-beta)
    }
  s <- 1/(exp(lambda*opth)-1)
  tau <- (1-(1+lambda*opth)*exp(-lambda*opth))/(lambda*(1-exp(-lambda*opth)))
  optXbar$FAR <- lambda*s/ARL1
  optXbar$ATS <- ARL2*opth - tau
  optXbar$call <- match.call()
  return(structure(optXbar,class="edcc"))
}
##' @rdname ecoXbar
##' @export
# calculate the ECH for given parameters
echXbar <- function(h,L,n,delta=2,lambda=.05,P0=NULL,P1=NULL,C0=NULL,C1=NULL,Cr=25,Cf=50,T0=0.0167,Tc=1,Tf=0,Tr=0,a=1,b=.1,d1=1,d2=1, sided = "two"){
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



##' Calculate the optimum parameters of n(sample size), h(sampling
##' interval), k(reference value) and H(decision interval) for
##' Economic Design of the CUSUM control chart. For more information about
##' the reference value see 'Details'.
##'
##' When parameter \code{par} is specified, optimization algorithms
##' would be used as default. \code{par} can be specified as:
##' \code{par = c(h, H)} where \code{h} and \code{H} are the initial
##' values of smapling interval and decision interval when \code{n} is
##' specified; or \code{par = c(h, H, n)}. Good inital values may lead
##' to good optimum results.
##' 
##' When parameters \code{h}, \code{H}, \code{n} are all undefined,
##' \code{ecoCusum} function will try to find the global optimum point
##' to minimize the ECH (Expected Cost per Hour) using optimization
##' algorithms (\code{optim} function), but in this case \code{n}
##' would not be integer. It is usually helpful for the experimenter
##' to find the region where the optimum point may exist quickly. When
##' \code{h} and \code{H} are undefined but \code{n} is given as an
##' integer vector, \code{ecoCusum} function will try to find the
##' optimum point for each \code{n} value using optimization
##' algorithms. When \code{h}, \code{H} and \code{n} are all given,
##' \code{ecoCusum} function will use a "grid method" way to calculate the
##' optimum point, that is ECH for all the combinations of the
##' parameters will be calculated. The "grid method" way is much slower
##' than using optimization algorithms, but it would be a good choice
##' when optimization algorithms fail to work well.
##' 
##' There is strong numerical and theoretical evidence that for given
##' L1, the value of L0 approaches its maximum when k(reference value)
##' is chosen mid-way the between AQL and the RQL: $k = mu0 +
##' 0.5*delta*sigma (Appl. Statist.(1974) 23, No. 3, p. 420). For this
##' reason we treat k as a constant value and optimize n, h and H.
##' For cost parameters either {P0, P1} or {C0, C1} is needed.  If P0
##' and P1 are given, they will be used first, else C0 and C1 will be
##' used. For economic design of the CUSUM chart, when \code{d1} and
##' \code{d2} are both 1, only if the difference between P0 and P1
##' keeps the same, the results are identical. If the difference
##' between C0 and C1 keeps the same, the optimum parameters are
##' almost the same but the ECH(Expected Cost per Hour) values will
##' change.

##'
##' \code{echCusum} is used to calculate the ECH (Expected Cost per
##' Hour) for one given design point.
##' 
##' @title Economic design for the CUSUM control chart
##' @param h sampling interval. It can be a numeric vector or left
##' undefined. See 'Details'
##' @param H decision interval. It can be a numeric vector or left
##' undefined. See 'Details'
##' @param n sample size. It can be an integer vector or left
##' undefined. See 'Details'
##' @param delta shift in process mean in standard deviation units
##' when assignable cause occurs (delta = |mu1 - mu0|/sigma), where
##' sigma is the standard deviation of observations; mu0 is the
##' in-control process mean; mu1 is the out-of-control process mean.
##' Default value is 2.
##' @param lambda we assume the in-control time follows a exponential
##' distribution with mean 1/lambda. Default value is 0.05.
##' @param P0 profit per hour earned by the process operating in
##' control. See 'Details'.
##' @param P1 profit per hour earned by the process operating out of
##' control 
##' @param C0 cost per hour due to nonconformities produced while the
##' process is in control.
##' @param C1 cost per hour due to nonconformities produced while the
##' process is out of control.(C1 > C0)
##' @param Cr cost for searching and repairing the assignable cause,
##' including any downtime.
##' @param Cf cost per false alarm, including the cost of searching
##' for the cause and the cost of downtime if production ceases during
##' search.
##' @param T0 time to sample and chart one item.
##' @param Tc expected time to discover the assignable cause.
##' @param Tf expected search time when false alarm occurs.
##' @param Tr expected time to repair the process.
##' @param a fixed cost per sample.
##' @param b cost per unit sampled.
##' @param d1 flag for whether production continues during searches
##' (1-yes, 0-no). Default value is 1.
##' @param d2 flag for whether production continues during repairs
##' (1-yes, 0-no). Default value is 1.
##' @param nlevels 
##' 30. It works only when \code{contour.plot} is TRUE.
##' @param sided distinguish between one-, two-sided and Crosier's
##' modified two-sided CUSUM scheme by choosing "one", "two", and
##' "Crosier", respectively. See details in \code{\link[spc]{xcusum.arl}}
##' @param par initial values for the parameters to be optimized
##' over. It can be a vector of length 2 or 3. See 'Details'
##' @param contour.plot a logical value indicating wether a contour
##' plot should be drawn. Default is FALSE.
##' @param call.print a logical value indicating whether the "call"
##' should be printed on the contour plot. Default is TRUE
##' @param ... other arguments to be passed to \code{optim} function.
##' @return The \code{ecoCusum} function returns an object of class
##' "edcc", which is a list of elements \code{optimum},
##' \code{cost.frame}, \code{FAR} and \code{ATS}. \code{optimum} is a
##' vector with the optimum parameters and the corresponding ECH
##' value; \code{cost.frame} is a dataframe with the optimum
##' parameters and the corresponding ECH values for all given
##' \code{n}(if \code{n} is not specified, \code{cost.frame} won't be
##' returned); \code{FAR} indicates the false alarm rate during the
##' in-control time, which is calculated as lambda*(average number of
##' false alarm); \code{ATS} indicates the average time to signal
##' after the occurrence of an assignable cause, calculated as h*ARL2
##' - tau, where tau is the expected time of occurrence of the
##' assignable cause given it occurs between the i-th and (i+1)st
##' samples.
##' The \code{echCusum} function returns the calculated ECH value only.
##' @seealso \code{\link{ecoXbar}}, \code{\link{ecoEwma}},
##' \code{\link[spc]{xcusum.arl}}, \code{\link[stats]{optim}},
##' \code{\link{update.edcc}}, \code{\link{contour}}
##' @references Weicheng Zhu, Changsoon Park (2013), {edcc: An R
##' Package for the Economic Design of the Control
##' Chart}. \emph{Journal of Statistical Software}, 52(9),
##' 1-24. \url{http://www.jstatsoft.org/v52/i09/}
##'
##' Lorenzen and Vance (1986). The economic design of
##' control charts: a unified approach, \emph{Technometrics}, 28. 3-10.
##' 
##' Chiu, W.K. (1974). The economic design of CUSUM charts for
##' controlling normal means, \emph{Journal of the Royal Statistical
##' Society. Series C (Applied Statistics)}, 23(3), 420-433.
##' @examples
##' #Chiu, W.K. (1974). Applied Statistics, 23, p427 Table3, row 1-2,14
##' ## LINE 1
##' ## global optimization to h, H and n, when lambda = 0.01, "Nelder-Mead" optimization algorithm doesn't work
##' #(y <- ecoCusum( P0=150,P1=50,Cr=30,d1=0,d2=0))
##' ## we can try other algorithms:
##' (y1 <- ecoCusum( P0=150,P1=50,Cr=30,d1=0,d2=0,method="BFGS"))
##' # Based on the global optimum above, we specify the range of the
##' # parameters like this
##' (yy1 <- ecoCusum( h=seq(1.3,1.45,by=.01), H=seq(.5,0.6,by=.01),n=4:6,
##' P0=150,P1=50,Cr=30,d1=0,d2=0))
##' ## LINE 2
##' (y2 <- ecoCusum( P0=150,P1=50,Cr=30,d1=0,d2=0,lambda=0.05))
##' (yy2 <- ecoCusum( h=seq(.6,0.7,by=.01), H=seq(.5,0.6,by=.01),n=3:6,
##' P0=150,P1=50,Cr=30,d1=0,d2=0,lambda=0.05))
##' contour(yy2)
##' ## LINE 14
##' (y14 <- ecoCusum(n=30,P0=150,P1=50,Cr=30,delta=0.5,d1=0,d2=0,method="L-BFGS-B"))
##' (yy14 <- ecoCusum(h=seq(2.55,2.65,by=0.01),H=seq(0.3,0.4,by=0.01),
##' n=28:30,P0=150,P1=50,Cr=30,delta=0.5,d1=0,d2=0))
##' #Douglas (2009). Statistical quality control: a modern introduction, sixth edition, p470.
##' ecoCusum(lambda=.05,P0=110,P1=10,Cr=25,Cf=50,Tr=0,Tf=0,Tc=1,T0=.0167,a=1)
##' ecoCusum(h=seq(0.75,0.85,by=.01),H=seq(.55,0.65,by=.01),n=4:6,lambda=.05,
##' P0=110,P1=10,Cr=25,Cf=50,Tr=0,Tf=0,Tc=1,T0=.0167,a=1)
##' @export
ecoCusum <- function( h, H,n,delta = 2,lambda = .01, P0 = NULL, P1 = NULL,C0 = NULL,C1 = NULL, Cr = 20, Cf = 10,T0 = 0, Tc = .1,Tf = .1,Tr = 0.2, a = .5, b = .1, d1 = 1, d2 = 1, nlevels=30, sided = "one", par=NULL, contour.plot=FALSE, call.print=TRUE,...){
    if(!is.null(par)){
      if(!is.vector(par)) stop("par should be a vector of length 2 or 3!") else
      if(!(length(par)==3 || length(par)==2)) stop("par should be a vector of length 2 or 3!") else
      if(length(par)==3){
        opt <- optim(par,.echCusum1,delta = delta,lambda = lambda, P0 = P0, P1 = P1, C0 = C0, C1 = C1, Cr = Cr, Cf = Cf,T0 = T0, Tc = Tc,Tf = Tf,Tr = Tr, a = a, b = b, d1 = d1, d2 = d2,sided = sided,...)
#      cat('The optimum parameters maybe around h=', opt$par[1],' H=', opt$par[2],' n=', opt$par[3],' and the minimum ECH is about', opt$value,'\n')
        optimum <- as.data.frame(t(c(opt$par,opt$value)))
        gn <- round(as.numeric(optimum[3])) # global optimum of n
        cost.frame <- NULL
        for(k in c(gn-1, gn, gn+1)){
          opt <- optim(par[1:2],.echCusum2,n=k,delta = delta,lambda = lambda, P0 = P0, P1 = P1, C0 = C0, C1 = C1, Cr = Cr, Cf = Cf,T0 = T0, Tc = Tc,Tf = Tf,Tr = Tr, a = a, b = b, d1 = d1, d2 = d2,sided = sided,...)
          cost.frame <- rbind(cost.frame,c(opt$par,k,opt$value))
        }
        optimum <- cost.frame[which(cost.frame[,4]==min(cost.frame[,4])),]
        names(optimum) <- c("Optimum h","Optimum H","Optimum n","ECH")
        optCusum <- list(optimum=optimum)
      } else{
      cost.frame <- NULL
      for(k in n){
        opt <- optim(par,.echCusum2,n=k,delta = delta,lambda = lambda, P0 = P0, P1 = P1, C0 = C0, C1 = C1, Cr = Cr, Cf = Cf,T0 = T0, Tc = Tc,Tf = Tf,Tr = Tr, a = a, b = b, d1 = d1, d2 = d2,sided = sided,...)
        cost.frame <- rbind(cost.frame,c(opt$par,k,opt$value))
      }
      optimum <- cost.frame[which(cost.frame[,4]==min(cost.frame[,4])),]
      names(optimum) <- c("Optimum h","Optimum H","Optimum n","ECH")
      colnames(cost.frame) <- c("Optimum h","Optimum H","Optimum n","ECH")
      rownames(cost.frame) <- rep("",length(n))
      optCusum <- list(optimum=optimum,cost.frame=cost.frame)
    }
    }else
    if(missing(h) && missing(H) && missing(n)){
      opt <- optim(c(1,1,5),.echCusum1,delta = delta,lambda = lambda, P0 = P0, P1 = P1, C0 = C0, C1 = C1, Cr = Cr, Cf = Cf,T0 = T0, Tc = Tc,Tf = Tf,Tr = Tr, a = a, b = b, d1 = d1, d2 = d2,sided = sided,...)
#    cat('The optimum parameters maybe around h=', opt$par[1],' H=', opt$par[2],' n=', opt$par[3],' and the minimum ECH is about', opt$value,'\n')
      optimum <- as.data.frame(t(c(opt$par,opt$value)))
      gn <- round(as.numeric(optimum[3])) # global optimum of n
      cost.frame <- NULL
      for(k in c(gn-1, gn, gn+1)){
        opt <- optim(c(1,1),.echCusum2,n=k,delta = delta,lambda = lambda, P0 = P0, P1 = P1, C0 = C0, C1 = C1, Cr = Cr, Cf = Cf,T0 = T0, Tc = Tc,Tf = Tf,Tr = Tr, a = a, b = b, d1 = d1, d2 = d2,sided = sided,...)
        cost.frame <- rbind(cost.frame,c(opt$par,k,opt$value))
      }
      optimum <- cost.frame[which(cost.frame[,4]==min(cost.frame[,4])),]
      names(optimum) <- c("Optimum h","Optimum H","Optimum n","ECH")
      optCusum <- list(optimum=optimum)
    }else
  if(missing(h) && missing(H)){
    cost.frame <- NULL
    for(k in n){
      opt <- optim(c(1,1),.echCusum2,n=k,delta = delta,lambda = lambda, P0 = P0, P1 = P1, C0 = C0, C1 = C1, Cr = Cr, Cf = Cf,T0 = T0, Tc = Tc,Tf = Tf,Tr = Tr, a = a, b = b, d1 = d1, d2 = d2,sided = sided,...)
      cost.frame <- rbind(cost.frame,c(opt$par,n=k,opt$value))
    }
    optimum <- cost.frame[which(cost.frame[,4]==min(cost.frame[,4])),]
    names(optimum) <- c("Optimum h","Optimum H","Optimum n","ECH")
    colnames(cost.frame) <- c("Optimum h","Optimum H","Optimum n","ECH")
    rownames(cost.frame) <- rep("",length(n))
    optCusum <- list(optimum=optimum,cost.frame=cost.frame)
  }else{
  cost.frame <- NULL
  opt.mat <- matrix(NA,length(h),length(H))
  for(i in 1:length(h))
    for(j in 1:length(H))
      opt.mat[i,j] <- echCusum(h[i],H[j],n[1],delta = delta,lambda = lambda, P0 = P0, P1 = P1, C0 = C0, C1 = C1, Cr = Cr, Cf = Cf,T0 = T0, Tc = Tc,Tf = Tf,Tr = Tr, a = a, b = b, d1 = d1, d2 = d2,sided = sided)
  aa <- which(opt.mat==min(opt.mat),arr.ind=TRUE)
  opt.ech <- min(opt.mat)
  cost.frame <- rbind(cost.frame,c(h[aa[1,][1]],H[aa[1,][2]],n[1],opt.ech))
  for(nk in n[-1]){
    mat <- matrix(NA,length(h),length(H))
    for(i in 1:length(h))
      for(j in 1:length(H))
        mat[i,j] <- echCusum(h[i],H[j],nk,delta = delta,lambda = lambda, P0 = P0, P1 = P1, C0 = C0, C1 = C1, Cr = Cr, Cf = Cf,T0 = T0, Tc = Tc,Tf = Tf,Tr = Tr, a = a, b = b, d1 = d1, d2 = d2,sided = sided)
    aa <- which(mat==min(mat),arr.ind=TRUE)
    cost.frame <- rbind(cost.frame,c(h[aa[1,][1]],H[aa[1,][2]],nk, ech <- min(mat)))
    if(ech < opt.ech){
      opt.ech <- ech; opt.mat <- mat
    }
  }
  optimum <- cost.frame[which(cost.frame[,4]==min(cost.frame[,4])),]
  names(optimum) <- c("Optimum h","Optimum H","Optimum n","ECH")
  colnames(cost.frame) <- c("Optimum h","Optimum H","Optimum n","ECH")
  rownames(cost.frame) <- rep("",length(n))
  if(contour.plot){
    if(call.print){
      ca <- match.call()
      ca[[1]] <- NULL
    par(mar=c(6.1,4.1,4.1,2.1))
    contour(h,H,opt.mat,main=strwrap(paste("ecoCusum(",paste(names(ca),"=",unlist(ca),sep="",collapse=", "),")")),xlab="h",ylab="H",nlevels=nlevels)
    } else{
      par(mar=c(6.1,4.1,2.1,2.1))
      contour(h,H,opt.mat,xlab="h",ylab="H",nlevels=nlevels)
    }
    mtext(sprintf('Opt h=%s   Opt H=%s   n=%s   ECH=%s',optimum[1], optimum[2], optimum[3],round(optimum[4],digits=4)),side=1,line=4.5)
    points(optimum[1],optimum[2],pch=3)
  }
  optCusum <- list(optimum=optimum, cost.frame=cost.frame, opt.mat=opt.mat)
}
    opth <- as.numeric(optCusum$optimum[1])
    optH <- as.numeric(optCusum$optimum[2])
    optn <- as.numeric(optCusum$optimum[3])
    delta.std <- sqrt(as.numeric(optn))*delta
    ARL1 <- as.numeric(xcusum.arl(0.5*delta.std, optH, 0, sided=sided))
    ARL2 <- as.numeric(xcusum.arl(0.5*delta.std, optH, delta.std, sided=sided))
    s <- 1/(exp(lambda*opth)-1)
    tau <- (1-(1+lambda*opth)*exp(-lambda*opth))/(lambda*(1-exp(-lambda*opth)))
    optCusum$FAR <- lambda*s/ARL1
    optCusum$ATS <- ARL2*opth - tau
    optCusum$call <- match.call()
    return(structure(optCusum,class="edcc"))
}

##' @rdname ecoCusum
##' @export
# calculate the ECH for given parameters
echCusum <- function(h,H,n,delta = 2,lambda = .01, P0 = NULL, P1 = NULL,C0 = NULL,C1 = NULL, Cr = 20, Cf = 10,T0 = 0, Tc = .1,Tf = .1,Tr = 0.2, a = .5, b = .1, d1 = 1, d2 = 1,sided = "one"){
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

##' Calculate the optimum parameters, n(sample size), h(sampling
##' interval), w(weight to the present sample) and k(number of
##' s.d. from control limits to center line) for economic Design of
##' the EWMA control chart .
##' 
##' Parameter \code{w} should always be given, because the range of
##' \code{w} is so restricted that optimization algorithms usually
##' don't converge.
##' 
##' When parameter \code{par} is specified, optimization algorithms
##' would be used as default. \code{par} can be specified as:
##' \code{par = c(h, k)} where \code{h} and \code{k} are the initial
##' values of smapling interval and control limit when \code{n} is
##' specified; or \code{par = c(h, k, n)}. Good inital values may lead
##' to good optimum results.
##' 
##' When parameters \code{h}, \code{k}, \code{n} are
##' all undefined, \code{ecoEwma} function will try to find the global
##' optimum point to minimize the ECH (Expected Cost per Hour) using
##' optimization algorithms (\code{optim} function), but in this case
##' \code{n} would not be integer. It is usually helpful for the
##' experimenter to find the region where the optimum point may exist
##' quickly. When \code{h} and \code{k} are undefined but \code{n} is
##' given as an integer vector, \code{ecoEwma} function will try to
##' find the optimum point for each \code{n} value using optimization
##' algorithms. When \code{h}, \code{k} and \code{n} are all given,
##' \code{ecoEwma} function will use a "grid method" way to calculate the
##' optimum point, that is ECH for all the combinations of the
##' parameters will be calculated. The "grid method" way is much slower
##' than using optimization algorithms, but it would be a good choice
##' when optimization algorithms fail to work well.
##' 
##' For cost parameters either {P0, P1} or {C0, C1} is needed.  If P0
##' and P1 are given, they will be used first, else C0 and C1 will be
##' used.  For economic design of the EWMA chart, when \code{d1} and
##' \code{d2} are both 1, only if the difference between P0 and P1
##' keeps the same, the results are identical. If the difference
##' between C0 and C1 keeps the same, the optimum parameters are
##' almost the same but the ECH(Expected Cost per Hour) values will
##' change.

##'
##' \code{echEwma} is used to calculate the ECH (Expected Cost per
##' Hour) for one given design point.
##' 
##' @title Economic design for the EWMA control chart
##' @param h sampling interval. It can be a numeric vector or left
##' undefined. See 'Details'
##' @param w the weight value between 0 and 1 given to the latest
##' sample. It must be specified.
##' @param k control limit coefficient. It can be a numeric vector or
##' left undefined. See 'Details'
##' @param n sample size. It can be an integer vector or left
##' undefined. See 'Details'
##' @param delta shift in process mean in standard deviation units
##' when assignable cause occurs (delta = |mu1 - mu0|/sigma), where
##' sigma is the standard deviation of observations; mu0 is the
##' in-control process mean; mu1 is the out-of-control process mean.
##' Default value is 2.
##' @param lambda we assume the in-control time follows a exponential
##' distribution with mean 1/lambda. Default value is 0.05.
##' @param P0 profit per hour earned by the process operating in
##' control. See 'Details'.
##' @param P1 profit per hour earned by the process operating out of
##' control 
##' @param C0 cost per hour due to nonconformities produced while the
##' process is in control.
##' @param C1 cost per hour due to nonconformities produced while the
##' process is out of control.(C1 > C0)
##' @param Cr cost for searching and repairing the assignable cause,
##' including any downtime.
##' @param Cf cost per false alarm, including the cost of searching
##' for the cause and the cost of downtime if production ceases during
##' search.
##' @param T0 time to sample and chart one item. 
##' @param Tc expected time to discover the assignable cause.
##' @param Tf expected search time when false alarm occurs. 
##' @param Tr expected time to repair the process.
##' @param a fixed cost per sample.
##' @param b cost per unit sampled.
##' @param d1 flag for whether production continues during searches
##' (1-yes, 0-no). Default value is 1.
##' @param d2 flag for whether production continues during repairs
##' (1-yes, 0-no). Default value is 1.
##' @param nlevels number of contour levels desired. Default value is
##' 30. It works only when \code{contour.plot} is TRUE.
##' @param sided distinguish between one- and two-sided EWMA control
##' chart by choosing "one" and "two", respectively. See details in
##' \code{\link[spc]{xewma.arl}}
##' @param par initial values for the parameters to be optimized
##' over. It can be a vector of length 2 or 3. See 'Details'
##' @param contour.plot a logical value indicating whether a contour
##' plot should be drawn. Default is FALSE.
##' @param call.print a logical value indicating whether the "call"
##' should be printed on the contour plot. Default is TRUE
##' @param ... other arguments to be passed to contour function.
##' @return The \code{ecoEwma} function returns an object of class
##' "edcc", which is a list of elements \code{optimum},
##' \code{cost.frame}, \code{FAR} and \code{ATS}. \code{optimum} is a
##' vector with the optimum parameters and the corresponding ECH
##' value; \code{cost.frame} is a dataframe with the optimum
##' parameters and the corresponding ECH values for all given
##' \code{n}(if \code{n} is not specified, \code{cost.frame} won't be
##' returned); \code{FAR} indicates the false alarm rate during the
##' in-control time, which is calculated as lambda*(average number of
##' false alarm); \code{ATS} indicates the average time to signal
##' after the occurrence of an assignable cause, calculated as h*ARL2
##' - tau, where tau is the expected time of occurrence of the
##' assignable cause given it occurs between the i-th and (i+1)st
##' samples.
##' The \code{echEwma} function returns the calculated ECH value only.
##' @seealso \code{\link{ecoXbar}}, \code{\link{ecoCusum}},
##' \code{\link[spc]{xewma.arl}}, \code{\link{update.edcc}},
##' \code{\link[stats]{optim}},\code{\link{contour}}
##' @references Weicheng Zhu, Changsoon Park (2013), {edcc: An R
##' Package for the Economic Design of the Control
##' Chart}. \emph{Journal of Statistical Software}, 52(9),
##' 1-24. \url{http://www.jstatsoft.org/v52/i09/}
##'
##' Lorenzen and Vance (1986). The economic design of
##' control charts: a unified approach, \emph{Technometrics}, 28. 3-10.
##' @examples
##' #Douglas (2009). Statistical quality control: a modern introduction, sixth edition, p470.
##' ## Set w from 0.1 to 1 by 0.1 to catch the trend.
##' ecoEwma(w=seq(0.1,1,by=0.1),P0=110,P1=10,Cf=50)
##' #yy = ecoEwma(h = seq(.7,1,by=.01), w = seq(0.8,1,by=.01),k = seq(2.9,3.3, by = 0.01), n = 4:5, P0 = 110, P1 = 10, Cf = 50, contour.plot = TRUE)
##' 
##' ##$optimum
##' ##Optimum h Optimum k Optimum n Optimum w       ECH 
##' ##  0.81000   2.99000   5.00000   0.95000  10.36482 
##' #contour(yy)
##' @export
ecoEwma <- function( h, w, k, n, delta = 2,lambda = .05, P0 = NULL, P1 = NULL,C0 = NULL,C1 = NULL, Cr = 25, Cf = 10,T0 = 0.0167,Tc = 1, Tf = 0, Tr = 0, a = 1, b = .1,d1=1,d2=1, nlevels=30, sided="two",par=NULL,contour.plot=FALSE,call.print=TRUE,...){
  if(!is.null(par)){
    if(!is.vector(par)) stop("par should be a vector of length 2 or 3!") else
    if(!(length(par)==2 || length(par)==3)) stop("par should be a vector of length 2 or 3!") else
    if(length(par)==3){
      cost.frame <- NULL
      for(wi in w){
        opt <- optim(par,.echEwma1,w=wi,delta = delta,lambda = lambda, P0 = P0, P1 = P1, C0 = C0, C1 = C1, Cr = Cr, Cf = Cf,T0 = T0, Tc = Tc,Tf = Tf,Tr = Tr, a = a, b = b, d1 = d1, d2 = d2,sided = sided,...)
        cost.frame <- rbind(cost.frame,c(opt$par,wi,opt$value))
      }
      optimum <- cost.frame[which(cost.frame[,5]==min(cost.frame[,5])),]
      names(optimum) <- c("Optimum h","Optimum k","Optimum n", "Optimum w","ECH")
      colnames(cost.frame) <- c("Optimum h","Optimum k","Optimum n", "Optimum w","ECH")
      rownames(cost.frame) <- rep("",length(w))
      optEwma <- list(optimum=optimum, cost.frame=cost.frame)
    } else{
      cost.frame <- NULL
      for(wi in w){
        for(nk in n){
          opt <- optim(par,.echEwma2,w=wi,n=nk,delta = delta,lambda = lambda, P0 = P0, P1 = P1, C0 = C0, C1 = C1, Cr = Cr, Cf = Cf,T0 = T0, Tc = Tc,Tf = Tf,Tr = Tr, a = a, b = b, d1 = d1, d2 = d2,sided = sided,...)
          cost.frame <- rbind(cost.frame,c(opt$par,nk,wi,opt$value))
        }
      }
      optimum <- cost.frame[which(cost.frame[,5]==min(cost.frame[,5])),]
      names(optimum) <- c("Optimum h","Optimum k","Optimum n", "Optimum w","ECH")
      colnames(cost.frame) <- c("Optimum h","Optimum k","Optimum n","Optimum w","ECH")
      rownames(cost.frame) <- rep("",length(n)*length(w))
      optEwma <- list(optimum=optimum, cost.frame=cost.frame)
    }
  }else
   if( missing(w))
     stop("At least 'w' should be given!") else
   if(missing(h) && missing(k) && missing(n)){
     cost.frame <- NULL
     for(wi in w){
       opt <- optim(c(1,3,5),.echEwma1,w=wi,delta = delta,lambda = lambda, P0 = P0, P1 = P1, C0 = C0, C1 = C1, Cr = Cr, Cf = Cf,T0 = T0, Tc = Tc,Tf = Tf,Tr = Tr, a = a, b = b, d1 = d1, d2 = d2,sided = sided,...)
       cost.frame <- rbind(cost.frame,c(opt$par,wi,opt$value))
     }
     optimum <- cost.frame[which(cost.frame[,5]==min(cost.frame[,5])),]
     names(optimum) <- c("Optimum h","Optimum k","Optimum n", "Optimum w","ECH")
     colnames(cost.frame) <- c("Optimum h","Optimum k","Optimum n", "Optimum w","ECH")
     rownames(cost.frame) <- rep("",length(w))
     optEwma <- list(optimum=optimum, cost.frame=cost.frame)
   }else
   if(missing(h) && missing(k)){
     cost.frame <- NULL
     for(wi in w){
       for(nk in n){
         opt <- optim(c(1,3),.echEwma2,w=wi,n=nk,delta = delta,lambda = lambda, P0 = P0, P1 = P1, C0 = C0, C1 = C1, Cr = Cr, Cf = Cf,T0 = T0, Tc = Tc,Tf = Tf,Tr = Tr, a = a, b = b, d1 = d1, d2 = d2,sided = sided,...)
      cost.frame <- rbind(cost.frame,c(opt$par,nk,wi,opt$value))
      }
    }
     optimum <- cost.frame[which(cost.frame[,5]==min(cost.frame[,5])),]
     names(optimum) <- c("Optimum h","Optimum k","Optimum n", "Optimum w","ECH")
    colnames(cost.frame) <- c("Optimum h","Optimum k","Optimum n","Optimum w","ECH")
    rownames(cost.frame) <- rep("",length(n)*length(w))
    optEwma <- list(optimum=optimum,cost.frame=cost.frame)
  }else{
   cost.frame <- NULL
   opt.mat <- array(NA,c(length(h),length(w),length(k)))
   for(i in 1:length(k))
     for(j in 1:length(w))
       opt.mat[,j,i] = sapply(h,FUN=echEwma,w=w[j],k=k[i],n=n[1],delta = delta,lambda = lambda, P0 = P0, P1 = P1,C0 = C0,C1 = C1, Cr = Cr, Cf = Cf,T0 = T0,Tc = Tc, Tf = Tf, Tr = Tr, a = a, b = b,d1=d1,d2=d2,sided=sided)
   aa <- which(opt.mat==min(opt.mat),arr.ind=TRUE)
   opt.ech <- min(opt.mat)
   cost.frame <- rbind(cost.frame,c(h[aa[1,1]],k[aa[1,3]],n[1],w[aa[1,2]],min(opt.mat)))
   for(nk in n[-1]){
     mat <- array(NA,c(length(h),length(w),length(k)))
     for(i in 1:length(k))
       for(j in 1:length(w))
         mat[,j,i] = sapply(h,FUN=echEwma,w=w[j],k=k[i],n=nk,delta = delta,lambda = lambda, P0 = P0, P1 = P1,C0 = C0,C1 = C1, Cr = Cr, Cf = Cf,T0 = T0,Tc = Tc, Tf = Tf, Tr = Tr, a = a, b = b,d1=d1,d2=d2,sided=sided)
     aa <- which(mat==min(mat),arr.ind=TRUE)
     cost.frame <- rbind(cost.frame,c(h[aa[1,1]],k[aa[1,3]],nk, w[aa[1,2]],ech <- min(mat)))
     if(ech < opt.ech){
       opt.ech <- ech; opt.mat <- mat
     }
   }
   optimum <- cost.frame[which(cost.frame[,5]==min(cost.frame[,5])),]
   names(optimum) <- c("Optimum h","Optimum k","Optimum n", "Optimum w","ECH")
   colnames(cost.frame) <- c("Optimum h","Optimum k","n","Optimum w","ECH")
   rownames(cost.frame) <- rep("",length(n))
   mat.hk <- opt.mat[,which(w==optimum[4]),]
   mat.hw <- opt.mat[,,which(k==optimum[2])]
   mat.wk <- opt.mat[which(h==optimum[1]),,]
   if(contour.plot){
     if(call.print){
       ca <- match.call()
       ca[[1]] <- NULL
       par(mfrow=c(2,2),mar=c(5.1,4.1,1.1,1.1),oma=c(2,0,3,0))
       contour(h,k,mat.hk,xlab="h",ylab="k",nlevels=nlevels)
       points(optimum[1],optimum[3],pch=3)
       contour(h,w,mat.hw,xlab="h",ylab="w",nlevels=nlevels)
       points(optimum[1],optimum[2],pch=3)
       contour(w,k,mat.wk,xlab="w",ylab="k",nlevels=nlevels)
       points(optimum[2],optimum[3],pch=3)
       mtext(sprintf('Opt h=%s   Opt w=%s   Opt k=%s   n=%s   ECH=%s',optimum[1], optimum[4], optimum[2], optimum[3], round(optimum[5], digits=4)),outer=T,line=.3,side=1)
       title(strwrap(paste("ecoEwma(",paste(names(ca),"=",unlist(ca),sep="",collapse=", "),")")),outer=T)
     } else{
       par(mfrow=c(2,2),mar=c(5.1,4.1,1.1,1.1),oma=c(2,0,1,0))
       contour(h,k,mat.hk,xlab="h",ylab="k",nlevels=nlevels)
       points(optimum[1],optimum[2],pch=3)
       contour(h,w,mat.hw,xlab="h",ylab="w",nlevels=nlevels)
       points(optimum[1],optimum[4],pch=3)
       contour(w,k,mat.wk,xlab="w",ylab="k",nlevels=nlevels)
       points(optimum[4],optimum[2],pch=3)
       mtext(sprintf('Opt h=%s   Opt w=%s   Opt k=%s   n=%s   ECH=%s',optimum[1], optimum[4], optimum[2],optimum[3],round(optimum[5],digits=4)),outer=T,line=.3,side=1)
     }
   }
   optEwma <- list(optimum=optimum,cost.frame=cost.frame,opt.mat=list(mat.hk=mat.hk,mat.hw=mat.hw,mat.wk=mat.wk))
 }
  opth <- as.numeric(optEwma$optimum[1])
  optk <- as.numeric(optEwma$optimum[2])
  optn <- as.numeric(optEwma$optimum[3])
  optw <- as.numeric(optEwma$optimum[4])
  delta.std <- sqrt(as.numeric(optn))*delta
  ARL1 <- as.numeric(xewma.arl(optw,optk,0,sided=sided))
  ARL2 <- as.numeric(xewma.arl(optw,optk,delta.std,sided=sided))
  s <- 1/(exp(lambda*opth)-1)
  tau <- (1-(1+lambda*opth)*exp(-lambda*opth))/(lambda*(1-exp(-lambda*opth)))
  optEwma$FAR <-  lambda*s/ARL1
  optEwma$ATS <- ARL2*opth - tau
  optEwma$call <- match.call()
  return(structure(optEwma,class="edcc"))
}


##' @rdname ecoEwma
##' @export
# calculate the ECH for given parameters
echEwma <- function(h,w,k,n,delta = 2,lambda = .05, P0 = NULL, P1 = NULL,C0 = NULL,C1 = NULL, Cr = 25, Cf = 10,T0 = 0.0167,Tc = 1, Tf = 0, Tr = 0, a = 1, b = .1,d1=1,d2=1,sided="two"){
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

##' @S3method print edcc
##' @method print edcc
print.edcc <- function(x,...){
  xx <- x
  xx$opt.mat <- NULL
  xx$call <- NULL
  class(xx) <- NULL
  print.default(xx,...)
}

##' 'update' will update and (by default) re-fit a model.  It does
##' this by extracting the call stored in the object, updating the
##' call and (by default) evaluating that call. 
##'
##' S3 method for update.
##' @title Update for an "edcc" class object
##' @param object an object of "edcc" class
##' @param ... additional arguments to the call, or arguments with
##'anged values.
##' @param evaluate If true evaluate the new call else return the
##' call.
##' @return the fitted object
##' @S3method update edcc
##' @method update edcc
##' @examples x <- ecoXbar(P0=110,P1=10)
##' update(x,P0=NULL,P1=NULL,C0=10,C1=110)
update.edcc <- function(object,...,evaluate=TRUE){
    if (is.null(call <- object$call)) 
        stop("need an object with call component")
    extras <- match.call(expand.dots = FALSE)$...
    if (length(extras)) {
        existing <- !is.na(match(names(extras), names(call)))
        for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
        if (any(!existing)) {
            call <- c(as.list(call), extras[!existing])
            call <- as.call(call)
        }
    }
    if (evaluate) 
        eval(call, parent.frame())
    else call
}


##' contour plot of an "edcc" class object
##'
##' S3 method of contour plot for "edcc" class object
##' @title Contour plot of an "edcc" class object
##' @param x an object of "edcc" class
##' @param call.print a logical value indicating whether the the R
##' command should be printed on the contour plot. Default is TRUE
##' @param ... arguments to be passed to contour plot, see
##' \code{\link[graphics]{contour}} for details
##' @return a contour plot
##' @S3method contour edcc
##' @method contour edcc
##' @seealso \code{\link{ecoXbar}}, \code{\link{ecoCusum}},
##' \code{\link{ecoEwma}}, \code{\link{update.edcc}},
##' \code{\link[graphics]{contour}}
##' @examples
##' x <- ecoXbar(h=seq(0.7,0.9,by=.01),L=seq(2.8,3.3,by=.01),n=4:6,P0=110,
##' P1=10,nlevels=50,contour.plot=TRUE)
##' contour(x,nlevels=20,lty=2,col=2,call.print=FALSE)
contour.edcc <- function(x,call.print=TRUE,...){
  if (is.null(call <- x$call)) 
    stop("need an object with call component")
  extras <- match.call(expand.dots = FALSE)$...
  if(call[[1]]==quote(ecoXbar)){
    h <- call[["h"]]
    L <- call[["L"]]
    nlevels <- call[["nlevels"]]
    if(is.null(nlevels)) nlevels <- 30
    if(missing(call.print)){
      call.print <- call[["call.print"]]
      if(is.null(call.print))
        call.print <- TRUE
    }
    opt.mat <- x$opt.mat
    call1 <- list(quote(contour.default),h,L,quote(opt.mat),nlevels=nlevels,xlab="h",ylab="L")
    if (length(extras)) {
      existing <- !is.na(match(names(extras), c("nlevels")))
      call1[["nlevels"]] <- extras[[names(extras)[existing]]]
      if (any(!existing)) {
        call1 <- c(call1, extras[!existing])
      }
    }
    call1 <- as.call(call1)
    if(call.print){
      par(mar=c(6.1,4.1,4.1,2.1))
      eval(call1)
      call[[1]] <- NULL
      title(strwrap(paste("ecoXbar(",paste(names(call),"=",unlist(call),sep="",collapse=", "),")")))
    } else{
      par(mar=c(6.1,4.1,2.1,2.1))
      eval(call1)
    }
      optimum <- x$optimum
      mtext(sprintf('Opt h=%s   Opt L=%s   n=%s   ECH=%s',optimum[1], optimum[2], optimum[3],round(optimum[4],digits=4)),side=1,line=4.5)
    points(optimum[1],optimum[2],pch=3)
  }else
  if(call[[1]]==quote(ecoCusum)){
      h <- call[["h"]]
      H <- call[["H"]]
      nlevels <- call[["nlevels"]]
      if(is.null(nlevels)) nlevels <- 30
      if(missing(call.print)){
        call.print <- call[["call.print"]]
        if(is.null(call.print))
          call.print <- TRUE
      }
      opt.mat <- x$opt.mat
      call1 <- list(quote(contour.default),h,H,quote(opt.mat),nlevels=nlevels,xlab="h",ylab="H")
      if (length(extras)) {
        existing <- !is.na(match(names(extras), c("nlevels")))
        call1[["nlevels"]] <- extras[[names(extras)[existing]]]
        if (any(!existing)) {
          call1 <- c(call1, extras[!existing])
        }
      }
      call1 <- as.call(call1)
      if(call.print){
        par(mar=c(6.1,4.1,4.1,2.1))
        eval(call1)
        call[[1]] <- NULL
        title(strwrap(paste("ecoCusum(",paste(names(call),"=",unlist(call),sep="",collapse=", "),")")))
      } else{
        par(mar=c(6.1,4.1,2.1,2.1))
        eval(call1)
      }
      optimum <- x$optimum
      mtext(sprintf('Opt h=%s   Opt H=%s   n=%s   ECH=%s',optimum[1], optimum[2], optimum[3],round(optimum[4],digits=4)),side=1,line=4.5)
      points(optimum[1],optimum[2],pch=3)
    }else
  if(call[[1]]==quote(ecoEwma)){
    h <- call[["h"]]
    w <- call[["w"]]
    k <- call[["k"]]
    nlevels <- call[["nlevels"]]
    if(is.null(nlevels)) nlevels <- 30
    if(missing(call.print)){
      call.print <- call[["call.print"]]
      if(is.null(call.print))
        call.print <- TRUE
    }
    opt.mat <- x$opt.mat
    call1 <- list(quote(contour.default),h,k,quote(opt.mat$mat.hk),nlevels=nlevels,xlab="h",ylab="k")
    call2 <- list(quote(contour.default),h,w,quote(opt.mat$mat.hw),nlevels=nlevels,xlab="h",ylab="w")
    call3 <- list(quote(contour.default),w,k,quote(opt.mat$mat.wk),nlevels=nlevels,xlab="w",ylab="k")
    if (length(extras)) {
      existing <- !is.na(match(names(extras), c("nlevels")))
      call1[["nlevels"]] <- extras[[names(extras)[existing]]]
      call2[["nlevels"]] <- extras[[names(extras)[existing]]]
      call3[["nlevels"]] <- extras[[names(extras)[existing]]]
      if (any(!existing)) {
        call1 <- c(call1, extras[!existing])
        call2 <- c(call2, extras[!existing])
        call3 <- c(call3, extras[!existing])
      }
    }
    call1 <- as.call(call1)
    call2 <- as.call(call2)
    call3 <- as.call(call3)
    optimum <- x$optimum
    if(call.print){
      par(mfrow=c(2,2),mar=c(5.1,4.1,1.1,1.1),oma=c(2,0,3,0))
      eval(call1)
      points(optimum[1],optimum[3],pch=3)
      eval(call2)
      points(optimum[1],optimum[2],pch=3)
      eval(call3)
      points(optimum[2],optimum[3],pch=3)
      call[[1]] <- NULL
      title(strwrap(paste("ecoEwma(",paste(names(call),"=",unlist(call),sep="",collapse=", "),")")),outer=T)
    } else{
      par(mfrow=c(2,2),mar=c(5.1,4.1,1.1,1.1),oma=c(2,0,1,0))
      eval(call1)
      points(optimum[1],optimum[3],pch=3)
      eval(call2)
      points(optimum[1],optimum[2],pch=3)
      eval(call3)
      points(optimum[2],optimum[3],pch=3)
    }
    mtext(sprintf('Opt h=%s   Opt w=%s   Opt k=%s   n=%s   ECH=%s',optimum[1], optimum[2], optimum[3],optimum[4],round(optimum[5],digits=4)),outer=T,line=.3,side=1)
  }
}
