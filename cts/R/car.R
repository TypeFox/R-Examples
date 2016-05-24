#.First.lib <- function(lib, pkg)
.onLoad <- function(lib, pkg)
        library.dynam("cts", pkg, lib)

car_control <- function(fty=1, n.ahead=10, trace=FALSE, ari=TRUE, vri=FALSE, vr=0, pfi="MAPS",ccv="CTES", lpv=TRUE, scc=TRUE,  nit=40, opm=1, rgm=1, req=0.5, con=1.0e-5, rpe=1.0, ivl=1.0e-2, fac=1.0e1, stl=1.0e-5, sml=1.0e2, gtl=1.0e5, kst=TRUE, fct=TRUE){
  RET <- list(fty=fty, n.ahead=n.ahead, trace=trace, ari=ari, vri=vri, vr=vr, pfi=pfi,ccv=ccv, lpv=lpv, scc=scc, nit=nit, opm=opm, rgm=rgm, req=req, con=con, rpe=rpe, ivl=ivl, fac=fac, stl=stl, sml=sml, gtl=gtl, kst=kst, fct=fct)
  class(RET) <- c("car_control")
  RET
}

car <-
function(x, y=NULL, scale=1.5, order=3, ctrl=car_control())
  {
    call <- match.call()
 fty=ctrl$fty; n.ahead=ctrl$n.ahead; trace=ctrl$trace
    ari=ctrl$ari; phi=rep(0, order); vri=ctrl$vri; vr=ctrl$vr; pfi=ctrl$pfi; ccv=ctrl$ccv; lpv=ctrl$lpv; scc=ctrl$scc; nit=ctrl$nit; opm=ctrl$opm; rgm=ctrl$rgm; req=ctrl$req; con=ctrl$con; rpe=ctrl$rpe; ivl=ctrl$ivl; fac=ctrl$fac; stl=ctrl$stl; sml=ctrl$sml; gtl=ctrl$gtl; kst=ctrl$kst; fct=ctrl$fct;
    if (NCOL(x)==2){
      tim <- x[,1]
      ser <- x[,2]
   series <- deparse(substitute(x))}
    else {
      tim <- x
      ser <- y
      series <- deparse(substitute(y))
    }
    if (length(tim)!=length(ser)) stop("Number of time and observations are not equal")
    len <- length(tim)
    if (len > 5000) stop("Number of time beyond 5000 is not implemented")
    ### corresponding to setup.f
    csz <- 0
    if (pfi=="QLFA") pfi <- 1
    if (pfi=="QLFS") pfi <- 2
    if (pfi=="DIRA") pfi <- 3
    if (pfi=="DIRS") pfi <- 4
    if (pfi=="MAPA") pfi <- 5
    if (pfi=="MAPS") pfi <- 6
    if (pfi==0) stop("Invalid PFI option:",pfi)
#    scale <- sca
    if (scale <= csz) stop("Invalid SCALE value:",scale)
    if (order <= 0 || order > 20) stop("Model order value must between 1 and 20")

    if (ari==FALSE) ari <- 0
    if (ari==TRUE){
      if (order==length(phi)) ari <- 1
      else stop("Initial phi value length not equal to the model order")}
    if (ari < 0) stop("Invalid ARI option:", ari)
    if (ari==0) phi <- rep(0,order)

    if (vri==TRUE) vri <- 1
    if (vri==FALSE) vri <- 0
    if (vri < 0) stop("Invalid VRI option:", vri)
    if (vri==0) vr <- csz
    if (vr < 0) stop("Invalid VR value:", vr)
    if (is.null(ccv)) ccv <- 0
#    if (ccv=="NULL") ccv <- 0
    if (ccv=="MNCT") ccv <- 1
    if (ccv=="CTES") ccv <- 2
    if (ccv < 0) stop("Invalid CCV option:", ccv)

    if (lpv==FALSE) lyap <- 0
    if (lpv==TRUE) lyap <- 1
    if (lyap < 0) stop("Invalid LPV option:", lpv)
    if (lyap==1) prdg <- csz
    if (prdg < 0) stop("Invalid PRDG value:", prdg)
    
    if (scc==FALSE) scc <- 0
    if (scc==TRUE) scc <- 1
    if (scc < 0) stop("Invalid SCC option:", scc)
    if (scc==0 && lyap==1)
      {
        lyap <- 0
        prdg <- 1.0e4
        if(trace)        
        cat("\nSCC=N CAUSES RESET OF LPV=Y TO LPV=N AND PRDG=1.0D4","\n")
      }
    if (n.ahead < 0 || n.ahead > 5000) stop("Invalid forecasting steps:", n.ahead)
    if (trace) cat("\nREADING OF MODEL PARAMETER PARAMETER SUCCESSFUL")
    
    if (nit < 0 || nit > 100) nit <- 25
    if (opm < 0) stop("Invalid OPM value:", opm)
    if (rgm < 0) stop("Invalid RGM value:", rgm)
    if (req < 0) stop("Invalid REQ value:", req)
    concrit <- con
    if (con < 0) stop("Invalid CONCRIT value:", concrit)
    rpert <- rpe
    if (rpe < 0) stop("Invalid RPERT value:", rpert)
    ivlam <- ivl
    if (ivlam < 0) stop("Invalid IVLAM value:", ivlam)
    if (fac < 0) stop("Invalid FAC value:", fac)
    stlam <- stl
    if (stlam < 0) stop("Invalid STLAM value:", stlam)
    smlam <- sml
    if (smlam < stlam) stop("Invalid SMLAM value:", smlam)
    gtlam <- gtl
    if (gtlam < 0) stop("Invalid GTLAM value:", gtlam)
    
    if (kst==FALSE) kst <- 0
    if (kst==TRUE) kst <- 1
    if (kst < 0) stop("Invalid KST value:", kst)

    if (fct==FALSE) fct <- 0
    if (fct==TRUE) fct <- 1
    if (fct < 0) stop("Invalid FCT value:", fct)

    if(fty < 1 || fty > 3) stop("Invalid fty value:", fty)
    tra <- 0
    if(trace){
    tra <- 1
    cat("\nREADING OF CONTROL PARAMETER SUCCESSFUL","\n")
    }
    np1 <- 0
    z <-.Fortran("setup",
                 as.integer(pfi),
                 as.integer(order),
                 as.integer(vri),
                 as.integer(ccv),
                 as.double(scale),
                 as.integer(ari),
                 as.double(vr),
                 as.double(phi),
                 as.integer(lyap),
                 as.double(prdg),
                 as.integer(scc),
                 as.integer(fct),
                 as.integer(fty),
                 as.integer(len),
                 as.integer(n.ahead),
                 as.double(tim),
                 as.double(ser),
                 as.integer(nit),
                 as.integer(opm),
                 as.integer(rgm),
                 as.double(req),
                 as.double(concrit),
                 as.double(rpert),
                 as.double(ivlam),
                 as.double(fac),
                 as.double(stlam),
                 as.double(smlam),
                 as.double(gtlam),
                 as.integer(kst),
                 np1=as.integer(np1),
                 as.integer(tra),
                 package="cts")
    if(trace)
    .Fortran("display",package="cts")
        zpar <- .Fortran("loop",
                ss=double(nit+1),
                bit=as.double(matrix(0, nit+1, 22)),
                errno=integer(1),  package="cts")
        if(zpar$errno == 1) stop("ROOT WITH POSITIVE REAL PART, WHICH WAS CALLED IN loop.f")
        else if(zpar$errno == 2) stop("ERROR IN LAPACK SUBROUTINE zgesv, WHICH WAS CALLED IN cinvert.f")
        else if(zpar$errno == 3) stop("PROGRAM FAILS IN SLICE ROUTINE LYBSC, WHICH WAS CALLED IN resg1d.f")
        else if(zpar$errno == 30) stop("PROGRAM FAILS IN SLICE ROUTINE LYBSC, WHICH WAS CALLED IN resg1dpre.f")
        else if(zpar$errno == 31) stop("PROGRAM FAILS IN SLICE ROUTINE LYBSC, WHICH WAS CALLED IN resg1dpre1.f")
        else if(zpar$errno == 4) stop("PROGRAM FAILS IN SLICE ROUTINE MEPAD, WHICH WAS CALLED IN resg1d.f")
        else if(zpar$errno == 40) stop("PROGRAM FAILS IN SLICE ROUTINE MEPAD, WHICH WAS CALLED IN resg1dpre.f")
        else if(zpar$errno == 51) stop("PROGRAM FAILS IN SUBROUTINE DPOTRF, WHICH WAS CALLED IN revg1.f")
        else if(zpar$errno == 52) stop("PROGRAM FAILS IN SUBROUTINE DPOTRI, WHICH WAS CALLED IN revg1.f")
        else if(zpar$errno == 6) stop("PROGRAM FAILS IN SUBROUTINE DIVC , WHICH WAS CALLED IN dpca.f")
        neff <- zpar$ss > 0
        zpar$tnit <- seq(0, sum(neff))
        zpar$ss <- zpar$ss[neff]
        zpar$bit <- matrix(zpar$bit, nrow=nit+1)
        zpar$bit <- zpar$bit[neff, 1:(order+1)]
        zfin <- .Fortran("complete",
                ss=double(1),
                bit=double(22),
                package="cts")
        zpar$ss <- c(zpar$ss, zfin$ss)
        zpar$bit <- rbind(zpar$bit, zfin$bit[1:(order+1)])
    Z <- .Fortran("setcom",
                  pfi1=integer(1),
                  arp1=integer(1),
                  np1=integer(1),
                  vri1=integer(1),
                  ccv1=integer(1),
                  len1=integer(1),
                  scale1=double(1),            
                  vr1=double(1),
                  sigsq1=double(1),
                  ##ESSP1 in Fortran has dimention (22,22)
                  essp1=as.double(matrix(0,22,22)),
                  ecov1=as.double(matrix(0,22,22)),
                  b1= double(z$np1),
                  delb1=double(z$np1),
                  rootr1=double(order),
                  rooti1=double(order),
                  package="cts")
    essp <- matrix(Z$essp1,byrow=FALSE,ncol=22)
    essp <- essp[1:Z$np1,1:Z$np1]
    ecov <- matrix(Z$ecov1,byrow=FALSE,ncol=22)
    ecov <- ecov[1:order,1:order]
    n <- length(tim)
  if (kst==1)
    {     
    .Fortran("kfilsm",package="cts")   
    Zkfil <- .Fortran("setkfilsm",
                      fser1=double(n),
                      fvar1=double(n),
                      sser1=double(n),
                      svar1=double(n),
                      sres1=double(n),
                      package="cts")
    filser <- Zkfil$fser1
    filvar <- Zkfil$fvar1
    sser <- Zkfil$sser1
    svar <- Zkfil$svar1
    sres <- Zkfil$sres1
  }
    else
      {
        filser <-  NULL
        filvar <- NULL
        sser <- NULL
        svar <- NULL
        sres <- NULL
      }

        if (pfi==1)
          ST1 <- 'QLFA'
        else if (pfi==2)
          ST1 <- 'QLFS'
        else if (pfi==3)
          ST1 <- 'DIRA'
        else if (pfi==4)
          ST1 <- 'DIRS'
        else if (pfi==5)
          ST1 <- 'MAPA'
    else
    ST1 <- 'MAPS'

    if (ccv==0)
      ST2 <- NULL
    else if (ccv==1) ST2 <- 'MNCT'
    else ST2 <- 'CTES'
    
    if (vri==1) ST3 <- 'N'
    else ST3 <- 'Y'

    .Fortran("update",package="cts")
    phi <- .Fortran("setupdate",
                    phi1=double(20),
    ###             phi1=double(order), ### changed to avoid array overrun, 4/4/2012
                    package="cts")$phi 
    ### I am not sure this phi is updated
    if (n.ahead > 0)
      {
        .Fortran("forecast",package="cts")
        if(fct==TRUE && fty==1)
          ntim <- len+n.ahead
        else ntim <- len
        Zfor <- .Fortran("setfor",
                         pre1=double(n.ahead),
                         prv1=double(n.ahead),
                         tim1=double(ntim),
                         pre2=double(5000),
                         prv2=double(5000),
                         package="cts")
        pre <- Zfor$pre1
        prv <- Zfor$prv1
        tim <- Zfor$tim1
        pre2 <- Zfor$pre2[Zfor$prv2 > 0]
        prv2 <- Zfor$prv2[Zfor$prv2 > 0]
      }
    else
      {
        pre <- NULL
        prv <- NULL
      }
    ntim <- length(ser); ntim1 <- length(tim)-length(pre)
    ntim2 <- length(tim)
    aic <- ntim * log(zpar$ss[length(zpar$tnit)]) + 2*Z$np1
    bic <- ntim * log(zpar$ss[length(zpar$tnit)]) + Z$np1*log(ntim)
    structure(list(call=call,series = series,order=Z$arp1,
                   np=Z$np1,scale=Z$scale1,
                   vri=Z$vri1, x.mean=ifelse(vri==1, 
                   Z$b1[order+2], Z$b1[order+1]), 
                   vob=ifelse(vri==1,Z$b1[order+1]^2, NA), vr=Z$vr1,
                   sigma2=Z$sigsq1,phi=Z$b1[1:order],
                   #sigma2=Z$sigsq1,phi=phi,
                   b=Z$b1,delb=Z$delb1,essp=essp,
                   ecov=ecov,
                   rootr=Z$rootr1,rooti=Z$rooti1,
                   tim=tim[1:ntim],ser=ser,n.used=Z$len1,
                   filser=filser,filvar=filvar,
                   sser=sser,svar=svar,
                   stdres=sres, pretime = tim[(ntim1+1):ntim2],
                   predict=pre,predict.var=prv*Z$sigsq1, pre2=pre2, prv2=prv2*Z$sigsq1, fty=fty, tnit=zpar$tnit, ss=zpar$ss, bit=zpar$bit, aic=aic, bic=bic),
              class="car")
  }

print.car <- function(x, digits = 3, ...)
{
  cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
  coef <- drop(round(x$phi, digits = digits))
  cat("\nEstimated coefficients\n")
  names(coef) <- paste("phi_",seq(length=x$order), sep="")
  print(coef)
  invisible(x)
}

summary.car <- function(object, ...)
{
  x <- object
  digits <- 3
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  cat(paste("\nOrder of model = ", x$order, ", sigma^2 = ",
      format(x$sigma2, digits = digits), sep=""),"\n")
  tmp <- coef <- drop(round(x$phi, digits = digits))
#  coef <- drop(round(c(x$phi, x$x.mean), digits = digits))
#  names(coef) <- c(seq(length=x$order), "mean")
  if(x$vri==1){
#  tmp <- rbind(coef, x$delb[-(x$order+1)])
  cat("\nObservation error variance:",format(x$vob*x$sigma2, digits=digits),"\n")
}
  tmp <- rbind(coef, x$delb[1:x$order])
  cat("\nEstimated coefficients (standard errors):\n")
  colnames(tmp) <- paste("phi_",seq(length=x$order), sep="")
  rownames(tmp) <- c("coef", "S.E.")
  print(round(tmp, digits))
  if(!is.na(x$x.mean)){
  cat("\nEstimated mean (standard error):\n")
  print(round(x$x.mean, digits=digits))
  print(round(x$delb[x$np], digits))
}
  corout <- FALSE
  if(corout){
  cat("\nCorrelations:\n")
  x$essp <- round(x$essp, digits)
  x$essp[upper.tri(x$essp)] <- ""
  tmp <- as.data.frame(x$essp)
  colnames(tmp) <- 1:dim(x$essp)[1]
  print(tmp)
  }
  invisible(x)
}

predict.car <- function(object, se.fit=TRUE, digits = 3, plot.it=TRUE, ...)
{
      ar <- object$phi
      p <- object$order
      pretime <- object$pretime
      names(pretime) <- seq(length=length(object$predict))
      pred <- drop(round(object$predict,digits = digits))
      names(pred) <- seq(length=length(object$predict))
      prv <- drop(round(object$predict.var,digits = digits))
      names(prv) <- seq(length=length(object$predict.var))
      cat("\nCall:\n", deparse(object$call), "\n\n", sep = "")
      tmp <- rbind(pretime, pred)
      rownames(tmp) <- c("Time", "Predict") 
      print(tmp)
      if(plot.it)
      plot.predict.car(object, ...) 
      RET <- list(pretime=pretime, pred=pred)
}

plot.predict.car <- function(object,xlab = "time",
    ylab = "", type = "l", main = NULL, sub = NULL,...)
  {
    if(class(object) != "car")
    stop("object must be car\n")
    if(object$fty==1){ ### forecast type is beyond of the observation
    plot(object$tim, object$ser,
         xlab=xlab, ylab=ylab, type=type, main=main,sub=sub, ...)
    points(object$pretime, object$predict, lty=2, col="red", ...)
    }
    else{
    npre <- length(object$pretime)
    ntim <- length(object$tim)
    plot(object$tim, object$ser,
         xlab=xlab,ylab=ylab,type=type,main=main,sub=sub, ...)
#    plot(object$tim[1:length(object$ser)],object$ser,ylim=c(min(object$ser)-2*max(sqrt(object$predict.var)),max(object$ser)+2*max(sqrt(object$predict.var))),xlab=xlab,ylab=ylab,type=type,main=main,sub=sub, ...)
    points(object$pretime, object$predict,
           lty=2,col="red", ...)
}
#    lines(object$tim[(length(object$tim) - length(object$predict) +
#                      1):length(object$tim)],object$predict+2*sqrt(object$predict.var),lty=2,col="blue")
#    lines(object$tim[(length(object$tim) - length(object$predict) +
#                      1):length(object$tim)],object$predict-2*sqrt(object$predict.var),lty=2,col="blue")
  }

tsdiag <- function(object, ...)
    UseMethod("tsdiag")

tsdiag.car <- function(object, gof.lag = 10, ...)
{
    ## plot standardized residuals, acf of residuals,
  ## cumulative periodogram and Ljung-Box p-values
    if(class(object) != "car")
    stop("object must be 'car' class\n") 
    oldpar<- par(mfrow = c(2, 2))
    on.exit(par(oldpar))
    stdres <- object$stdres
    plot(stdres, type = "h", main = "Standardized Residuals", ylab = "")
    abline(h = 0)
    acf(object$stdres, plot = TRUE, main = "ACF of Standardized Residuals",
        na.action = na.pass)
    cpgram(stdres,main="Cumulative periodogram")
    nlag <- gof.lag
    pval <- numeric(nlag)
    for(i in 1:nlag) pval[i] <- Box.test(stdres, i, type="Ljung-Box")$p.value
    plot(1:nlag, pval, xlab = "lag", ylab = "p value", ylim = c(0,1),
         main = "p values for Ljung-Box statistic")
    abline(h = 0.05, lty = 2, col = "blue")
    invisible(object)
}

