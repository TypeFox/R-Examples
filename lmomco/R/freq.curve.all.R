"freq.curve.all" <-
function(lmom,aslog10=FALSE,asprob=TRUE,
              no2para=FALSE,no3para=FALSE,no4para=FALSE,no5para=FALSE,
              step=FALSE,show=FALSE,
              xmin=NULL,xmax=NULL,xlim=NULL,
              ymin=NULL,ymax=NULL,ylim=NULL,
              aep4=FALSE,exp=TRUE,gam=TRUE,gev=TRUE,gld=FALSE,
              glo=TRUE,gno=TRUE,gpa=TRUE,gum=TRUE,kap=TRUE,
              nor=TRUE,pe3=TRUE,wak=TRUE,wei=TRUE,...) {

    if(! are.lmom.valid(lmom)) {
      warning("L-moments are invalid")
      return()
    }

    if(no2para) {
     exp <- FALSE
     gam <- FALSE
     gum <- FALSE
     nor <- FALSE
    }
    if(no3para) {
     gev <- FALSE
     glo <- FALSE
     gno <- FALSE
     gpa <- FALSE
     pe3 <- FALSE
     wei <- FALSE
    }
    if(no4para) {
     aep4 <- FALSE
     gld  <- FALSE
     kap  <- FALSE
    }
    if(no5para) {
     wak <- FALSE
    }

    F <- nonexceeds()
    n <- length(F)
    AEP4 <- EXP <- GAM <- GEV <- GLD <- GLO <- GNO <- vector(mode="numeric",length=n)
    GPA  <- GUM <- KAP <- NOR <- PE3 <- WAK <- WEI <- vector(mode="numeric",length=n)

    if(aep4 == TRUE) {
      if(show == TRUE) cat("4-p Asymmetric Exponential Power distribution (takes awhile)--")
      P <- paraep4(lmom,...)
      if(show == TRUE) cat("parameters--")
      if(P$type == "kap") {
         message("   For AEP4 logic and AEP4 solution does not exist, no quantiles will be computed")
      } else {
         AEP4 <- quaaep4(F,P)
         if(aslog10 == TRUE) AEP4 <- log10(AEP4)
         if(show == TRUE) cat("quantiles\n")
      }
    }
    if(exp == TRUE) {
      if(show == TRUE) cat("Exponential distribution--")
      P <- parexp(lmom)
      if(show == TRUE) cat("parameters--")
      EXP <- quaexp(F,P)
      if(aslog10 == TRUE) EXP <- log10(EXP)
      if(show == TRUE) cat("quantiles\n")
    }
    if(gam == TRUE) {
      if(show == TRUE) cat("Gamma distribution--")
      P <- pargam(lmom)
      if(show == TRUE) cat("parameters--")
      GAM <- quagam(F,P)
      if(aslog10 == TRUE) GAM <- log10(GAM)
      if(show == TRUE) cat("quantiles\n")
    }
    if(gev == TRUE) {
      if(show == TRUE) cat("Generalized Extreme Value distribution--")
      P <- pargev(lmom)
      if(show == TRUE) cat("parameters--")
      GEV <- quagev(F,P)
      if(aslog10 == TRUE) GEV <- log10(GEV)
      if(show == TRUE) cat("quantiles\n")
    }
    if(gld == TRUE) {
      if(show == TRUE) cat("Generalized Lambda distribution (takes awhile)--")
      P <- pargld(lmom,...)
      if(show == TRUE) cat("parameters--")
      GLD <- quagld(F,P)
      if(aslog10 == TRUE) GLD <- log10(GLD)
      if(show == TRUE) cat("quantiles\n")
    }
    if(glo == TRUE) {
      if(show == TRUE) cat("Generalized Logistic distribution--")
      P <- parglo(lmom)
      if(show == TRUE) cat("parameters--")
      GLO <- quaglo(F,P)
      if(aslog10 == TRUE) GLO <- log10(GLO)
      if(show == TRUE) cat("quantiles\n")
    }
    if(gno == TRUE) {
      if(show == TRUE) cat("Generalized Normal distribution--")
      P <- pargno(lmom)
      if(show == TRUE) cat("parameters--")
      GNO <- quagno(F,P)
      if(aslog10 == TRUE) GNO <- log10(GNO)
      if(show == TRUE) cat("quantiles\n")
    }
    if(gpa == TRUE) {
      if(show == TRUE) cat("Generalized Pareto distribution--")
      P <- pargpa(lmom)
      if(show == TRUE) cat("parameters--")
      GPA <- quagpa(F,P)
      if(aslog10 == TRUE) GPA <- log10(GPA)
      if(show == TRUE) cat("quantiles\n")
    }
    if(gum == TRUE) {
      if(show == TRUE) cat("Generalized Gumbel distribution--")
      P <- pargum(lmom)
      if(show == TRUE) cat("parameters--")
      GUM <- quagum(F,P)
      if(aslog10 == TRUE) GUM <- log10(GUM)
      if(show == TRUE) cat("quantiles\n")
    }
    if(kap == TRUE) {
      if(show == TRUE) cat("Kappa distribution--")
      P <- parkap(lmom)
      if(show == TRUE) cat("parameters--")
      KAP <- quakap(F,P)
      if(aslog10 == TRUE) KAP <- log10(KAP)
      if(show == TRUE) cat("quantiles\n")
    }
    if(nor == TRUE) {
      if(show == TRUE) cat("Normal distribution--")
      P <- parnor(lmom)
      if(show == TRUE) cat("parameters--")
      NOR <- quanor(F,P)
      if(aslog10 == TRUE) NOR <- log10(NOR)
      if(show == TRUE) cat("quantiles\n")
    }
    if(pe3 == TRUE) {
      if(show == TRUE) cat("Pearson Type III distribution--")
      P <- parpe3(lmom)
      if(show == TRUE) cat("parameters--")
      PE3 <- quape3(F,P)
      if(aslog10 == TRUE) PE3 <- log10(PE3)
      if(show == TRUE) cat("quantiles\n")
    }
    if(wak == TRUE) {
      if(show == TRUE) cat("Wakeby distribution--")
      P <- parwak(lmom)
      if(show == TRUE) cat("parameters--")
      if(P$ifail == 2) { # GPA instead but likely GPA is triggered elsewhere
         message("   For WAK logic a WAK solution does not exist, no quantiles will be computed")     
      } else {
         WAK <- quawak(F,P)
         if(aslog10 == TRUE) WAK <- log10(WAK)
         if(show == TRUE) cat("quantiles\n")
      }
    }
    if(wei == TRUE) {
      if(show == TRUE) cat("Weibull distribution--")
      P <- parwei(lmom)
      if(show == TRUE) cat("parameters--")
      WEI <- quawei(F,P)
      if(aslog10 == TRUE) WEI <- log10(WEI)
      if(show == TRUE) cat("quantiles\n")
    }

    if(all(AEP4 == 0)) AEP4 <- rep(NA, n)
    if(all(EXP ==  0))  EXP <- rep(NA, n)
    if(all(GAM ==  0))  GAM <- rep(NA, n)
    if(all(GEV ==  0))  GEV <- rep(NA, n)
    if(all(GLD ==  0))  GLD <- rep(NA, n)
    if(all(GLO ==  0))  GLO <- rep(NA, n)
    if(all(GNO ==  0))  GNO <- rep(NA, n)
    if(all(GPA ==  0))  GPA <- rep(NA, n)
    if(all(GUM ==  0))  GUM <- rep(NA, n)
    if(all(KAP ==  0))  KAP <- rep(NA, n)
    if(all(NOR ==  0))  NOR <- rep(NA, n)
    if(all(PE3 ==  0))  PE3 <- rep(NA, n)
    if(all(WAK ==  0))  WAK <- rep(NA, n)
    if(all(WEI ==  0))  WEI <- rep(NA, n)

    Q <- data.frame(nonexceeds = F,
                    aep4 = AEP4, exp = EXP, gam = GAM, gev = GEV, glo = GLO,
                    gld = GLD, gno = GNO, gpa = GPA, gum = GUM,
		            kap = KAP, nor = NOR, pe3 = PE3, wak = WAK, wei = WEI)
    if(show == TRUE) {
      xlab <- "NONEXCEEDANCE PROBABILITY"
      if(asprob == TRUE) {
        F <- qnorm(F)
        xlab <- "STANDARD NORMAL VARIATE"
      }

      limx <- range(F,finite=TRUE)
      if(length(xmin) == 1) limx[1] <- xmin
      if(length(xmax) == 1) limx[2] <- xmin
      if(length(xlim) == 2) limx    <- xlim

      limy <- range(AEP4,EXP,GAM,GEV,GLO,GLD,GNO,GPA,GUM,KAP,NOR,PE3,WAK,WEI,finite=TRUE)
      if(length(ymin) == 1) limy[1] <- ymin
      if(length(ymax) == 1) limy[2] <- ymin
      if(length(ylim) == 2) limy    <- ylim
      plot(F,F,type='n',xlim=c(limx[1],limx[2]), ylim=c(limy[1],limy[2]),
                        xlab=xlab, ylab="QUANTILES")
      lines(F,Q$aep4,col=6,lty=2, lwd=2) # mimic of GLD
      lines(F,Q$exp, col=2)
      lines(F,Q$gam, col=2)
      lines(F,Q$gev, col=3)
      lines(F,Q$glo, col=3)
      lines(F,Q$gld, col=4, lty=2, lwd=2) # mimic of AEP4
      lines(F,Q$gno, col=3)
      lines(F,Q$gpa, col=3)
      lines(F,Q$gum, col=2)
      lines(F,Q$kap, col=4, lwd=2)
      lines(F,Q$nor, col=1)
      lines(F,Q$pe3, col=3)
      lines(F,Q$wak, col=6)
      lines(F,Q$wei, col=3, lty=2) # same color as GEV
    }
    return(Q)
}

