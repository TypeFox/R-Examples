# Functions for the book "Practices of Medical and Health Data Analysis using R"
# written by Minato Nakazawa, 2007-2014.
# rev. 0.1, 29 Mar 2010
# rev. 0.2, 24 Aug 2010, combined with demogjpn.R
# rev. 0.2.1, 7 May 2011, fix the exceptional treatment of radarchart() concerning "left"
# rev. 0.2.2, 11 Dec 2011, mhchart function was added.
# rev. 0.3.2, 6 Feb 2012, fix the exceptional treatment of radarchart() for too many NA's
# rev. 0.3.4, 27 Apr 2012, add new axistype options of radarchart()
# rev. 0.3.6, 8 Jan 2013, add new pdensity and pfcol options of radarchart()
# rev. 0.3.7, 7 Feb 2013, add a new centerzero option of radarchart()
# rev. 0.4.3, 27 January 2014, bug fix of pvalueplot().  Important!!
# rev. 0.4.4, 3 May 2014, bug fix of and adding an option to radarchart().
# rev. 0.5.0, 4 August 2014, label size option was added to radarchart().
# rev. 0.5.1, 15 September 2014, Jvital data was updated and Jvital2013byPref was added.

SIQR <- function(X, mode=1) { 
 if (mode==1) { ret <- (fivenum(X)[4]-fivenum(X)[2])/2 }
 else { ret <- IQR(X)/2 } # equals to (quantile(X)[4]-quantile(X)[2])/2
 return (ret)
}

truemedian <- function(X, h=1) { 
   YY <- rep(0,length(X))
   XX <- table(X)
   q <- length(XX)
   k <- 0
   for (i in 1:q) {
      L <- as.numeric(names(XX)[i])-h/2
      for (j in 1:XX[[i]]) {
         k <- k+1
         YY[k] <- L+h*(2*j-1)/(2*XX[[i]])
      }
   }
   return (median(YY))
}

geary.test <- function(X) {
 m.X <- mean(X)
 l.X <- length(X)
 G <- sum(abs(X-m.X))/sqrt(l.X*sum((X-m.X)^2))
 p <- pnorm((G-sqrt(2/pi))/sqrt(1-3/pi)*sqrt(l.X))
 RVAL <- list(statistic = c(G = G), p.value = p, method = "Geary's normality test", 
        data.name = deparse(substitute(X)))
    class(RVAL) <- "htest"
    return(RVAL)
}

gstem <- function (X, scale=1) {
 .stem.out <- capture.output(stem(X,scale))
 .stem.len <- length(.stem.out)
 plot(c(1,2), c(1,.stem.len), type="n", axes=FALSE, xlab="", ylab="")
 text(rep(1,.stem.len), .stem.len:1, .stem.out, pos=4)
}

ratedifference <- function(a, b, PT1, PT0, CRC=FALSE, conf.level=0.95) {
 .M <- a+b
 .T <- PT1+PT0
 .IR1 <- a/PT1
 .IR0 <- b/PT0
 .IRT <- .M/.T
 norm.pp <- qnorm(1-(1-conf.level)/2)
 if (CRC) {
  .IRC1 <- norm.pp*sqrt(a/(PT1*PT1))
  .IRC0 <- norm.pp*sqrt(b/(PT0*PT0))
  .IRCT <- norm.pp*sqrt(.M/(.T*.T))
  .MAT <- matrix(c(a, b, .M, PT1, PT0, .T, .IR1, .IR0, .IRT, 
   .IR1-.IRC1, .IR0-.IRC0, .IRT-.IRCT, .IR1+.IRC1, .IR0+.IRC0, .IRT+.IRCT), 3, 5)
  colnames(.MAT) <- c("Cases","Person-time","Incidence rates","Lower CL","Upper CL")
  rownames(.MAT) <- c("Exposed","Unexposed","Total")
 } else {
  .MAT <- matrix(c(a, b, .M, PT1, PT0, .T, .IR1, .IR0, .IRT), 3, 3)
  colnames(.MAT) <- c("Cases","Person-time","Incidence rates")
  rownames(.MAT) <- c("Exposed","Unexposed","Total")
 }
 class(.MAT) <- "table"
 print(.MAT)
 ESTIMATE <- .IR1-.IR0
 .CHI <- ESTIMATE/sqrt(a/(PT1*PT1)+b/(PT0*PT0))
 p.v <- 2*(1-pnorm(abs(.CHI)))
 RDL <- ESTIMATE-norm.pp*sqrt(a/(PT1*PT1)+b/(PT0*PT0))
 RDU <- ESTIMATE+norm.pp*sqrt(a/(PT1*PT1)+b/(PT0*PT0))
 CINT <- c(RDL,RDU)
 attr(CINT, "conf.level") <- conf.level
 RVAL <- list(p.value=p.v, conf.int=CINT, estimate=ESTIMATE, 
  method="Incidence rate difference and its significance probability (H0: The difference equals to zero)",
  data.name=paste(deparse(substitute(a)), deparse(substitute(b)),
   deparse(substitute(PT1)), deparse(substitute(PT0))))
 class(RVAL) <- "htest"
 return(RVAL)
}

rateratio <- function(a, b, PT1, PT0, conf.level=0.95) {
 .M <- a+b
 .T <- PT1+PT0
 .MAT <- matrix(c(a, b, .M, PT1, PT0, .T), 3, 2)
 colnames(.MAT) <- c("Cases","Person-time")
 rownames(.MAT) <- c("Exposed","Unexposed","Total")
 class(.MAT) <- "table"
 print(.MAT)
 ESTIMATE <- (a/PT1)/(b/PT0)
 norm.pp <- qnorm(1-(1-conf.level)/2)
 .CHI <- (a-(PT1/.T)*.M)/sqrt(.M*(PT1/.T)*(PT0/.T))
 p.v <- 2*(1-pnorm(abs(.CHI)))
 RRL <- ESTIMATE*exp(-norm.pp*sqrt(1/a+1/b))
 RRU <- ESTIMATE*exp(norm.pp*sqrt(1/a+1/b))
 CINT <- c(RRL,RRU)
 attr(CINT, "conf.level") <- conf.level
 RVAL <- list(p.value=p.v, conf.int=CINT, estimate=ESTIMATE, 
  method="Incidence rate ratio estimate and its significance probability",
  data.name=paste(deparse(substitute(a)), deparse(substitute(b)),
   deparse(substitute(PT1)), deparse(substitute(PT0))))
 class(RVAL) <- "htest"
 return(RVAL)
}

riskdifference <- function(a, b, N1, N0, CRC=FALSE, conf.level=0.95) {
 .M <- a+b
 .T <- N1+N0
 .R1 <- a/N1
 .R0 <- b/N0
 .RT <- .M/.T
 norm.pp <- qnorm(1-(1-conf.level)/2)
 if (CRC) {
  .RC1 <- norm.pp*sqrt(a*(N1-a)/(N1^3))
  .RC0 <- norm.pp*sqrt(b*(N0-b)/(N0^3))
  .RCT <- norm.pp*sqrt(.M*(.T-.M)/(.T^3))
  .MAT <- matrix(c(a, b, .M, N1, N0, .T, .R1, .R0, .RT, 
   .R1-.RC1, .R0-.RC0, .RT-.RCT, .R1+.RC1, .R0+.RC0, .RT+.RCT), 3, 5)
  colnames(.MAT) <- c("Cases","People at Risk","Risk","Lower CL","Upper CL")
  rownames(.MAT) <- c("Exposed","Unexposed","Total")
 } else {
  .MAT <- matrix(c(a, b, .M, N1, N0, .T, .R1, .R0, .RT), 3, 3)
  colnames(.MAT) <- c("Cases","People at risk","Incidence rates")
  rownames(.MAT) <- c("Exposed","Unexposed","Total")
 }
 class(.MAT) <- "table"
 print(.MAT)
 ESTIMATE <- .R1-.R0
 .CHI <- ESTIMATE/sqrt(a*(N1-a)/(N1^3)+b*(N0-b)/(N0^3))
 p.v <- 2*(1-pnorm(abs(.CHI)))
 RDL <- ESTIMATE-norm.pp*sqrt(a*(N1-a)/(N1^3)+b*(N0-b)/(N0^3))
 RDU <- ESTIMATE+norm.pp*sqrt(a*(N1-a)/(N1^3)+b*(N0-b)/(N0^3))
 CINT <- c(RDL,RDU)
 attr(CINT, "conf.level") <- conf.level
 RVAL <- list(p.value=p.v, conf.int=CINT, estimate=ESTIMATE, 
  method="Risk difference and its significance probability (H0: The difference equals to zero)",
  data.name=paste(deparse(substitute(a)), deparse(substitute(b)),
   deparse(substitute(N1)), deparse(substitute(N0))))
 class(RVAL) <- "htest"
 return(RVAL)
}

riskratio <- function(X, Y, m1, m2, conf.level=0.95) {
 .MAT <- matrix(c(X, Y, m1-X, m2-Y, m1, m2), 2, 3)
 colnames(.MAT) <- c("Disease","Nondisease","Total")
 rownames(.MAT) <- c("Exposed","Nonexposed")
 class(.MAT) <- "table"
 print(.MAT)
 ESTIMATE <- (X/m1)/(Y/m2)
 n1 <- X+Y
 Total <- m1+m2
 n2 <- Total-n1
 norm.pp <- qnorm(1-(1-conf.level)/2)
 p.v <- 2*(1-pnorm(abs((X-n1*m1/Total)/sqrt(n1*n2*m1*m2/Total/Total/(Total-1)))))
 RRL <- ESTIMATE*exp(-norm.pp*sqrt(1/X-1/m1+1/Y-1/m2))
 RRU <- ESTIMATE*exp(norm.pp*sqrt(1/X-1/m1+1/Y-1/m2))
 CINT <- c(RRL,RRU)
 attr(CINT, "conf.level") <- conf.level
 RVAL <- list(p.value=p.v, conf.int=CINT, estimate=ESTIMATE, 
  method="Risk ratio estimate and its significance probability",
  data.name=paste(deparse(substitute(X)), deparse(substitute(Y)),
   deparse(substitute(m1)), deparse(substitute(m2))))
 class(RVAL) <- "htest"
 return(RVAL)
}

oddsratio <- function(a, b, c, d, conf.level=0.95) {
 .MAT <- matrix(c(a, b, M1<-a+b, c, d, M0<-c+d, N1<-a+c, N0<-b+d, Total<-a+b+c+d), 3, 3)
 colnames(.MAT) <- c("Disease","Nondisease","Total")
 rownames(.MAT) <- c("Exposed","Nonexposed","Total")
 class(.MAT) <- "table"
 print(.MAT)
 ESTIMATE <- (a*d)/(b*c)
 norm.pp <- qnorm(1-(1-conf.level)/2)
 p.v <- 2*(1-pnorm(abs((a-N1*M1/Total)/sqrt(N1*N0*M1*M0/Total/Total/(Total-1)))))
 ORL <- ESTIMATE*exp(-norm.pp*sqrt(1/a+1/b+1/c+1/d))
 ORU <- ESTIMATE*exp(norm.pp*sqrt(1/a+1/b+1/c+1/d))
 CINT <- c(ORL,ORU)
 attr(CINT, "conf.level") <- conf.level
 RVAL <- list(p.value=p.v, conf.int=CINT, estimate=ESTIMATE, 
  method="Odds ratio estimate and its significance probability",
  data.name=paste(deparse(substitute(a)), deparse(substitute(b)),
   deparse(substitute(c)), deparse(substitute(d))))
 class(RVAL) <- "htest"
 return(RVAL)
}

Kappa.test <- function(x, y=NULL, conf.level=0.95) {
 DNAME <- deparse(substitute(x))
 METHOD <- paste("Estimate Cohen's kappa statistics and test the null hypothesis ",
  "that the extent of agreement is same as random (kappa=0)",sep="")
 if (is.data.frame(x)) 
  x <- as.matrix(x)
 if (is.matrix(x)) {
  if (any(dim(x) < 2)) 
   stop("'x' must have at least 2 rows and columns")
  if (!is.numeric(x) || any(x < 0) || any(is.na(x))) 
   stop("all entries of 'x' must be nonnegative and finite")
 }
 else {
  if (is.null(y)) 
   stop("if 'x' is not a matrix, 'y' must be given")
  if (length(x) != length(y)) 
   stop("'x' and 'y' must have the same length")
  DNAME <- paste(DNAME, "and", deparse(substitute(y)))
  OK <- complete.cases(x, y)
  x <- factor(x[OK])
  y <- factor(y[OK])
  if ((nlevels(x) < 2) || (nlevels(y) < 2)) 
   stop("'x' and 'y' must have at least 2 levels")
  x <- table(x, y)
 }
 nr <- nrow(x)
 nc <- ncol(x)
 if (nr != nc) {
  stop("levels for 2 dimensions are different")
 }
 N <- sum(x)
 Po <- sum(diag(x))/N
 Pe <- sum(rowSums(x)*colSums(x)/N)/N
 kappa <- (Po-Pe)/(1-Pe)
 JUDGEMENT <- c("No agreement","Slight agreement","Fair agreement",
  "Moderate agreement","Substantial agreement","Almost perfect agreement")
 # This criterion is given by Landis JR, Koch GG (1977) Biometrics, 33: 159-174.
 judge <- JUDGEMENT[min(which(kappa<seq(0,1,0.2)))]
 seK0 <- sqrt(Pe/(N*(1-Pe)))
 seK <- sqrt(Po*(1-Po)/(N*(1-Pe)^2))
 norm.pp <- qnorm(1-(1-conf.level)/2)
 Z <- kappa/seK0
 p.v <- 1-pnorm(Z)
 kappaL <- kappa-norm.pp*seK
 kappaU <- kappa+norm.pp*seK
 CINT <- c(kappaL,kappaU)
 attr(CINT, "conf.level") <- conf.level
 RVAL <- list(statistic=c(Z=Z), estimate=kappa, conf.int=CINT, p.value=p.v, 
  method=METHOD, data.name=DNAME)
 class(RVAL) <- "htest"
 RVAL2 <- list(Result=RVAL,Judgement=judge)
 return(RVAL2)
}

roc <- function(values,iscase) {
 cutoffs <- unique(sort(values))
 cutoffs <- c(cutoffs,max(values)+1)
 NSERIES <- length(cutoffs)
 sensitivity <- rep(0,NSERIES)
 falsepositive <- rep(0,NSERIES)
 dist <- rep(0,NSERIES)
 aucp <- rep(0,NSERIES)
 DIS <- sum(iscase==1)
 HLT <- sum(iscase==0)
 for (i in 1:NSERIES) {
  cutoffpoint <- cutoffs[i]
  positives <- ifelse(values >= cutoffpoint, 1, 0)
  PinD <- sum(iscase==1 & positives==1)
  NinH <- sum(iscase==0 & positives==0)
  sensitivity[i] <- PinD/DIS
  falsepositive[i] <- 1-NinH/HLT
  dist[i] <- sqrt((PinD/DIS)^2+(NinH/HLT)^2)
  aucp[i] <- ifelse(i==1,(1-falsepositive[i])*sensitivity[i],
             (falsepositive[i-1]-falsepositive[i])*sensitivity[i])
 }
 RVAL <- list(cutoff=cutoffs,sens=sensitivity,falsepos=falsepositive,distance=dist,aucpiece=aucp)
 class(RVAL) <- "roc"
 return(RVAL)
}

print.roc <- function(x, ...) {
 cat("cutoff\tsens\t1-spec\tdist\n")
 cat(sprintf("%5.3f\t%5.3f\t%5.3f\t%5.3f\n", x[[1]],x[[2]],x[[3]],x[[4]]))
 mlcs <- paste("Fittest Cut Off:%5.3f, Sensitivity:%5.3f, Specificity:%5.3f\nAUC=%5.3f\n", ..., sep="")
 mlcc <- which.max(x[[4]])
 cat(sprintf(mlcs,x[[1]][mlcc],x[[2]][mlcc],1-x[[3]][mlcc],sum(x[[5]])))
 invisible(x)
}

plot.roc <- function(x, ...) {
 plot(x[[3]],x[[2]],type="l",lwd=2,xlab="1-specificity",ylab="sensitivity", ...)
 lines(c(0,1),c(0,1),lwd=1,lty=2)
 invisible(x)
}

CronbachAlpha <- function(X) { # X must be matrix or data.frame with more than 2 columns
 dim(X)[2]/(dim(X)[2]-1)*(1-sum(apply(X,2,var))/var(rowSums(X)))
}

NagelkerkeR2 <- function(rr) { # rr must be the result of lm/glm
 n <- nrow(rr$model)
 R2 <- (1-exp((rr$dev-rr$null)/n))/(1-exp(-rr$null/n))
 RVAL <- list(N=n, R2=R2)
 return(RVAL)
}

VIF <- function(X) { 1/(1-summary(X)$r.squared) }

percentile <- function(dat) { # convert numeric vector into percentiles
 pt1 <- quantile(dat, probs=seq(0, 1, by=0.01), type=7) # set minimum 0 percentile.
 pt2 <- unique(as.data.frame(pt1), fromLast=TRUE)
 pt3 <- rownames(pt2)
 pt4 <- as.integer(strsplit(pt3, "%"))
 datp <- pt4[as.integer(cut(dat, c(0, pt2$pt1), labels=1:length(pt3)))]
 return(datp)
}

radarchart <- function(df, axistype=0, seg=4, pty=16, pcol=1:8, plty=1:6, plwd=1,
                       pdensity=NULL, pangle=45, pfcol=NA, cglty=3, cglwd=1,
                       cglcol="navy", axislabcol="blue", title="", maxmin=TRUE,
                       na.itp=TRUE, centerzero=FALSE, vlabels=NULL, vlcex=NULL,
                       caxislabels=NULL, calcex=NULL,
                       paxislabels=NULL, palcex=NULL, ...) {
  if (!is.data.frame(df)) { cat("The data must be given as dataframe.\n"); return() }
  if ((n <- length(df))<3) { cat("The number of variables must be 3 or more.\n"); return() }
  if (maxmin==FALSE) { # when the dataframe does not include max and min as the top 2 rows.
    dfmax <- apply(df, 2, max)
    dfmin <- apply(df, 2, min)
    df <- rbind(dfmax, dfmin, df)
  }
  plot(c(-1.2, 1.2), c(-1.2, 1.2), type="n", frame.plot=FALSE, axes=FALSE, 
       xlab="", ylab="", main=title, asp=1, ...) # define x-y coordinates without any plot
  theta <- seq(90, 450, length=n+1)*pi/180
  theta <- theta[1:n]
  xx <- cos(theta)
  yy <- sin(theta)
  CGap <- ifelse(centerzero, 0, 1)
  for (i in 0:seg) { # complementary guide lines, dotted navy line by default
    polygon(xx*(i+CGap)/(seg+CGap), yy*(i+CGap)/(seg+CGap), lty=cglty, lwd=cglwd, border=cglcol)
    if (axistype==1|axistype==3) CAXISLABELS <- paste(i/seg*100,"(%)")
    if (axistype==4|axistype==5) CAXISLABELS <- sprintf("%3.2f",i/seg)
    if (!is.null(caxislabels)&(i<length(caxislabels))) CAXISLABELS <- caxislabels[i+1]
    if (axistype==1|axistype==3|axistype==4|axistype==5) {
     if (is.null(calcex)) text(-0.05, (i+CGap)/(seg+CGap), CAXISLABELS, col=axislabcol) else
     text(-0.05, (i+CGap)/(seg+CGap), CAXISLABELS, col=axislabcol, cex=calcex)
    }
  }
  if (centerzero) {
    arrows(0, 0, xx*1, yy*1, lwd=cglwd, lty=cglty, length=0, col=cglcol)
  }
  else {
    arrows(xx/(seg+CGap), yy/(seg+CGap), xx*1, yy*1, lwd=cglwd, lty=cglty, length=0, col=cglcol)
  }
  PAXISLABELS <- df[1,1:n]
  if (!is.null(paxislabels)) PAXISLABELS <- paxislabels
  if (axistype==2|axistype==3|axistype==5) {
   if (is.null(palcex)) text(xx[1:n], yy[1:n], PAXISLABELS, col=axislabcol) else
   text(xx[1:n], yy[1:n], PAXISLABELS, col=axislabcol, cex=palcex)
  }
  VLABELS <- colnames(df)
  if (!is.null(vlabels)) VLABELS <- vlabels
  if (is.null(vlcex)) text(xx*1.2, yy*1.2, VLABELS) else
  text(xx*1.2, yy*1.2, VLABELS, cex=vlcex)
  series <- length(df[[1]])
  SX <- series-2
  if (length(pty) < SX) { ptys <- rep(pty, SX) } else { ptys <- pty }
  if (length(pcol) < SX) { pcols <- rep(pcol, SX) } else { pcols <- pcol }
  if (length(plty) < SX) { pltys <- rep(plty, SX) } else { pltys <- plty }
  if (length(plwd) < SX) { plwds <- rep(plwd, SX) } else { plwds <- plwd }
  if (length(pdensity) < SX) { pdensities <- rep(pdensity, SX) } else { pdensities <- pdensity }
  if (length(pangle) < SX) { pangles <- rep(pangle, SX)} else { pangles <- pangle }
  if (length(pfcol) < SX) { pfcols <- rep(pfcol, SX) } else { pfcols <- pfcol }
  for (i in 3:series) {
    xxs <- xx
    yys <- yy
    scale <- CGap/(seg+CGap)+(df[i,]-df[2,])/(df[1,]-df[2,])*seg/(seg+CGap)
    if (sum(!is.na(df[i,]))<3) { cat(sprintf("[DATA NOT ENOUGH] at %d\n%g\n",i,df[i,])) # for too many NA's (1.2.2012)
    } else {
      for (j in 1:n) {
        if (is.na(df[i, j])) { # how to treat NA
          if (na.itp) { # treat NA using interpolation
            left <- ifelse(j>1, j-1, n)
            while (is.na(df[i, left])) {
              left <- ifelse(left>1, left-1, n)
            }
            right <- ifelse(j<n, j+1, 1)
            while (is.na(df[i, right])) {
              right <- ifelse(right<n, right+1, 1)
            }
            xxleft <- xx[left]*CGap/(seg+CGap)+xx[left]*(df[i,left]-df[2,left])/(df[1,left]-df[2,left])*seg/(seg+CGap)
            yyleft <- yy[left]*CGap/(seg+CGap)+yy[left]*(df[i,left]-df[2,left])/(df[1,left]-df[2,left])*seg/(seg+CGap)
            xxright <- xx[right]*CGap/(seg+CGap)+xx[right]*(df[i,right]-df[2,right])/(df[1,right]-df[2,right])*seg/(seg+CGap)
            yyright <- yy[right]*CGap/(seg+CGap)+yy[right]*(df[i,right]-df[2,right])/(df[1,right]-df[2,right])*seg/(seg+CGap)
            if (xxleft > xxright) {
              xxtmp <- xxleft; yytmp <- yyleft;
              xxleft <- xxright; yyleft <- yyright;
              xxright <- xxtmp; yyright <- yytmp;
            }
            xxs[j] <- xx[j]*(yyleft*xxright-yyright*xxleft)/(yy[j]*(xxright-xxleft)-xx[j]*(yyright-yyleft))
            yys[j] <- (yy[j]/xx[j])*xxs[j]
          } else { # treat NA as zero (origin)
            xxs[j] <- 0
            yys[j] <- 0
          }
        }
        else {
          xxs[j] <- xx[j]*CGap/(seg+CGap)+xx[j]*(df[i, j]-df[2, j])/(df[1, j]-df[2, j])*seg/(seg+CGap)
          yys[j] <- yy[j]*CGap/(seg+CGap)+yy[j]*(df[i, j]-df[2, j])/(df[1, j]-df[2, j])*seg/(seg+CGap)
        }
      }
      if (is.null(pdensities)) {
        polygon(xxs, yys, lty=pltys[i-2], lwd=plwds[i-2], border=pcols[i-2], col=pfcols[i-2])
      } else {
        polygon(xxs, yys, lty=pltys[i-2], lwd=plwds[i-2], border=pcols[i-2], 
         density=pdensities[i-2], angle=pangles[i-2], col=pfcols[i-2])
      }
      points(xx*scale, yy*scale, pch=ptys[i-2], col=pcols[i-2])
    }
  }
}

pvalueplot <- function(XTAB, plot.OR=FALSE, plot.log=FALSE, xrange=c(0.01,5)) {
# XTAB must be 2x2 cross table.
# ref. Rothman KJ (2002) Epidemiology: An introduction. Oxford Univ. Press
# ref. Rothman KJ (2012) Epidemiology: An introduction. 2nd Ed.  Oxford Univ. Press
# According to 2nd ed. p.156, p-value function must match with nested confidence intervals.
# Limitation: p-values less than 0.0005 are not calculated.
  x.a <- XTAB[1,1]
  x.b <- XTAB[1,2]
  x.c <- XTAB[2,1]
  x.d <- XTAB[2,2]
  x.N1 <- sum(XTAB[,1])
  x.N0 <- sum(XTAB[,2])
  cp <- c(1:9/1000, 1:9/100, 10:90/100, 0.9+1:9/100, 0.99+1:9/1000)
  cpx <- c(cp, 1, rev(cp))
  cpy <- c(cp/2, 0.5, 0.5+cp/2)
  cRR <- exp(log(x.a*x.N0/x.b/x.N1)+qnorm(cpy)*sqrt(1/x.a-1/x.N1+1/x.b-1/x.N0))
  cOR <- exp(log(x.a*x.d/x.b/x.c)+qnorm(cpy)*sqrt(1/x.a+1/x.b+1/x.c+1/x.d))
  
  if (plot.OR) { rval <- data.frame(OR=cOR, p.value=cpx) } else {
    rval <- data.frame(RR=cRR, p.value=cpx) }
  OpLog <- ifelse(plot.log, "x", "")
  plot(rval, type="l", xlim=xrange, log=OpLog)  
  return(rval)
}

pairwise.fisher.test <- function (x, n, p.adjust.method = p.adjust.methods, ...) 
{
    # contribution from Dr. Shigenobu Aoki.
    # 28 Dec. 2009
    # exact version of pairwise.prop.test()
    p.adjust.method <- match.arg(p.adjust.method)
    METHOD <- "Pairwise comparison of proportions (Fisher)"
    DNAME <- deparse(substitute(x))
    if (is.matrix(x)) {
        if (ncol(x) != 2) 
            stop("'x' must have 2 columns")
        n <- rowSums(x)
        x <- x[, 1]
    }
    else {
        DNAME <- paste(DNAME, "out of", deparse(substitute(n)))
        if (length(x) != length(n)) 
            stop("'x' and 'n' must have the same length")
    }
    OK <- complete.cases(x, n)
    x <- x[OK]
    n <- n[OK]
    if (length(x) < 2) 
        stop("too few groups")
    compare.levels <- function(i, j) {
        fisher.test(cbind(x[c(i, j)], n[c(i, j)]-x[c(i, j)]), ...)$p.value
    }
    level.names <- names(x)
    if (is.null(level.names)) 
        level.names <- seq_along(x)
    PVAL <- pairwise.table(compare.levels, level.names, p.adjust.method)
    ANS <- list(method = METHOD, data.name = DNAME, p.value = PVAL, 
        p.adjust.method = p.adjust.method)
    class(ANS) <- "pairwise.htest"
    return(ANS)
}

mhchart <- function(LIST, XLIM=c(15,45), COL="black", FILL="white", BWD=1, ...) {
# maternity history chart
# inspired by Wood JW (1994) "Dynamics of Human Reproduction", Aldine de Gruyter, New York.
 NN <- length(LIST)
 BASE <- rep(0,NN)
 names(BASE) <- names(LIST)
 YPOS <- barplot(BASE, horiz=TRUE, xlim=XLIM, ...)
 for (i in 1:NN) {
  DAT <- LIST[[i]]
  NX <- length(DAT)
  rect(DAT[1:(NX-1)],YPOS[i]-0.5,DAT[2:NX],YPOS[i]+0.5,border=COL,col=FILL,lwd=BWD)
 }
 XSEG <- seq(XLIM[1], XLIM[2], by=5)[-1]
 segments(XSEG,0,XSEG,YPOS[NN]+1,col="navy",lty=3)
}
