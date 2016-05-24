# --------------------------------------------------------------------------
# S3 method for printing the results of pa.ABE(), pa.scABE()
# --------------------------------------------------------------------------
print.pwrA <- function(x, digits=4, plotit=TRUE, ...)
{
  
  if(interactive() && plotit) plot(x)
  
  min.pwr  <- x$minpower
  CV.max   <- max(x$paCV[,"CV"])
  CVmaxI   <- which(x$paCV[,"CV"]==CV.max)
  CV.min   <- min(x$paCV[,"CV"])
  CVminI   <- which(x$paCV[,"CV"]==CV.min)
  min.pwrN <- min(x$paN[,"pwr"])
  min.N    <- min(x$paN[,"N"])

  if (abs(x$paCV[CVmaxI,"pwr"]-min.pwr)>1e-4) CV.max <- NA
  if (abs(x$paCV[CVminI,"pwr"]-min.pwr)>1e-4) CV.min <- NA
  
  if (x$plan[1,"theta0"]<=1) {
    min.theta0 <- min(x$paGMR[,"theta0"])
  } else {
    min.theta0 <- max(x$paGMR[,"theta0"])
  }
  method <- x$method
  if (method=="scABE") {
    meth <- ifelse(x$regulator=="EMA", " (EMA/scABEL)", " (FDA/RSABE)")
    meth <- ifelse(x$regulator=="ANVISA", " (ANVISA/scABEL)", " (FDA/RSABE)")
    method <- paste0(method, meth)
  }
  cat("Sample size plan ", method, "\n",sep="")
  # without "nlast"
  print(x$plan[,-ncol(x$plan)], row.names = FALSE)
  cat("\nPower analysis\n")
  cat("CV, theta0 and number of subjects which lead to min. acceptable ",
      "power of at least ", round(min.pwr, digits), ":\n", sep="")
  #react to RSABE NTID where there may be a CV.min, CV.max which
  if (method!="RSABE NTID"){
    # let to power=minpower
    cat(" CV= ", round(CV.max, digits), ", theta0= ",
        round(min.theta0, digits),"\n", sep="")
  } else {
      # let to power=minpower
      cat(" CV = (", round(CV.min, digits), ", ", round(CV.max, digits),
          "), theta0= ", round(min.theta0, digits),"\n", sep="")
  }
  cat(" N = ", min.N, " (power= ", round(min.pwrN, digits), ")\n", sep="")
  cat("\n")
}

# ----------------------------------------------------------------------------
# S3 method for plotting the results of pa.ABE(), pa.scABE(), pa.NTIDFDA()
#
# Author D.Labes, from original code by H. Schuetz
# reworked by H. Schuetz to avoid overlay of target power line with legend
# ----------------------------------------------------------------------------
plot.pwrA <- function(x, pct=TRUE, cols=c("blue", "red"), ...)
{
  # make colors between the both given
  mkColors <- function(){
    # all varibles are seen if this function is inside plot
    colno <- cut(pwr, breaks=unique(c(1*fact, targetpower, pwr[pwr<=targetpower])),
                 labels=FALSE, include.lowest=TRUE, right=FALSE)
    clr1 <-colorRampPalette(rev(cols))(length(unique(colno)))
    clr  <- clr1[colno]
    return(clr)
  }
  
  if (pct) {
    fact    <- 100
    ylabtxt <- "% power"
    dec     <- 2
    pctsign <- "%"
  }   else {
    fact    <- 1
    ylabtxt <- "power"
    dec     <- 4
    pctsign <- ""
  } 

  n <- nrow(x$paCV)
  # Attention next functiones only if last value is minpower
  minpower    <- fact*x$minpower
  CV.max      <- fact*max(x$paCV[,"CV"])
  CV.min      <- fact*min(x$paCV[,"CV"])
  
  # CV of call of pa.XXX() function
  if (x$method=="ABE") {
    CV  <- fact*x$plan[1,"CV"]
  } else {
    # Attention this functions only if CVwT==CVwR!!!
    CV  <- fact*x$plan[1,"CVwR"]
  }
  GMR         <- x$plan[1,"theta0"]
  theta1      <- x$plan[1,"theta1"]
  theta2      <- x$plan[1,"theta2"]
  n.est       <- x$plan[1,"Sample size"]
  pwr.est     <- fact*x$plan[1,"Achieved power"]
  targetpower <- fact*x$plan[1,"Target power"]
  reg         <- x$regulator
  design      <- x$plan[1,"Design"]
  
  mklegend <- function(method){
    if (method=="scABE"){
      algo <- ifelse(reg == "FDA", "RSABE", "scABEL")
      legend("topright", legend=c(algo, paste0("(",reg,")")), cex=0.90,
      bg="white", box.lty=0)
    }
    if (method=="RSABE NTID"){
      legend("topright", legend=c("RSABE", "(NTID FDA)"), cex=0.90,
      bg="white", box.lty=0)
    }
  }
  
  op <- par(no.readonly=TRUE) # save par() options
  par(mar=c(c(4, 4, 2.5, 0.75))+0.1) # default for B, L, T, R: c(5, 4, 4, 2) + 0.1
  par(cex.main=0.95, cex.axis=0.95, cex.lab=0.95, mgp=c(2,0.75,0), tcl=-0.2)
  
  # plot at a paneel with 4 pieces
  split.screen(c(2, 2))
  
  screen(1) ### 'Sensitivity' of CV (GMR and n constant) ###
  pwr <- as.numeric(fact*x$paCV[,"pwr"])
  CVs <- as.numeric(fact*x$paCV[,"CV"])
  seg <- length(pwr); s <- seq(seg-1)
  clr <- mkColors()
  xlabtxt <- "CV"
  if (fact==100) xlabtxt <- "CV %"

  if (x$method=="ABE"){
    plot(CVs, pwr, type="n",
         main=paste0("Higher variability\n", "theta0 = ", GMR, ", n = ", n.est),
         lwd=2, xlab=xlabtxt, ylab=ylabtxt, las=1)
    box()
    abline(h=c(targetpower, fact*0.8, minpower), lty=3, col="grey50")
    segments(CVs[s], pwr[s], CVs[s+1], pwr[s+1], lwd=2, col=clr[s])
    points(CVs[1], pwr[1], col=clr[1], pch=16, cex=1.25)
    points(CVs[seg], pwr[seg], col=clr[seg], pch=16, cex=1.25)
    text(CV, (minpower+(pwr.est-minpower)*0.1), 
         labels=paste0("CV = ", signif(CV.max, 4), pctsign," (", 
                       round(minpower, dec), pctsign,")"), cex=0.9, pos=4)
  } else {
    # any scABE (including RSABE NTID)
    plot(CVs, pwr, type="n",
         main=paste0("Lower/higher variability\n", "theta0 = ", GMR, ", n = ", n.est),
         lwd=2, xlab=xlabtxt, ylab=ylabtxt, las=1)
    abline(h=c(targetpower, 0.8*fact, minpower), lty=3, col="grey50")
    mklegend(x$method)
    box()
    segments(CVs[s], pwr[s], CVs[s+1], pwr[s+1], lwd=2, col=clr[s])
    # mark the plan CV and power
    points(CV, pwr.est, col=cols[1], pch=16, cex=1.25)
    # mark the max. CV if pwr for it is = minpower
    if (abs(pwr[seg]-minpower)/minpower<=1e-4){
      points(CVs[seg], pwr[seg], col=clr[seg], pch=16, cex=1.25)
    }
    txt <- paste0("CV = ", signif(CV.max, 4), pctsign," (", round(minpower, dec)
                  , pctsign,")")
    if  (x$method=="RSABE NTID"){
      if(abs(pwr[1]-minpower)/minpower<=1e-4) {
        #we have also CV.min with power=minpower
        points(CV.min, pwr[1], col=clr[seg], pch=16, cex=1.1)
        txt <- paste0("CV = (", signif(CV.min, 4),", ", signif(CV.max, 4),
                      pctsign,") (",round(minpower, dec), pctsign,")")
      }
    }
    text(min(CVs), (min(pwr)+(pwr.est-min(pwr))*0.1), labels=txt,
         cex=0.8, pos=4)
  }
  
  screen(2) ### 'Sensitivity' of GMR (CV and n constant) ###
  pwr <- as.numeric(fact*x$paGMR[,"pwr"])
  GMRs <- as.numeric(x$paGMR[,"theta0"])
  #GMR.min <- ifelse(GMR<=1, GMRs[1], GMRs[length(GMRs)])
  GMR.min <- GMRs[1] # or better min(GMRs) or max(GMRs)?
  seg <- length(pwr); s <- seq(seg-1)
  clr  <- mkColors()
  plot(GMRs, pwr, type="n", 
       main=paste0("Larger deviation of theta0 from 1\n", "CV = ",
                   CV, pctsign,", n = ", n.est), 
       lwd=2, xlim=c(GMR, GMR.min), xlab="theta0", ylab=ylabtxt, las=1,
       cex.main=0.95, cex.axis=0.95)
  abline(h=c(targetpower, fact*0.8, minpower), lty=3, col="grey50")
  mklegend(x$method)
  box()
  segments(GMRs[s], pwr[s], GMRs[s+1], pwr[s+1], lwd=2, col=clr[s])
  # the next assumes that the values start at GMR and end on GMR.min (maybe also max!)
  # TODO rework if plan.GMR not at border
  points(GMRs[1], pwr[1], col=clr[1], pch=16, cex=1.25)
  points(GMRs[seg], pwr[seg], col=clr[seg], pch=16, cex=1.25)
  text(GMR, (minpower+(pwr.est-minpower)*0.1), 
       labels=paste0("GMR = ",signif(GMR.min, 4), " (",
                     round(minpower, dec), pctsign, ")"),
       cex=0.85, pos=4)
  
  screen(3) ### Sensitivity of n (GMR and CV constant) ###
  pwr <- as.numeric(fact*x$paN[,"pwr"])
  Ns  <- as.numeric(x$paN[,"N"])
  clr <- mkColors()
  if(length(clr)==1) clr <- cols[1]
  xticks <- NULL
  nNs    <- length(Ns)
  if(nNs<5 & nNs>1) xticks <- c(max(Ns), min(Ns), nNs-1)
  plot(Ns, pwr, type="n", 
       main=paste0("Drop-outs\n", "theta0 = ", GMR, ", CV = ", CV, pctsign),
       lwd=2, xlim=c(max(Ns), min(Ns)), ylim=c(minpower, pwr.est),
       xlab="n", xaxp=xticks,
       ylab=ylabtxt, las=1, cex.main=0.95)
  abline(h=c(targetpower, fact*0.8, minpower), lty=3, col="grey50")
  mklegend(x$method)
  box()
  points(Ns, pwr, pch=16, cex=0.8, col=clr)
  points(Ns[length(Ns)], pwr[length(Ns)], col=clr[length(Ns)],
         pch=16, cex=1.25)
  points(n.est, pwr.est, col=clr[1], pch=16, cex=1.25)
  text(max(Ns), (minpower+(pwr.est-minpower)*0.1), 
       labels=paste0("n = ", min(Ns), " (", signif(min(pwr), 4), pctsign, ")"),
       cex=0.85, pos=4)
  
  screen(4) ### Some basic information ###
  if (x$method!="RSABE NTID"){
      CVtxt <- sprintf("  %s %+5.1f%%", "CV =",  100*(CV.max-CV)/CV)
  } else {
      CVtxt <- ""
      if(abs(fact*x$paCV[1,"pwr"]-minpower)/minpower<=1e-4) {
        #we have also CV.min with power=minpower
        CVtxt <- sprintf("  %s %+5.1f%%", "CVmin =", 100*(CV.min-CV)/CV)
        CVtxt <- c(CVtxt, sprintf("  %s %+5.1f%%", "CVmax =", 100*(CV.max-CV)/CV))
      } else {
        #we have only CV.max with power=minpower
        CVtxt <- sprintf("  %s %+5.1f%%", "CV =", 100*(CV.max-CV)/CV)
      }
  }
  if(x$method=="ABE"){
    BEARtxt <- paste("BE AR: ", round(theta1,4),
                     " ... ", round(theta2,4), sep="")
  }
  if(x$method=="scABE"){
    # (widened) acceptance range
    if(x$regulator=="FDA"){
      Ltxt <-"implied BE AR: "
      wtheta1 <- min(theta1,exp(CV2se(CV/fact)*log(theta1)/0.25))
      wtheta2 <- max(theta2,exp(CV2se(CV/fact)*log(theta2)/0.25))
    } else { #EMA
      Ltxt <- "(widened) BE AR: " 
      CVV <- min(0.5,CV/fact)      # cap
      wtheta1 <- min(theta1,exp(-CV2se(CVV)*0.76))
      wtheta2 <- max(theta2,exp(CV2se(CVV)*0.76))
    }
    BEARtxt <- paste0(Ltxt, round(wtheta1,4),
                     " ... ", round(wtheta2,4))
  }
  if(x$method=="RSABE NTID"){
    Ltxt <-"implied BE AR: "
    wtheta1 <- max(theta1,exp(CV2se(CV/fact)*log(0.9)/0.1))
    wtheta2 <- min(theta2,exp(-CV2se(CV/fact)*log(0.9)/0.1))
    BEARtxt <- paste0(Ltxt, round(wtheta1,4),
                      " ... ", round(wtheta2,4))
  }  
  plot(1, type="n", axes=F, xlab="", ylab="")
  if (fact==100){
    # percent
    legend("topleft", 
           legend=c(paste0(design, " design", "; assumed:"),
                    sprintf("  %s %1.0f%%%s %5.4f", "CV =", CV, ", theta0 =", GMR),
                    BEARtxt,
                    "power:",
                    sprintf("  %s %2.0f%%", "target =", targetpower),
                    sprintf("  %s %5.2f%% %s %i%s", "estimated =", pwr.est,
                            "(n =", n.est, ")"),
                    sprintf("  %s %2.0f%%", "minimum acceptable =", minpower),
                    "acceptable rel. deviations:",
                    #TODO:react to RSABE NTID where there may be also a CVmin
                    CVtxt,
                    sprintf("  %s %+5.2f%%", "GMR =", 100*(GMR.min-GMR)/GMR),
                    sprintf("  %s %+5.1f%%", "n =",   100*(min(Ns)-n.est)/n.est)),
           bty="n", cex=0.80)
  } else {
    # ratios
    legend("topleft", 
           legend=c(paste0(design, " design", "; assumed:"),
                    sprintf("  %s %5.4f%s %5.4f", "CV =", CV, ", theta0 =", GMR),
                    "power:",
                    sprintf("  %s %5.4f", "target =", targetpower),
                    sprintf("  %s %5.4f %s %i%s", "estimated =", pwr.est,
                            "(n =", n.est, ")"),
                    sprintf("  %s %5.4f", "minimum acceptable =", minpower),
                    "acceptable rel. deviations:",
                    CVtxt,
                    sprintf("  %s %+5.2f%%", "GMR =", 100*(GMR.min-GMR)/GMR),
                    sprintf("  %s %+5.1f%%", "n =",   100*(min(Ns)-n.est)/n.est)),
           bty="n", cex=0.85)
  }
  
  close.screen(all.screens=TRUE)
  par(op) #reset options
}