###
### R routines for the R package mvmeta (c) Antonio Gasparrini 2012-2014
#
print.summary.mvmeta <-  
function(x, digits=4, ...) {
#
################################################################################
#
  # CREATE USEFUL OBJECTS
  methodname <- c("reml","ml","fixed","mm","vc")
  methodlabel <- c("REML","ML","Fixed","Method of moments",
    "Variance components")
  bscovname <- c("unstr","diag","id","cs","hcs","ar1","prop","cor","fixed")
  bscovlabel <- c("General positive-definite","Diagonal",
    "Multiple of identity","Compound symmetry","Heterogeneous compound symmetry",
    "Autoregressive of first order","Proportional to fixed matrix",
    "Fixed correlation","Fixed")
  int <- attr(x$terms,"intercept")==1L
#
################################################################################
# HEADING AND SUBHEADING
#
  # HEADING
  cat("Call:  ",paste(deparse(x$call),sep="\n",collapse="\n"),"\n\n",sep="")
  #
  # SUB-HEADING
  cat(if(x$dim$k==1L)"Uni" else "Multi","variate ",
      ifelse(x$method=="fixed","fixed","random"),"-effects meta-",
      ifelse(x$dim$p-int>0,"regression","analysis"),"\n",sep="")
  # CHECK LATER FOR CHOICE META-ANALYSIS OR METAREGRESSION
  cat("Dimension: ",x$dim$k,"\n",sep="")
  if(x$method!="fixed") {
    cat("Estimation method: ",methodlabel[which(x$method==methodname)],
      "\n",sep="")
  }
  cat("\n")
#
###########################################################################
# FIXED EFFECTS ESTIMATES
#
  cat("Fixed-effects coefficients","\n",sep="")
#
  # COMPUTE SIGNIFICANCE
  signif <- symnum(x$coefficients[,"Pr(>|z|)"],corr=FALSE,na=FALSE,cutpoints=c(0, 
    0.001,0.01,0.05,0.1,1),symbols=c("***","**","*","."," "))
#
  # PRODUCE THE TABLE
  tabletot <- formatC(x$coefficients,digits=digits,format="f")
  tabletot <- cbind(tabletot,signif)
  colnames(tabletot)[7] <- ""
#
  # FOR SIMPLE META-ANALYSIS
  if(x$dim$p-int==0L) {
    if(x$dim$k>1L) rownames(tabletot) <- x$lab$k
    print(tabletot,quote=FALSE,right=TRUE,print.gap=2)
  }  
#
  # FOR META-REGRESSION
  if(x$dim$p-int>0L) {
    p <- x$dim$p
    for(i in seq(x$dim$k)) {
      ind <- seq((i-1)*p+1,(i-1)*p+p)
      table <- tabletot[ind,,drop=FALSE]
      rownames(table) <- x$lab$p
      if(x$dim$k>1) cat(x$lab$k[i],":","\n")
      print(table,quote=FALSE,right=TRUE,print.gap=2)
    }
  }
  cat("---\nSignif. codes: ",attr(signif,"legend"),"\n\n")
#
###########################################################################
# RANDOM COMPONENTS
#
  if(!x$method=="fixed") {  
    cat("Between-study random-effects (co)variance components","\n",sep="")
    if(x$dim$k>1L) cat("\t","Structure: ",bscovlabel[which(x$bscov==bscovname)],
      "\n",sep="")
#
    # STANDARD DEVIATIONS
    sd <- formatC(cbind("Std. Dev"=sqrt(diag(x$Psi))),digits=digits,format="f")
    if(x$dim$k==1L) rownames(sd) <- ""
#
    # CORRELATIONS
    if(x$dim$k>1L) {
      corRan <- x$corRandom
      corRan[upper.tri(x$corRan,diag=TRUE)] <- NA
      dimnames(corRan) <- NULL
      corRan <- format(corRan[-1,-ncol(corRan),drop=FALSE],digits=digits,
        format="f")
      corRan <- rbind(x$lab$k[-x$dim$k],corRan)
      colnames(corRan) <- c("Corr",rep("",x$dim$k-2))
      corRan[grep("NA",corRan)] <- ""
    } else corRan <- NULL
#
    # PRODUCE THE TABLE
    print(cbind(sd,corRan),quote=FALSE,right=TRUE,na.print="",print.gap=2)
#
    if(!is.null(x$negeigen) && x$negeigen>0) {
      cat("(Note: Truncated estimate - ",x$negeigen,
        " negative eigenvalues set to 0)","\n",sep="")
    }
    cat("\n")
  }
#
###########################################################################
# OVERALL QTEST AND I-SQUARE
#
  Q <- formatC(x$qstat$Q,digits=digits,format="f")
  pvalue <- formatC(x$qstat$pvalue,digits=digits,format="f")
  i2 <- formatC(pmax((x$qstat$Q-x$qstat$df)/x$qstat$Q*100,1),digits=1,format="f")
  cat(if(x$qstat$k==1) "Uni" else "Multi","variate ","Cochran Q-test for ",
    if(x$qstat$residual) "residual ", "heterogeneity:","\n",sep="")
  cat("Q = ",Q[1]," (df = ",x$qstat$df[1],"), p-value = ",pvalue[1],"\n",sep="")
  cat("I-square statistic = ",i2[1],"%","\n\n",sep="")
  
###########################################################################
# FIT STATS
#
  cat(x$dim$m," studies, ",x$df$nall," observations, ",x$df$fixed," fixed and ",
    x$df$random," random-effects parameters","\n",sep="")
  if(na <- length(x$na.action)) cat("(",na," stud",ifelse(na>1L,"ies","y"),
    " removed due to missingness",")\n",sep="")
  if(!x$method%in%c("mm","vc")) {
    table <- c(x$logLik,x$AIC,x$BIC)
    names(table) <- c("logLik","AIC","BIC")
    table <- formatC(table,digits=digits,format="f")
    print(table,quote=FALSE,right=TRUE,print.gap=2)
  }
  cat("\n")
#
}

