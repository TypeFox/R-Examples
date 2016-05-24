ttestPower <- 
function(n=NULL, s=NULL, n1=NULL, n2=NULL, s1=NULL, s2=NULL, 
         mmd=NULL, msmd=NULL, mdp=.8, mu0=NULL, 
         pdf.file=NULL, pdf.width=5, pdf.height=5, ...) {
      
  cat("\n")
  
  # for all null arguments, pick up values from previous smd.t
  if (sum(sapply(list(s, n1, n2), is.null)) == 3) {  # all are NULL
    if (exists("n1", 1, inherits=FALSE) && exists("n2", 1, inherits=FALSE)) {
      n1 <- get("n1", 1, inherits=FALSE)
      n2 <- get("n2", 1, inherits=FALSE)
    } 
    else {
      if (is.null(n)) {
        cat("\n"); stop(call.=FALSE, "\n","------\n",
        "Need to specify sample size, either n or n1 and n2.\n\n")
      }
    }
    if (is.null(s1) && is.null(s2)) { 
          cat("\n"); stop(call.=FALSE, "\n","------\n",
          "Need to specify a sample standard deviation, s,\n",
          "or two standard deviations, s1 and s2.\n\n")
    }
  }
  
  if ( (is.null(n)) && (is.null(n1)) && (is.null(n2)) ) { 
        cat("\n"); stop(call.=FALSE, "\n","------\n",
        "Need to specify either a common sample size for both groups, n,\n",
        "or specify both group sample sizes, n1 and n2, from which\n",
        "the within-group or pooled standard deviation, sw, is computed.\n\n")
  }

  if ( (is.null(s)) && (is.null(s1)) && (is.null(s2)) ) {
        cat("\n"); stop(call.=FALSE, "\n","------\n",
        "Need to specify either a single standard deviation, s, \n",
        "or specify two group standard deviations, s1 and s2, plus a\n",
        "common sample size n or individual sample sizes n1 and n2, from which\n",
        "the within-group standard deviation is computed.\n\n")
  }

  if ( (mdp < 0)  || (mdp > 1) ) {
        cat("\n"); stop(call.=FALSE, "\n","------\n",
        "Minimum desired power, mdp, must be between 0 and 1.\n\n")
  }

  if ( !is.null(mmd) && !is.null(msmd) ) {
        cat("\n"); stop(call.=FALSE, "\n","------\n",
        "Specify only one of mmd and msmd as one implies the other.\n\n")
  }

  if ( (!is.null(n1) || !is.null(n2)) && !is.null(mu0) ) {
        cat("\n"); stop(call.=FALSE, "\n","------\n",
        "Indicated two samples with n1 and n2 but only a single sample with mu0.\n\n")
  }
  
  
  cat("------------------------------------------------------------\n")
  if (is.null(mu0)) {
    mytype <- "two.sample"
    cat("Power Curve Analysis for Independent Groups t-test\n")
  }
  else {
    mytype <- "one.sample"
    cat("Power Curve Analysis for One Sample t-test\n")
    cat("------------------------------------------------------------\n")
    cat("mu0 =", mu0, "\n")
  }
  
  # power curve for two groups, assuming mean difference of 0
  if (mytype == "two.sample") {
    cat("------------------------------------------------------------\n")
    
    # n1, n2 and n all should be set
    if (!is.null(n)) { 
      cat("\n", "Sample size: n = ", n, sep="")
      n1 <- n;  n2 <- n 
    }
    else {
      if ( !is.null(n1) && !is.null(n2) ) {
        n = 2 / ( 1/n1 + 1/n2 )  # harmonic mean
        cat("\n", "Sample sizes: n1 = ", n1, " and n2 = ", n2, sep="")
        cat("\n", "Harmonic mean of the two sample sizes: n = ", n, sep="", "\n")
      }
    }
        
    # within-group standard deviation	(need n1 and n2)
    if (is.null(s)) {
      df1 <- n1 - 1
      df2 <- n2 - 1
      ssq <- (df1*s1^2 + df2*s2^2) / (df1 + df2)
      s <- sqrt(ssq)
      cat("\n", "Equal Group Variances Assumed", sep="", "\n")
    }
    cat("\n", "Within-group (pooled) Standard Deviation:  sw = ", s, sep="", "\n")
    
    mytitle <- "Power Curve for Independent Groups t-test"
    myxlab <- bquote(paste("Alternative Values of ", mu[1] - mu[2]))
    H0 <- 0
  }
  
  # power values for a single sample, triggered by nonzero mu0
  else {
    mytitle <- "Power Curve for One Sample t-test"
    myxlab <- bquote(paste("Alternative Values of ", mu))
    H0 <- mu0
  }
  cat("------------------------------------------------------------\n")
  
  # get value of mmd if msmd supplied
  if ( !is.null(mmd) | !is.null(msmd) ) {
    if (!is.null(mmd)) msmd <- mmd / s   # null, because msmd only informs mmd
    if (!is.null(msmd)) mmd <- msmd * s
  }
  

  # set up graphics system
  .opendev(pdf.file, pdf.width, pdf.height)
  
  .ttp2graph(myxlab, mytitle, n, s, mdp, mmd, msmd, mytype, H0, ...)

  if (!is.null(pdf.file)) .showfile(pdf.file, "power curve")

  cat("\n")
    
}
