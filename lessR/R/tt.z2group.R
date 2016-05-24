.TwoGroup <-
function(YA, YB, n1, n2, m1, m2, s1, s2, from.data,
         Ynm, Xnm, X1nm, X2nm, brief, digits.d,
         conf.level, alternative, mmd, msmd, Edesired, bw1, bw2, graph,
         line.chart, show.title, pdf.file, pdf.width, pdf.height, ...)  {        
 
  if ( brief  &&  (!is.null(mmd) || !is.null(msmd)) ) { 
    cat("\n"); stop(call.=FALSE, "\n","------\n",
      "mmd and msmd do not work with the brief version.\n\n")
  }

  cat("Compare", Ynm, "across", Xnm, "levels", X1nm, "and", X2nm, "\n")
  cat("------------------------------------------------------------\n\n")

  # get variable labels if exist
  options(xname = Xnm)
  options(yname = Ynm)
  gl <- .getlabels()
  x.lbl <- gl$xl
  y.lbl <- gl$yl

  if ( (!is.null(x.lbl)) || (!is.null(y.lbl)) ) {
    cat("Response Variable:  ", Ynm, ", ", as.character(y.lbl), sep="", "\n")
    cat("Grouping Variable:  ", Xnm, ", ", as.character(x.lbl), sep="", "\n")
    cat("\n")
  }
  if (is.null(y.lbl)) y.lbl <- Ynm

  if (!brief) cat("\n------ Description ------\n\n")

  if (from.data) {
    n1 <- sum(!is.na(YA))
    n2 <- sum(!is.na(YB))
    n1.miss <- sum(is.na(YA))
    n2.miss <- sum(is.na(YB))
    YA <- na.omit(YA)
    YB <- na.omit(YB)
    m1 <- mean(YA)
    m2 <- mean(YB)
    s1 <- sd(YA)
    s2 <- sd(YB)
    v1 <- var(YA)
    v2 <- var(YB)
  }
  else {
    v1 <- s1^2
    v2 <- s2^2
  }

  if (from.data) dig.smr.d  <- digits.d  else dig.smr.d <- digits.d - 1

  clpct <- paste(toString(round((conf.level)*100, 2)), "%", sep="")
  Xnmval <- paste(Xnm, X1nm)
  cat(Ynm, " for ", Xnmval, ":  ", sep="")
  if (from.data) cat("n.miss = ", n1.miss, ",  ", sep="")
  cat("n = ", n1, sep="")
  cat(",  mean = ", .fmt(m1,dig.smr.d), ",  sd = ", .fmt(s1,dig.smr.d), sep="", "\n")
  Xnmval <- paste(Xnm, X2nm)
  cat(Ynm, " for ", Xnmval, ":  ", sep="")
  if (from.data) cat("n.miss = ", n2.miss, ",  ", sep="")
  cat("n = ", n2, sep="")
  cat(",  mean = ", .fmt(m2,dig.smr.d), ",  sd = ", .fmt(s2,dig.smr.d), sep="", "\n")
  if (!brief) cat("\n")

# sw
  df1 <- n1 - 1
  df2 <- n2 - 1
  swsq <- (df1*v1 + df2*v2) / (df1 + df2)
  sw <- sqrt(swsq)
  if (!brief) cat("Within-group Standard Deviation:  ", .fmt(sw), "\n")

  # sample mean difference
  mdiff <- m1 - m2


  if (!brief) {

    cat("\n\n------ Assumptions ------\n\n")

    cat("Note: These hypothesis tests can perform poorly, and the", "\n")
    cat("      t-test is typically robust to violations of assumptions.", "\n")
    cat("      Use as heuristic guides instead of interpreting literally.", "\n\n")


    if (from.data) {

  # Normality
      cat("Null hypothesis, for each group, is a normal distribution of ", sep="")
      cat(Ynm, ".", sep="", "\n")
      if (n1 > 30) {
        cat("Group " , X1nm, ": ", sep="")
        cat("Sample mean assumed normal because n>30, so no test needed.", sep="", "\n")
      }
      else {
        cat("Group", X1nm, " ")
        if (n1 > 2 && n1 < 5000) {
          nrm1 <- shapiro.test(YA)
          W.1 <- nrm1$statistic
          p.val1 <- nrm1$p.value
          cat(nrm1$method, ":  W = ", .fmt(W.1,3), ",  p-value = ",
              .fmt(p.val1,3), sep="", "\n")
        }
      else
      cat("Sample size out of range for Shapiro-Wilk normality test.", "\n")
    }  
    if (n2 > 30) {
      cat("Group " , X2nm, ": ", sep="")
      cat("Sample mean assumed normal because n>30, so no test needed.", sep="", "\n")
    }
    else {
      cat("Group", X2nm, " ")
      if (n2 > 2 && n2 < 5000) {
        nrm2 <- shapiro.test(YB)
        W.2 <- nrm2$statistic
        p.val2 <- nrm2$p.value
        cat(nrm2$method, ":  W = ", .fmt(W.2,3), ",  p-value = ",
            .fmt(p.val2,3), sep="", "\n")
      }
      else
      cat("Sample size out of range for Shapiro-Wilk normality test.", "\n")
    }  
    cat("\n")
  } 

  # Homogeneity of Variance
  # Var Ratio
  if (v1 >= v2) {
    vratio <- v1/v2
      vr <- paste(.fmt(v1), "/", .fmt(v2), sep="")
      df.num <- df1
      df.den <- df2
  }
  else {
    vratio <- v2/v1
      vr <- paste(.fmt(v2), "/", .fmt(v1), sep="")
      df.num <- df2
      df.den <- df1
  }

  p.var <- pf(vratio, df1=df.num, df2=df.den)
  # adjust for two-sided test, results same as var.test{stats}
  p.var <- 2 * min(p.var, 1-p.var)

  cat("Null hypothesis is equal variances of ")
  cat(Ynm, ", i.e., homogeneous.", sep="", "\n")

  cat("Variance Ratio test:  F = ", vr, " = ", .fmt(vratio),
      ",  df = ", df.num, ";", 
      df.den, ",  p-value = ",  .fmt(p.var,3), sep="", "\n")

  if (from.data) { # Levene
    YAm <- abs(YA - median(YA))
      YBm <- abs(YB - median(YB))
      t.bf <- t.test(YAm, YBm, var.equal=TRUE)
      tvalue.bf <- t.bf$statistic
      df.bf <- t.bf$parameter
      pvalue.bf <- t.bf$p.value
      cat("Levene's test, Brown-Forsythe:  t = ", .fmt(tvalue.bf,3),
          ",  df = ", df.bf, sep="")
      cat(",  p-value = ", .fmt(pvalue.bf,3), sep="", "\n")
  }
}


  if (!brief) cat("\n\n------ Inference ------\n\n") else cat("\n---\n")

  # t-test
  sterr <- sw * sqrt(1/n1 + 1/n2)
  df <- df1 + df2
  if (alternative == "two.sided")
    tcut <- qt((1-conf.level)/2, df=df, lower.tail=FALSE)
  else if (alternative == "less")
    tcut <- qt(1-conf.level, df=df, lower.tail=FALSE)
  else if (alternative == "greater")
    tcut <- qt(1-conf.level, df=df, lower.tail=TRUE)
  if (from.data) {
    ttest <- t.test(YA, YB, var.equal=TRUE, conf.level=conf.level,
                    alternative=alternative)
    lb <- ttest$conf[1]
    ub <- ttest$conf[2]
    E <- (ub-lb) / 2
    tvalue <- ttest$statistic
    pvalue <- ttest$p.value
  }
  else {
    sterr <- sw * sqrt(1/n1 + 1/n2)
    E <- tcut*sterr
    lb <- mdiff-E
    ub <- mdiff+E
    tvalue <- mdiff/sterr
    if (alternative == "two.sided")
      pvalue <- 2 * pt(abs(tvalue), df=df, lower.tail=FALSE)
    else if (alternative == "less")
      pvalue <- pt(abs(tvalue), df=df, lower.tail=FALSE)
    else if (alternative == "greater")
      pvalue <- pt(abs(tvalue), df=df, lower.tail=TRUE)
  }

  if (!brief)
    cat("--- Assume equal population variances of", Ynm, "for each", Xnm, "\n\n")
  cat("t-cutoff: tcut = ", .fmt(tcut,3), "\n") 
  cat("Standard Error of Mean Difference: SE = ", .fmt(sterr), "\n")
 
  mytitle <- "\nHypothesis Test of 0 Mean Diff:  t = "
  cat(mytitle, .fmt(tvalue,3), ",  df = ", df, ",  p-value = ", .fmt(pvalue,3),
      sep="", "\n\n")
  cat("Margin of Error for ", clpct, " Confidence Level:  ", .fmt(E), sep="", "\n")
  cat(clpct," Confidence Interval for Mean Difference:  ", .fmt(lb), " to ", .fmt(ub), 
      sep="", "\n\n")

  if (!brief) {
    k1 <- v1/n1
    k2 <- v2/n2
    df.ne <- ((k1 + k2)^2) / ((k1^2)/(n1-1) + (k2^2)/(n2-1))
    sterr.ne <- sqrt(k1 + k2)
    if (alternative == "two.sided")
      tcut.ne <- qt((1-conf.level)/2, df=df.ne, lower.tail=FALSE)
    else if (alternative == "less")
      tcut.ne <- qt(1-conf.level, df=df.ne, lower.tail=FALSE)
    else if (alternative == "greater")
      tcut.ne <- qt(1-conf.level, df=df.ne, lower.tail=TRUE)

    if (from.data) {
      ttne <- t.test(YA, YB, var.equal=FALSE, conf.level=conf.level,
                     alternative=alternative)
      df.ne <- ttne$parameter
      lb.ne <- ttne$conf[1]
      ub.ne <- ttne$conf[2]
      E.ne <- (ub.ne-lb.ne)/2
      tvalue.ne <- ttne$statistic
      pvalue.ne <- ttne$p.value
    }
    else {
      E.ne <- tcut.ne*sterr.ne
      lb.ne <- mdiff-E.ne
      ub.ne <- mdiff+E.ne
      tvalue.ne <- mdiff / sterr.ne
    if (alternative == "two.sided")
      pvalue.ne <- 2 * pt(abs(tvalue.ne), df=df.ne, lower.tail=FALSE)
    else if (alternative == "less")
      pvalue.ne <- pt(abs(tvalue.ne), df=df.ne, lower.tail=FALSE)
    else if (alternative == "greater")
      pvalue.ne <- pt(abs(tvalue.ne), df=df.ne, lower.tail=TRUE)
    }

    cat("\n--- Do not assume equal population variances of", Ynm, "for each", Xnm, "\n\n")
    cat("t-cutoff: tcut = ", .fmt(tcut.ne,3), "\n") 
    cat("Standard Error of Mean Difference: SE = ", .fmt(sterr.ne), "\n")
    mytitle <- "\nHypothesis Test of 0 Mean Diff:  t = "
    cat(mytitle, .fmt(tvalue.ne,3), ",  df = ", .fmt(df.ne,3), 
        ", p-value = ", .fmt(pvalue.ne,3), sep="", "\n\n")
    cat("Margin of Error for ", clpct, " Confidence Level:  ", .fmt(E.ne), sep="", "\n")
    cat(clpct," Confidence Interval for Mean Difference:  ", .fmt(lb.ne), " to ",
        .fmt(ub.ne), sep="", "\n")

  }

  # mean difference and standardized mean difference
  smd <- mdiff/sw
  if (!brief) {
    cat("\n\n------ Effect Size ------\n\n")
    cat("--- Assume equal population variances of", Ynm, "for each", Xnm, "\n\n")
  }
  cat("Sample Mean Difference of ", Ynm, ":  " , .fmt(mdiff), sep="", "\n")
  cat("Standardized Mean Difference of ", Ynm, ", ",
      "Cohen's d:  ", .fmt(smd), sep="", "\n")

  #cid <- ci.smd(smd=smd, n.1=n1, n.2=n2, conf.level=conf.level)  # MBESS function
  #deltaL <-cid$Lower.Conf.Limit.smd
  #deltaU <- cid$Upper.Conf.Limit.smd
  #if (!brief) cat("\n")
  #cat(clpct," Confidence Interval for smd:  ",
      #.fmt(deltaL), " to ", .fmt(deltaU), sep="", "\n")

  if (!brief) {
    cat("\n\n------ Practical Importance ------\n\n")
    cat("Minimum Mean Difference of practical importance: mmd\n")
    if ( !is.null(mmd) | !is.null(msmd) ) {
      if (!is.null(mmd)) msmd <- mmd / sw
      if (!is.null(msmd)) mmd <- msmd * sw
      cat("Compare mmd =", .fmt(mmd,digits.d), " to the obtained value of md = ",
        .fmt(mdiff), "\n")
      cat("Compare mmd to the confidence interval for md: ", .fmt(lb), " to ", 
          .fmt(ub), "\n\n")
      cat("Minimum Standardized Mean Difference of practical importance: msmd\n")
      cat("Compare msmd = ", .fmt(msmd,digits.d),
          " to the obtained value of smd = ", .fmt(smd),"\n")
      #if (!is.null(deltaL)) cat("Compare smd to the confidence interval for smd: ", 
          #.fmt(deltaL), " to ", .fmt(deltaU), "\n")
    }
    else {
      cat("Minimum Standardized Mean Difference of practical importance: msmd\n")
      cat("Neither value specified, so no analysis\n")
    }
  }

  # needed sample size from Edesired
  if (!is.null(Edesired)) {
    zcut <- qnorm((1-conf.level)/2)
    ns <- 2*((zcut*sw)/Edesired)^2 
    n.needed <- ceiling(1.099*ns + 4.863) 
    cat("\n\n------ Needed Sample Size ------\n\n")
    if (Edesired > E) {
      cat("Note: Desired margin of error,", .fmt(Edesired), 
         "is worse than what was obtained,", .fmt(E), "\n\n") 
    }
    cat("Desired Margin of Error: ", .fmt(Edesired), "\n")
    cat("\n")
    cat("For the following sample size there is a 0.9 probability of obtaining\n")
    cat("the desired margin of error for the resulting 95% confidence interval.\n")
    cat("-------\n")
    cat("Needed sample size per group: ", n.needed, "\n")
    cat("\n")
    cat("Additional data values needed Group 1: ", n.needed-n1, "\n")
    cat("Additional data values needed Group 2: ", n.needed-n2, "\n")
  }

  # graphs
  if (from.data && graph) {

    # keep track of the number of plots in this routine, see if manage graphics
    plt.i <- 0
    plt.title  <- character(length=0)
    manage.gr <- .graphman()

    if (is.null(pdf.file)) {
      if (manage.gr) {
        if (!line.chart) n.win <- 1  else n.win <- 3
        .graphwin(n.win)
        dev.set(which=3)
        orig.params <- par(no.readonly=TRUE)
        on.exit(par(orig.params))
      }
    }

    if (line.chart) { 

      if (!is.null(pdf.file))
        pdf(file=paste("LineChart_",X1nm,".pdf",sep=""), width=pdf.width, height=pdf.height)

      plt.i <- plt.i + 1
      plt.title[plt.i] <- paste("Ordered Data:", paste(Xnm, X1nm))

     .lc.main(YA, type=NULL,
       col.line=getOption("col.stroke.pt"), col.area=NULL, col.box="black",
       col.stroke=getOption("col.stroke.pt"), 
       col.fill=getOption("col.fill.bar"), shape.pts=21,
       col.grid=getOption("col.grid"), col.bg=getOption("col.bg"),
       cex.axis=0.75, col.axis="gray30", rotate.values=0, offset=.5,
       xy.ticks=TRUE, line.width=1.1,
       xlab=NULL, ylab=paste(Ynm,": ",X1nm, sep=""),
       main=plt.title[plt.i], sub=NULL,
       cex=NULL, time.start=NULL, time.by=NULL, time.reverse=FALSE,
       center.line="default", quiet=TRUE)

     if (is.null(pdf.file)) {
       if (manage.gr) dev.set(which=4)
     }
     else {
       dev.off()
       pdf(file=paste("LineChart_",X2nm,".pdf",sep=""), width=pdf.width, height=pdf.height)
       .showfile(paste("LineChart_", X2nm, ".pdf", sep=""),
            paste("line chart of", Ynm, "for Group", X2nm))
     }

      plt.i <- plt.i + 1
      plt.title[plt.i] <- paste("Ordered Data:", paste(Xnm, X2nm))
 
     .lc.main(YB, type=NULL,
       col.line=getOption("col.stroke.pt"), col.area=NULL, col.box="black",
       col.stroke=getOption("col.stroke.pt"), 
       col.fill=getOption("col.fill.bar"), shape.pts=21,
       col.grid=getOption("col.grid"), col.bg=getOption("col.bg"),
       cex.axis=0.85, col.axis="gray30", rotate.values=0, offset=.5,
       xy.ticks=TRUE, line.width=1.1,
       xlab=NULL, ylab=paste(Ynm,": ",X2nm, sep=""),
       main=plt.title[plt.i], sub=NULL,
       cex=NULL, time.start=NULL, time.by=NULL, time.reverse=FALSE,

       center.line="default", quiet=TRUE)

      if (is.null(pdf.file)) {
        if (manage.gr) dev.set(which=5)
      }
      else {
        dev.off()
        .showfile(paste("LineChart_", X2nm, ".pdf", sep=""),
            paste("line chart of", Ynm, "for Group", X2nm))
      }
    }

    if (!is.null(pdf.file))
      pdf(file=pdf.file, width=pdf.width, height=pdf.height)

      plt.i <- plt.i + 1
      plt.title[plt.i] <- "Two-Group Plot"

    .TwoGraph(YA, YB, bw1, bw2, Ynm, Xnm, X1nm, X2nm, y.lbl, digits.d, brief,
              n1, m1, s1, n2, m2, s2, df, mdiff, sw, smd, mmd, msmd,
              clpct, tvalue, pvalue, ub, lb, show.title)
              #clpct, tvalue, pvalue, ub, lb, deltaL, deltaU, show.title)

    if (!is.null(pdf.file)) {
      dev.off()
      .showfile(pdf.file, paste("density plots of", Ynm, "for both groups"))
    }

    return(list(i=plt.i, ttl=plt.title))
  }

} # End Two Group
