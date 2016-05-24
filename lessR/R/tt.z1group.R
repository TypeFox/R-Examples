.OneGroup  <-
function(Y, Ynm, mu0=NULL, n=NULL, m=NULL, s=NULL, brief, bw1,
         from.data, conf.level, alternative, digits.d, mmd, msmd,
         Edesired, paired, graph, line.chart, show.title,
         pdf.file, pdf.width, pdf.height) { 

  # get variable label if exists
  gl <- .getlabels()
  y.lbl <- gl$yl

  if ( (!is.null(y.lbl)) ) {
    cat("Response Variable:  ", Ynm, ", ", as.character(y.lbl), sep="", "\n")
    cat("\n")
  }
  else y.lbl <- Ynm

  if (!brief) cat("\n------ Description ------\n\n")

  if (from.data) {
    n <- sum(!is.na(Y))
    n.miss <- sum(is.na(Y))
    Y <- na.omit(Y)
    m <- mean(Y)
    s <- sd(Y)
    v <- var(Y)
  }
  else {
    v <- s^2
  }

  clpct <- paste(toString(round((conf.level)*100, 2)), "%", sep="")

  if (from.data)
    dig.smr.d  <- digits.d
  else
    dig.smr.d  <- digits.d - 1
    

  if (Ynm != "Y") cat(Ynm,  ": ", sep="")
  if (from.data) cat(" n.miss = ", n.miss, ",  ", sep="") 
  cat("n = ", n, ",   mean = ", .fmt(m, dig.smr.d), 
      ",  sd = ", .fmt(s, dig.smr.d), sep="", "\n")

  if (brief) cat("\n")
  else  {
    if (from.data) {

      cat("\n\n------ Normality Assumption ------\n\n")
      # Normality
      if (n > 30) {
        cat("Sample mean assumed normal because n>30, so no test needed.", sep="", "\n")
      }
      else {
        cat("Null hypothesis is a normal distribution", sep="")
        if (Ynm != "Y") cat(" of ", Ynm, ".", sep="") else cat(".")
        cat("\n")
        if (n > 2 && n < 5000) {
          nrm1 <- shapiro.test(Y)
          W.1 <- round(nrm1$statistic,min(4,digits.d+1))
          p.val1 <- round(nrm1$p.value,min(4,digits.d+1))
          cat(nrm1$method, ":  W = ", W.1, ",  p-value = ", p.val1, sep="", "\n")
        }
        else
          cat("Sample size out of range for Shapiro-Wilk normality test.", "\n")
      }  
    } 
  }

  if (!brief) cat("\n\n------ Inference ------\n\n")

  # t-test
  if (!is.null(mu0)) m.dist <- m - mu0
  df <- n - 1
  sterr <- s * sqrt(1/n)
  if (alternative == "two.sided")
    tcut <- qt((1-conf.level)/2, df=df, lower.tail=FALSE)
  else if (alternative == "less")
    tcut <- qt(1-conf.level, df=df, lower.tail=FALSE)
  else if (alternative == "greater")
    tcut <- qt(1-conf.level, df=df, lower.tail=TRUE)
  if (from.data) {
    if (!is.null(mu0)) mu.null <- mu0 else mu.null <- 0
    ttest <- t.test(Y, conf.level=conf.level, alternative=alternative, mu=mu.null)
    df <- ttest$parameter
    lb <- ttest$conf[1]
    ub <- ttest$conf[2]
    E <- (ub-lb)/2
    if (!is.null(mu0)) {
      tvalue <- ttest$statistic
      pvalue <- ttest$p.value
    }
  }
  else {
    E <- tcut * sterr
    lb <- m - E
    ub <- m + E
    if (!is.null(mu0)) {
      tvalue <- m.dist/sterr
      if (alternative == "two.sided")
        pvalue <- 2 * pt(abs(tvalue), df=df, lower.tail=FALSE)
      else if (alternative == "less")
        pvalue <- pt(abs(tvalue), df=df, lower.tail=FALSE)
      else if (alternative == "greater")
        pvalue <- pt(abs(tvalue), df=df, lower.tail=TRUE)
    }
  }     

  cat("t-cutoff: tcut = ", .fmt(tcut,3), "\n") 
  cat("Standard Error of Mean: SE = ", .fmt(sterr), "\n\n")

  if (!is.null(mu0)) {
    cat("Hypothesized Value H0: mu =", mu0, "\n")
    txt <- "Hypothesis Test of Mean:  t-value = "
    cat(txt, .fmt(tvalue,3), ",  df = ", df, ",  p-value = ", .fmt(pvalue,3), sep="", "\n\n")
  }
  cat("Margin of Error for ", clpct, " Confidence Level:  ", .fmt(E), sep="", "\n")
  txt <- " Confidence Interval for Mean:  "
  cat(clpct, txt, .fmt(lb), " to ", .fmt(ub), sep="", "\n")

  # difference from mu0 and standardized mean difference
  if (!is.null(mu0)) {
    mdiff <- m - mu0
    smd <- abs(mdiff/s)
    if (!brief) cat("\n\n------ Effect Size ------\n\n") else cat("\n")
    cat("Distance of sample mean from hypothesized:  " , .fmt(mdiff), "\n",
        "Standardized Distance, Cohen's d:  ", .fmt(smd),
        sep="", "\n")
  }

  # needed sample size from Edesired
  if (!is.null(Edesired)) {
    zcut <- qnorm((1-conf.level)/2)
    ns <- ((zcut*s)/Edesired)^2 
    n.needed <- ceiling(1.132*ns + 7.368) 
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
    cat("Needed sample size: ", n.needed, "\n")
    cat("\n")
    cat("Additional data values needed: ", n.needed-n, "\n")
  }


  # graphs
  if (graph) {

    # keep track of the number of plots in this routine, see if manage graphics
    plt.i <- 0
    plt.title  <- character(length=0)
    manage.gr <- .graphman()
 
    if (is.null(pdf.file)) {
      if (manage.gr) {
        n.win <- 0
        if (!is.null(mu0)) n.win <- n.win + 1
        if (paired) n.win <- n.win + 1
        if (line.chart) n.win <- n.win + 1
        if (n.win > 0) {
          .graphwin(n.win)
          i.win <- 2   # start first graphics window on 3
          orig.params <- par(no.readonly=TRUE)
          on.exit(par(orig.params))
        }
      }
    }

    if (line.chart) {
      if (!is.null(pdf.file))
        pdf(file=paste("LineChart_",Ynm,".pdf",sep=""), width=pdf.width, height=pdf.height)

      if (manage.gr) {
        i.win  <- i.win + 1 
        dev.set(which=i.win)
      }

      plt.i <- plt.i + 1
      plt.title[plt.i] <- paste("Ordered Data for", Ynm)

      .lc.main(Y, type=NULL,
        col.line=getOption("col.stroke.pt"), col.area=NULL, col.box="black",
        col.stroke=getOption("col.stroke.pt"), 
        col.fill=getOption("col.fill.bar"), shape.pts=21,
        col.grid=getOption("col.grid"), col.bg=getOption("col.bg"),
        cex.axis=0.75, col.axis="gray30", rotate.values=0, offset=.5,
        xy.ticks=TRUE, line.width=1.1,
        xlab=NULL, ylab=NULL, main=plt.title[plt.i], sub=NULL, cex=NULL,
        time.start=NULL, time.by=NULL, time.reverse=FALSE,
        center.line="default", quiet=TRUE)

      if (!is.null(pdf.file)) {
        dev.off()
        .showfile(paste("LineChart_", Ynm, ".pdf", sep=""), paste("line chart of", Ynm))
      }

    }


  if (!is.null(mu0)) {

    if (manage.gr) {
      i.win  <- i.win + 1 
      dev.set(which=i.win)
    }

    if (!is.null(pdf.file))
      pdf(file=pdf.file, width=pdf.width, height=pdf.height)

    plt.i <- plt.i + 1
    plt.title[plt.i] <- "One-Group Plot"

    .OneGraph(Y, bw1, Ynm, y.lbl, digits.d, brief,
         n, m, mu0, mdiff, s, smd, mmd, msmd,
         clpct, tvalue, pvalue, ub, lb,
         show.title, pdf.file, pdf.width, pdf.height)

      if (!is.null(pdf.file)) {
        dev.off()
        .showfile(pdf.file, paste("density plot of", Ynm))
      }
    }

    return(list(i=plt.i, ttl=plt.title))
  }  # end if graph

}  # End One Group

