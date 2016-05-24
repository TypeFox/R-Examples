.ANOVAz2 <- 
function(av.out, y.values, x1.values, x2.values, nm, digits.d, brief,
         delim, rb.points, graphics, pdf, pdf.width, pdf.height) {

  if (grepl("*", delim, fixed=TRUE)) bet.grp  <- TRUE else bet.grp <- FALSE
  if (grepl("+", delim, fixed=TRUE)) wth.grp  <- TRUE else wth.grp <- FALSE

  p <- length(unique(na.omit(x1.values)))
  if (bet.grp)
    q <- length(unique(na.omit(x2.values)))
  if (wth.grp) {
    n <- length(unique(na.omit(x2.values)))
    q <- 1
  }


  tx <- character(length = 0)

  tx[length(tx)+1] <- "" 
  if (bet.grp)
    tx[length(tx)+1] <- "Two-way Between Groups ANOVA"
  else if (wth.grp) {
    tx[length(tx)+1] <- "Randomized Blocks ANOVA"
    tx[length(tx)+1] <- paste("  Factor of Interest: ", nm[2])
    tx[length(tx)+1] <- paste("  Blocking Factor:    ", nm[3])
    tx[length(tx)+1] <- ""
    tx[length(tx)+1] <- paste(
        "Note: For the resulting F statistic for", nm[2], "to be distributed as F,\n",
        "     the population covariances of", nm[1], "must be spherical.")
  }

  txbck2 <- tx


  title_des <- "  DESCRIPTIVE STATISTICS "

  txcn <- ""
  txcm <- ""
  if (bet.grp) {
    l <-  tapply(y.values, 
          list(x1.values, x2.values), length)
    l <- as.table(l)

    if (!brief) {

      tx <- character(length = 0)
      tx[length(tx)+1] <- paste("Cell Sample Size:", l[1,1])
      txcn <- tx


      tx <- character(length = 0)

      if (is.null(options()$knitr.in.progress)) {
        tx[length(tx)+1] <- "Cell Means"
        tx[length(tx)+1] <- ""
      }

      m <-  tapply(y.values, 
            list(x1.values, x2.values ), mean, na.rm=TRUE)
      m <- as.table(m)
      #names(dimnames(m)) <- c(nm[2], nm[3])
      tx2 <- .prntbl(t(m), digits.d, cc=NULL, v1.nm=nm[2], v2.nm=nm[3])
      for (i in 1:length(tx2)) tx[length(tx)+1] <- tx2[i]

      txcm <- tx

    }  # !brief
  }  # bet.grp


    tx <- character(length = 0)

    if (is.null(options()$knitr.in.progress)) {
      tx[length(tx)+1] <- "Marginal Means"
      tx[length(tx)+1] <- ""
    }

    tx[length(tx)+1] <- nm[2]
    m1 <-  tapply(y.values, x1.values, mean, na.rm=TRUE)
    m1 <- data.frame(t(m1))
    tx2 <- .prntbl(m1, digits.d)  # 1st treatment horizontal dimension
    for (i in 1:length(tx2)) tx[length(tx)+1] <- tx2[i]
    tx[length(tx)+1] <- ""
    tx[length(tx)+1] <- nm[3]
    m2 <-  tapply(y.values, x2.values, mean, na.rm=TRUE)
    m2 <- data.frame(t(m2))
    tx2 <- .prntbl(m2, digits.d)  # 2nd treatment horizontal dimension
    for (i in 1:length(tx2)) tx[length(tx)+1] <- tx2[i]

    txmm <- tx


    tx <- character(length = 0)

    if (is.null(options()$knitr.in.progress)) {
      tx[length(tx)+1] <- "Grand Mean"
      tx[length(tx)+1] <- ""
    }

    mg <- mean(y.values, na.rm=TRUE)
    tx[length(tx)+1] <- round(mg, digits.d+1)

    txgm <- tx


  tx <- character(length = 0)

  if (is.null(options()$knitr.in.progress)) {
    tx[length(tx)+1] <- "Cell Standard Deviations"
    tx[length(tx)+1] <- ""
  }

  txcs <- ""
  if (bet.grp) {

    s <-  tapply(y.values, 
          list(x1.values, x2.values ), sd, na.rm=TRUE)
    s <- as.table(s)
    tx2 <- .prntbl(t(s), digits.d, cc=NULL, v1.nm=nm[2], v2.nm=nm[3])
    for (i in 1:length(tx2)) tx[length(tx)+1] <- tx2[i]

    txcs <- tx
  }  # end between groups



  title_basic <- "  BASIC ANALYSIS"

  tx <- character(length = 0)

  if (is.null(options()$knitr.in.progress)) {
    tx[length(tx)+1] <- paste("Analysis of Variance")
    tx[length(tx)+1] <- ""
  }

  if(bet.grp) n.vars <- 4
  if(wth.grp) n.vars <- 3
  smc <- anova(av.out)
  buf <- 0 
  for (i in 1:n.vars) {
    lng.lbl <- nchar(rownames(smc)[i])
    if (lng.lbl > buf) buf <- lng.lbl 
   }
  max.num <- integer(length=0)
  for (icol in 1:3) {
    max.num[icol] <- 0 
    for (i in 1:n.vars) {
      ln.nm <- nchar(as.character(trunc(smc[i,icol]))) + digits.d + 2
      if (ln.nm > max.num[icol]) max.num[icol] <- ln.nm
    }
    if (icol != 1) if (max.num[icol] < 9L) max.num[icol] <- 9L 
  }
  df.lbl <- .fmtc("     df", max.num[1]+2)
  SS.lbl <- .fmtc(" Sum Sq", max.num[2]+1)
  MS.lbl <- .fmtc("Mean Sq", max.num[3]+1)
  if (bet.grp) 
    max.nm <- max(nchar(nm[2])+nchar(nm[3])+1, nchar("Residuals"))
  else
    max.nm <- max(nchar(nm[2]), nchar(nm[3]), nchar("Residuals"))
    
  tx[length(tx)+1] <- paste(format("", width=max.nm-5), df.lbl, SS.lbl, MS.lbl,
                                     "   F-value", "   p-value", sep="")
  for (i in 1:(n.vars)) {
    rlb <- .fmtc(rownames(smc)[i], buf)
    df <- format(sprintf("%i", smc[i,1]), width=max.num[1]-4, justify="right")
    SS <- format(sprintf("%7.*f", digits.d, smc[i,2]), width=max.num[2], justify="right")
    MS <- format(sprintf("%7.*f", digits.d, smc[i,3]), width=max.num[3], justify="right")
    if (i < n.vars) {
      fv <- format(sprintf("%7.*f", digits.d, smc[i,4]), width=9, justify="right")
      pv <- format(sprintf("%6.4f", smc[i,5]), width=9, justify="right")
      tx[length(tx)+1] <- paste(rlb, df, SS, MS, fv, pv) 
    }
    else
      tx[length(tx)+1] <- paste(rlb, df, SS, MS) 
  }

  txanv <- tx


  tx <- character(length = 0)

  if (is.null(options()$knitr.in.progress)) {
    tx[length(tx)+1] <- paste("Association and Effect Size")
    tx[length(tx)+1] <- ""
  }

  sm <- summary(av.out)
  msA <- sm[[1]][1,3]
  msB <- sm[[1]][2,3]
  msAB <- sm[[1]][3,3]
  msE <- sm[[1]][4,3]
  FA <- sm[[1]][1,4]
  FB <- sm[[1]][2,4]
  if (bet.grp) {
    FAB <- sm[[1]][3,4]
    n <- l[1,1]
  }

  omsq.A <- ((p-1)*(FA-1)) / ((p-1)*(FA-1) + n*p*q)
  if (bet.grp) omsq.B <- ((q-1)*(FB-1)) / ((q-1)*(FB-1) + n*p*q)
  if (wth.grp) intra.B <- (FB-1) / ((p-1)+(FB))
  if (bet.grp) omsq.AB <- ((p-1)*(q-1)*(FAB-1)) / ((p-1)*(q-1)*(FAB-1) + n*p*q)

  tx[length(tx)+1] <- paste("Partial Omega Squared for ", nm[2],
      ": ", .fmt(omsq.A, 2), sep="")
  if (bet.grp)
    tx[length(tx)+1] <- paste("Partial Omega Squared for ", nm[3],
        ": ", .fmt(omsq.B, 2), sep="")
  if (wth.grp)
    tx[length(tx)+1] <- paste("Partial Intraclass Correlation for ", nm[3],
        ": ", .fmt(intra.B, 2), sep="")
  if (bet.grp)
    tx[length(tx)+1] <- paste("Partial Omega Squared for ", nm[2], " & ", nm[3],
        ": ", .fmt(omsq.AB, 2), sep="")
  tx[length(tx)+1] <- ""

  if (omsq.A > 0) {
    fA.cohen <- sqrt( (omsq.A/(1-omsq.A)) )
    tx[length(tx)+1] <- paste("Cohen's f for ", nm[2], ": ", .fmt(fA.cohen, 2), sep="")
  }
  if (bet.grp) {
    if(omsq.B > 0) {
      fB.cohen <- sqrt( (omsq.B/(1-omsq.B)) )
      tx[length(tx)+1] <- paste("Cohen's f for ", nm[3], ": ", .fmt(fB.cohen, 2), sep="")
    }
    if (omsq.AB > 0) {
      fAB.cohen <- sqrt( (omsq.AB/(1-omsq.AB)) )
      tx[length(tx)+1] <- paste("Cohen's f for ", nm[2], "_&_", nm[3], ": ", .fmt(fAB.cohen, 2), sep="")
    }
  }
  if (wth.grp) {
    if (intra.B > 0) {
      fB.cohen <- sqrt( (intra.B/(1-intra.B)) )
      tx[length(tx)+1] <- paste("Cohen's f for ", nm[3], ": ", .fmt(fB.cohen, 2), sep="")
    }
  }

  txeft <- tx


  if (!brief) {
    tx <- character(length = 0)

    if (is.null(options()$knitr.in.progress))
      tx[length(tx)+1] <- paste("Tukey Multiple Comparisons of Means")

    HSD <- TukeyHSD(av.out)
    tx[length(tx)+1] <- ""
    tx[length(tx)+1] <- paste("Family-wise Confidence Level:", attr(HSD, which="conf.level"))
    tx[length(tx)+1] <- paste("\nFactor:", nm[2])

    txHSD <- .prntbl(HSD[[1]], digits.d)
    for (i in 1:length(txHSD)) tx[length(tx)+1] <- txHSD[i]

    if (bet.grp) {
      tx[length(tx)+1] <- paste("\nFactor:", nm[3])
      txHSD <- .prntbl(HSD[[2]], digits.d)  # second factor
      for (i in 1:length(txHSD)) tx[length(tx)+1] <- txHSD[i]
      tx[length(tx)+1] <- paste("\nCell Means")
      txHSD <- .prntbl(HSD[[3]], digits.d)  # interaction
      for (i in 1:length(txHSD)) tx[length(tx)+1] <- txHSD[i]
    }
  }

    txhsd <- tx


  # ------------------------------------
  # keep track of the number of plots, see if manage graphics
  plt.i <- 0
  plt.title  <- character(length=0)
  if (graphics) {
    manage.gr <- .graphman()

    # interaction plots
    if (!pdf) {
      if (manage.gr) {
        pdf.file <- NULL
        if (bet.grp) .graphwin(1)
        if (wth.grp) .graphwin(2)
        dev.set(which=3)
      }
    }
    else { 
      pdf.file <- "ANOVA_Interaction.pdf"
      pdf(file=pdf.file, width=pdf.width, height=pdf.height)
    }

    plt.i <- plt.i + 1
    plt.title[plt.i] <- "Interaction Plot"

    interaction.plot(x1.values, x2.values, y.values,
                     xlab=nm[2], ylab=nm[1], trace.label=nm[3],
                     main=plt.title[plt.i])

    # pdf
    if (pdf) {
      dev.off()
      .showfile(pdf.file, "interaction plot")
    }

    # fitted interaction plots
    if (wth.grp) {

      if (!pdf) {
        if (manage.gr) {
          dev.set(which=4)
        }
      }
      else { 
        pdf.file <- "ANOVA_FitInter.pdf"
        pdf(file=pdf.file, width=pdf.width, height=pdf.height)
      }

    plt.i <- plt.i + 1
    plt.title[plt.i] <- "Fitted Values"

      mn.y <- min(y.values, av.out$fitted)
      mx.y <- max(y.values, av.out$fitted)
      interaction.plot(x1.values, x2.values, 
               av.out$fitted, main=plt.title[plt.i],
               xlab=nm[2], ylab=nm[1], trace.label=nm[3], ylim=c(mn.y, mx.y))
      if (rb.points) {
        points(x1.values, y.values, pch=21, 
               bg=rgb(.6, .6, .6, alpha=getOption("trans.pts"), maxColorValue = 1))
        segments(as.numeric(x1.values), av.out$fitted, as.numeric(x1.values), y.values)
      }

      # pdf
      if (pdf) {
        dev.off()
        .showfile(pdf.file, "fitted values plot")
      }
    }
  }

  return(list(
    txbck2=txbck2, 
    title_des=title_des, txcn=txcn, txcm=txcm, txmm=txmm, txgm=txgm, txcs=txcs,
    title_basic=title_basic,
    txanv=txanv, txeft=txeft, txhsd=txhsd,
    i=plt.i, ttl=plt.title))

} 
