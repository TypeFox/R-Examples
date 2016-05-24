.ANOVAz1 <- 
function(av.out, y.values, x.values, nm, n.obs, digits.d, brief,
         graphics, pdf, pdf.width, pdf.height) {

  p <- length(unique(na.omit(x.values)))

  # cell stats
  n <- tapply(y.values, x.values, length) 
  m <- tapply(y.values, x.values, mean) 
  s <- tapply(y.values, x.values, sd) 
  mn <- tapply(y.values, x.values, min) 
  mx <- tapply(y.values, x.values, max) 

  if (!brief) {
    # get maximum chars in 1st three columns
    max.lv <- 0; max.n <- 0; max.m <- 0; max.s <- 0; max.mn <- 0; max.mx <- 0
    for (i in 1:p) {
      nch.lv <- nchar(as.character(levels(x.values)[i]))
      nch.n <- nchar(as.character(n[i]))
      nch.m <- nchar(.fmt(m[i],digits.d))
      nch.s <- nchar(.fmt(s[i],digits.d))
      nch.mn <- nchar(.fmt(mn[i],digits.d))
      nch.mx <- nchar(.fmt(mx[i],digits.d))
      if (nch.lv > max.lv) max.lv <- nch.lv
      if (nch.n > max.n) max.n <- nch.n
      if (nch.m > max.m) max.m <- nch.m
      if (nch.s > max.s) max.s <- nch.s
      if (nch.mn > max.mn) max.mn <- nch.mn
      if (nch.mx > max.mx) max.mx <- nch.mx
    }

    title_des <- "  DESCRIPTIVE STATISTICS "

    tx <- character(length = 0)

    n.tx <- format("n", width=max.lv+max.n+2, justify="right", sep="")
    m.tx <- format("mean", width=max.m+2, justify="right", sep="")
    s.tx <- format("sd", width=max.s+2, justify="right", sep="")
    mn.tx <- format("min", width=max.mn+2, justify="right", sep="")
    mx.tx <- format("max", width=max.mx+2, justify="right", sep="")
    tx[length(tx)+1] <- paste(n.tx, m.tx, s.tx, mn.tx, mx.tx)

    for (i in 1:p) {
      xval <- paste(format(levels(x.values)[i], width=max.lv, justify="left", sep="")) 
      nn <- format(n[i], width=max.n+1, justify="right", sep="") 
      mm <- format(sprintf("%.*f", digits.d, m[i]), width=max.m+2, justify="right") 
      ss <- format(sprintf("%.*f", digits.d, s[i]), width=max.s+2, justify="right") 
      mnn <- format(sprintf("%.*f", digits.d, mn[i]), width=max.mn+2, justify="right") 
      mxx <- format(sprintf("%.*f", digits.d, mx[i]), width=max.mx+2, justify="right") 
      tx[length(tx)+1] <- paste(xval, nn, mm, ss, mnn, mxx) 
    }

      mg <-  mean(y.values, na.rm=TRUE)
      tx[length(tx)+1] <- ""
      tx[length(tx)+1] <- paste("Grand Mean:", round(mg, digits.d+1))
  }  # end !brief 

  txdes <- tx


  # set up graphics system for 2 windows
  plt.i <- 0
  plt.title  <- character(length=0)
  if (graphics) {

  # keep track of the number of plots, see if manage graphics
    manage.gr <- .graphman()

    if (!pdf) {
      if (manage.gr) {
        .graphwin(2)
        dev.set(which=3)
      }
    }
    else { 
      pdf.file <- "ANOVA_Means.pdf"
      pdf(file=pdf.file, width=pdf.width, height=pdf.height)
    }

    plt.i <- plt.i + 1
    plt.title[plt.i] <- "Scatterplot with Cell Means"

    # scatter plot
    plot(x.values, y.values, type="n", axes=FALSE, ann=FALSE)

    usr <- par("usr")
    col.bg <- getOption("col.bg")
    rect(usr[1], usr[3], usr[2], usr[4], col=col.bg, border="black")

    axT1 <- 1:length(unique(x.values))

    col.grid <- getOption("col.grid")
    abline(v=axT1, col=col.grid, lwd=.5)
    abline(h=axTicks(2), col=col.grid, lwd=.5)

    .axes(levels(x.values), NULL, axT1, axTicks(2),
          par("usr")[1], par("usr")[3], cex.axis=.8, col.axis="gray30")

    main.lab <- plt.title[plt.i]
    sub.lab <- NULL
    x.label <- nm[2]
    y.label <- nm[1]
    .axlabs(x.label, y.label, main.lab, sub.lab, max.lbl.y=3, cex.lab=0.85) 

    col.fill <- getOption("col.fill.pt")
    col.stroke <- getOption("col.stroke.pt")
    points(x.values, y.values, pch=21, col=col.stroke, bg=col.fill, cex=0.8)

    # plot cell means
    pch.avg <- ifelse(getOption("colors")!="gray", 21, 23)
    bck.g <- ifelse(getOption("colors")!="gray", "gray15", "gray30")
    if (grepl(".black", getOption("colors"), fixed=TRUE)) bck.g <- "gray85"

    m.lvl <- numeric(length = 0)
    for (i in (1:length(levels(x.values))))
      m.lvl[i] <- mean(y.values[which(x.values==levels(x.values)[i])], na.rm=TRUE)
    abline(h=m.lvl, col="gray50", lwd=.5)
    points(m.lvl, pch=pch.avg, bg=bck.g)

    if (pdf) {
      dev.off()
      .showfile(pdf.file, "means chart")
    }
  }


  title_basic <- "  BASIC ANALYSIS"

  tx <- character(length = 0)

  if (is.null(options()$knitr.in.progress)) {
    tx[length(tx)+1] <- paste("Analysis of Variance")
    tx[length(tx)+1] <- ""
  }

  n.vars <- 2
  smc <- anova(av.out)

  # width of column 1
  max.c1 <- max(nchar("Residuals"), nchar(nm[1]))
  for (i in 1:n.vars) {
    c1 <- nchar(rownames(smc)[i])
    if (c1 > max.c1) max.c1 <- c1 
  }

  # width of data columns
  max.ln <- integer(length=0)
  for (i in 1:4) {
    ln.nm <- nchar(colnames(smc)[i])
    max.ln[i] <- ln.nm + 1
    for (j in 1:nrow(smc)) {
      xjc <- .fmt(smc[j,i], d=digits.d)
      if (nchar(xjc) > max.ln[i]) max.ln[i] <- nchar(xjc)
    }
    max.ln[i] <- max.ln[i] + 1L
    if (max.ln[i] < 9L) max.ln[i] <- 9L
  }

  df.lbl <- .fmtc("     df", max.ln[1]+1)
  SS.lbl <- .fmtc(" Sum Sq", max.ln[2]+1)
  MS.lbl <- .fmtc("Mean Sq", max.ln[3]+1)
  fv.lbl <- .fmtc("F-value", max.ln[4]+1)
  tx[length(tx)+1] <- paste(format("", width=max.c1-4), df.lbl, SS.lbl,
                             MS.lbl, fv.lbl, "   p-value", sep="")

  for (i in 1:n.vars) {
    rlb <- .fmtc(rownames(smc)[i], max.c1, j="left")
    df <- format(sprintf("%i", smc[i,1]), width=max.ln[1]-4, justify="right")
    SS <- format(sprintf("%7.*f", digits.d, smc[i,2]), width=max.ln[2], justify="right")
    MS <- format(sprintf("%7.*f", digits.d, smc[i,3]), width=max.ln[3], justify="right")
    if (i < n.vars) {
      fv <- format(sprintf("%7.*f", digits.d, smc[i,4]), width=9, justify="right")
      pv <- format(sprintf("%6.4f", smc[i,5]), width=9, justify="right")
      tx[length(tx)+1] <- paste(rlb, df, SS, MS, fv, pv)
    }
    else
      tx[length(tx)+1] <- paste(rlb, df, SS, MS)
    }

  sm <- summary(av.out)
  ssb <- sm[[1]][1,2]
  ssw <- sm[[1]][2,2]
  sst <- ssb + ssw
  msw <- sm[[1]][2,3]

  txanv <- tx


  tx <- character(length = 0)

  if (is.null(options()$knitr.in.progress)) {
    tx[length(tx)+1] <- paste("Association and Effect Size")
    tx[length(tx)+1] <- ""
  }

  rsq <- ssb / sst
  tx[length(tx)+1] <- paste("R Squared:", .fmt(rsq, 2))
  rsq.adj <-  1 - ( ((n.obs-1)/(n.obs-p)) * (1-rsq) )
  tx[length(tx)+1] <- paste("R Sq Adjusted:", .fmt(rsq.adj, 2))
  omsq <- (ssb - ((p-1)*msw)) / ((ssb + p*(mean(n)-1)*msw) + msw)
  tx[length(tx)+1] <- paste("Omega Squared:", .fmt(omsq, 2))
  if (omsq > 0) {
    tx[length(tx)+1] <- ""
    f.cohen <- sqrt( (omsq/(1-omsq)) )
    tx[length(tx)+1] <- paste("Cohen's f:", .fmt(f.cohen, 2))
  }

  txeft <- tx


  if (!brief) {
    tx <- character(length = 0)

    if (is.null(options()$knitr.in.progress)) {
      tx[length(tx)+1] <- paste("Tukey Multiple Comparisons of Means")
      tx[length(tx)+1] <- ""
    }

    HSD <- TukeyHSD(av.out)
    HSD <- TukeyHSD(av.out, which=nm[2])
    tx[length(tx)+1] <- paste("Family-wise Confidence Level:", attr(HSD, which="conf.level"))
    txHSD <- .prntbl(HSD[[1]], digits.d)
    for (i in 1:length(txHSD)) tx[length(tx)+1] <- txHSD[i]

    txhsd <- tx


    if (graphics) {
      if (!pdf) { 
        if (manage.gr) dev.set(which=4) 
      }
      else { 
        pdf.file <- "ANOVA_HSD.pdf"
        pdf(file=pdf.file, width=pdf.width, height=pdf.height)
      }

      plt.i <- plt.i + 1
      plt.title[plt.i] <- "95% family-wise confidence level"

      orig.params <- par(no.readonly=TRUE)
      on.exit(par(orig.params))
      par(mar=c(5.1,6.1,4.1,1.5))
      cex.axis <- .8; col.axis <- "gray30"; 
      plot(HSD, cex.axis=cex.axis, col.axis=col.axis, las=1)

      if (pdf) {
        dev.off()
        .showfile(pdf.file, "Tukey HSD chart")
        tx[length(tx)+1] <- ""
      }
    }
  }  # !brief

  return(list(
    title_des=title_des, txdes=txdes,
    title_basic=title_basic,
    txanv=txanv, txeft=txeft, txhsd=txhsd,
    i=plt.i, ttl=plt.title))

}
