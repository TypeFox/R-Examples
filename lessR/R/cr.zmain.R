.cr.main <-
function(x, y, brief, ...) {

  # get variable labels if exist
  gl <- .getlabels()
  x.name <- gl$xn; x.lbl <- gl$xl;
  y.name <- gl$yn; y.lbl <- gl$yl

  if (!is.factor(x)) {

    ct <- cor.test(x,y, ...)

    if (ct$method == "Pearson's product-moment correlation") {
      c.type <- "pearson" 
      sym <- "r"
    }
    else if (ct$method == "Spearman's rank correlation rho") {
      c.type <- "spearman" 
      sym <- names(ct$estimate)
    }
    else if (ct$method == "Kendall's rank correlation tau") {
      c.type <- "kendall" 
      sym <- names(ct$estimate)
    }
    if (c.type == "pearson") sym.pop <- "correlation" else sym.pop <- sym

    if (ct$alternative == "two.sided")
      h.txt <- "not equal to"
    else if (ct$alternative == "less")
      h.txt <- "less than"
    else if (ct$alternative == "greater")
      h.txt <- "greater than"

    # background
    tx <- character(length = 0)
    if (!brief) {
      txt <- paste("Correlation Analysis for Variables", x.name, "and", y.name)
      tx[length(tx)+1] <- txt
      tx[length(tx)+1] <- .dash2(nchar(txt))
    }

    tx[length(tx)+1] <- ""
    tx[length(tx)+1] <- paste(">>> ",ct$method, sep="")

    if (!is.null(x.lbl) || !is.null(y.lbl)) {
      tx[length(tx)+1] <- ""
      if (!is.null(x.lbl))
        tx[length(tx)+1] <- paste(x.name, ": ", as.character(x.lbl), sep="")
      else
        tx[length(tx)+1] <- paste(x.name)
      if (!is.null(y.lbl))
        tx[length(tx)+1] <- paste(y.name, ": ", as.character(y.lbl), sep="")
      else
        tx[length(tx)+1] <- paste(y.name)
    }
   
      n.pair <- sum(!is.na(x - y))  # number of points after pairwise deletion
      tx[length(tx)+1] <- ""
      tx[length(tx)+1] <- paste("Number of paired values with neither",
        "missing, n =", n.pair)

    if (!brief) {
      n.del <- sum(is.na(x - y))  # number of pairwise deleted observations
      tx[length(tx)+1] <- paste("Number of cases (rows of data) deleted:",
        n.del)
    }

    txb <- tx

    # descriptive
    tx <- character(length = 0)
    if (!brief) {
      if (c.type == "pearson") {
        covr <- cov(x, y, use="pairwise.complete.obs")
        tx[length(tx)+1] <- paste("Sample Covariance: s =", .fmt(covr,3))
      tx[length(tx)+1] <- ""
      }
      tx[length(tx)+1] <- paste("Sample Correlation: ", sym, " = ",
        .fmt(ct$estimate,3), sep="")
    }
    else
      tx[length(tx)+1] <- paste("Sample Correlation of ", x.name, " and ",
        y.name, ": ", sym, " = ", .fmt(ct$estimate,3), sep="")

    txd <- tx


    # inferential
    tx <- character(length = 0)

    tx[length(tx)+1] <- paste("Alternative Hypothesis: True", sym.pop,
      "is", h.txt, "0")
    tx[length(tx)+1] <- paste("  ", names(ct$statistic), "-value: ",
      .fmt(ct$statistic,3), sep="")
    if (c.type == "pearson")
        tx[length(tx)] <- paste(tx[length(tx)], ",  df: ",
          ct$parameter, sep="") 
    tx[length(tx)] <- paste(tx[length(tx)], ",  p-value: ",
      .fmt(ct$p.value,3), sep="")

    if (c.type == "pearson") {
      tx[length(tx)+1] <- ""
      tx[length(tx)+1] <- "95% Confidence Interval of Population Correlation"
      tx[length(tx)+1] <- paste("  Lower Bound:", .fmt(ct$conf.int,3)[1])
      tx[length(tx)] <- paste(tx[length(tx)], "     Upper Bound:",
        .fmt(ct$conf.int,3)[2])
      df <- round(ct$parameter,0)
      lb <- round(ct$conf.int,3)[1]
      ub <- round(ct$conf.int,3)[2]
    }
    else {
      df <- NA; lb <- NA; ub <- NA
    }

    txi <- tx

  }

  return(list(txb=txb, txd=txd, txi=txi,
    r=round(ct$estimate,3), tvalue=round(ct$statistic,3),
    df=df, pvalue=round(ct$p.value,3), lb, ub))

}
