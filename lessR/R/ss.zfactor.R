.ss.factor <-
function(x, by=NULL, brief=FALSE, digits.d=NULL, x.name, y.name=NULL,
         x.lbl=NULL, y.lbl=NULL, ...)  {


# construct a cross-tabs 
.prnfreq <- function(x, type, max.ln, max.c1, n.dash, ttl) {
  tx <- character(length = 0)

  # title
  tx[length(tx)+1] <- ttl 
  tx[length(tx)+1] <- .dash2(n.dash)
  tx[length(tx)+1] <- ""

  # col labels
  tx[length(tx)+1] <-  .fmtc(x.name, w=max.c1+3)
  tx[length(tx)+1] <-  format(y.name, width=max.c1, justify="left")
  w <- nchar(as.character(sum(x)))
  for (i in 1:ncol(x))
    tx[length(tx)] <- paste(tx[length(tx)], .fmtc(colnames(x)[i], w=max.ln[i]),
      sep="")

  # values
  for (i in 1:nrow(x)) {
    rwnm <- paste(" ", rownames(x)[i])
    tx[length(tx)+1] <-  format(rwnm, width=max.c1, justify="left")
    for (j in 1:ncol(x)) {
      if (type=="r") {
        tx[length(tx)] <- paste(tx[length(tx)], .fmt(x[i,j], d=3, w=max.ln[j]),
          sep="")
      }
      else if (type=="i")
        tx[length(tx)] <- paste(tx[length(tx)], .fmti(x[i,j], w=max.ln[j]),
          sep="")
    }
  }

  return(tx)
}  # end .prnfreq


  # begin
  # ---------------------------------

  # save ordered status before converting x to a table
  if (is.ordered(x) && is.null(by)) order.x <- TRUE else order.x <- FALSE
  if (is.ordered(by)) order.y <- TRUE else order.y <- FALSE

  # convert to table, with variable names, if needed
  if (!is.table(x) && !is.matrix(x)) {  # bc yields a table or matrix
    if (!is.null(by)) 
      x <- table(by,x, dnn=c(y.name,x.name)) 
    else
      x <- table(x, dnn=NULL)
  }

  # no title if two vars and no labels
  txttl <- ""
  dims <- length(dim(x))
  if (dims == 1 || (!is.null(x.lbl) || !is.null(y.lbl))) {  #  one var or labels
    txttl <- .title2(x.name, y.name, x.lbl, y.lbl, is.null(by), new.ln=FALSE)
    txttl[1] <- paste("\n", txttl[1])
  }

  # print table, chi-square analysis
  # -------------------------------------
  # two variables 
  if (!is.null(by) || is.matrix(x)) { 
    n.dim <- 2

    xx <- addmargins(x)

    # width of column 1
    if (!is.null(y.name))
      max.c1 <- nchar(y.name)
    else
      max.c1 <- 0
    for (i in 1:nrow(xx)) {
      c1 <- nchar(rownames(xx)[i])
      if (c1 > max.c1) max.c1 <- c1
    }
    max.c1 <- max.c1 + 2
    if (max.c1 < 5) max.c1 <- 5

    # width of data columns
    max.ln <- integer(length=0)
    for (i in 1:ncol(xx)) {
        ln.nm <- nchar(colnames(xx)[i])
      for (j in 1:nrow(xx)) {
        ln.vl <- nchar(as.character(xx[j,i]))
      }
        max.ln[i] <- max(ln.nm, ln.vl) + 1
        if (max.ln[i] < 4) max.ln[i] <- 4
    }

    # cell frequencies
    txfrq <- .prnfreq(xx, "i", max.ln, max.c1, n.dash=30,
                      ttl="Joint and Marginal Frequencies")

    tx <- character(length = 0)
    ch <- summary(as.table(x))
    if (!is.nan(ch$statistic)) {
      min.rc <- min(nrow(x)-1, ncol(x)-1)
      V <- sqrt(ch$statistic / (min.rc * ch$n.cases))
      txt <- ifelse(ch$parameter == 1, " (phi)", "") 
      txt <- paste("Cramer\'s V", txt, ":", sep="")
      tx[length(tx)+1] <- paste(txt, .fmt(V,3))
      tx[length(tx)+1] <- ""
      tx[length(tx)+1] <- paste("Chi-square Test", 
          ":  Chisq = ", .fmt(ch$statistic,3), ", df = ", ch$parameter,
          ", p-value = ", .fmt(ch$p.value,3), sep="")
      if (!ch$approx.ok) 
        tx[length(tx)+1] <- paste(">>> Low cell expected frequencies,",
            "so chi-squared approximation may not be accurate")
      tx[length(tx)+1] <- ""
    }
    else
      tx[length(tx)+1] <- paste(
          "Cross-tabulation table not well-formed, usually too many zeros\n",
          "Cramer's V and the chi-squared analysis not possible\n\n", sep="")
    txXV <- tx

    if (brief)
      return(list(n.dim=n.dim, txttl=txttl, txfrq=txfrq, txXV=txXV))


    # full analysis
    nan.flag <- FALSE

    for (i in 1:ncol(xx)) {
      if (max.ln[i] < 6) max.ln[i] <- 6
    }

    # cell proportions and marginals
    xx <- round(addmargins(prop.table(x)),3)
    txprp <- .prnfreq(xx, "r", max.ln, max.c1, n.dash=30,
                      ttl="Cell Proportions and Marginals")

    # cell proportions within each column
    x.col <- prop.table(x, margin=2)
    Sum <- numeric(ncol(x.col))
    for (i in 1:ncol(x.col)) {
      Sum[i] <- sum(x.col[,i])
      if (is.nan(Sum[i])) nan.flag <- TRUE
    }
    x.col2 <- round(rbind(x.col,Sum),3)
    names(dimnames(x.col2)) <- names(dimnames(x.col))

    txcol <- .prnfreq(x.col2, "r", max.ln, max.c1, n.dash=35,
                      ttl="Cell Proportions within Each Column")

    # cell proportions within each row
    x.row <- prop.table(x, margin=1)
    Sum <- numeric(nrow(x.row))
    for (i in 1:nrow(x.row)) {
      Sum[i] <- sum(x.row[i,])
      if (is.nan(Sum[i])) nan.flag <- TRUE
    }
    x.row2 <- round(cbind(x.row,Sum),3)
    names(dimnames(x.row2)) <- names(dimnames(x.row))

    txrow <- .prnfreq(x.row2, "r", max.ln, max.c1, n.dash=32,
                      ttl="Cell Proportions within Each Row")

    if (nan.flag)
      cat("\nNote: NaN results from all values missing for that cell or margin.\n",
                 "     so any division to compute a proportion is undefined.\n")

    return(list(n.dim=n.dim, txttl=txttl, txfrq=txfrq, txXV=txXV, txprp=txprp,
                txcol=txrow, txrow=txcol))  # back to ss or ss data frame

    # end full analysis

  }  # end two variable


  else {  # one variable
    n.dim <- 1
    if (length(names(x)) == sum(x)) {  # x is a vector of the counts
      if (length(x) > 100)
        cat("\nOnly the first 100 values listed.  To see all, use\n",
               "the  values  function.\n\n")
      nms <- character(length=0)
      for (i in 1:min(length(x), 100)) nms[i] <- names(x)[i]
      cat("\n")
      print(nms)
      cat("\n"); stop(call.=FALSE, "\n","------\n",
          "All values for ", x.name, " are unique, so analysis not meaningful\n\n",
          "Perhaps a row ID instead of a variable\n",
          "If so, use  row.names=n  option for Read, where n refers to the ",
          "nth column\n\n", sep="")
    }

    else {

      max.ln <- integer(length=0)
      for (i in 1:length(x)) {
        ln.nm <- nchar(names(x[i]))
        ln.vl <- nchar(as.character(x[i]))
        max.ln[i] <- max(ln.nm, ln.vl) + 1
        if (max.ln[i] < 6) max.ln[i] <- 6
      }

      tx <- character(length = 0)

      tx[length(tx)+1] <- format("", width=13)
      w <- nchar(as.character(sum(x)))
      for (i in 1:length(x))
        tx[length(tx)] <- paste(tx[length(tx)], .fmtc(names(x[i]), w=max.ln[i]))
      tx[length(tx)] <- paste(tx[length(tx)], .fmtc("Total", w=w+6))

      tx[length(tx)+1] <- "Frequencies: "
      for (i in 1:length(x))
        tx[length(tx)] <- paste(tx[length(tx)], .fmti(x[i], w=max.ln[i]))
      tx[length(tx)] <- paste(tx[length(tx)], .fmti(sum(x), w=w+6))

      tx[length(tx)+1] <- "Proportions: "
      sum.x <- sum(x)
      xp <- numeric(length=0)
      xp <- x/sum.x
      for (i in 1:length(x))
        tx[length(tx)] <- paste(tx[length(tx)], .fmt(xp[i], 3, max.ln[i]))
      tx[length(tx)] <- paste(tx[length(tx)], .fmtc("1.000", w=w+6))

      txcnt <- tx

      tx <- character(length = 0)
      ch <- suppressWarnings(chisq.test(x))  # provide own warning of small n
      tx[length(tx)+1] <- "Chi-squared test of null hypothesis of equal probabilities"
      tx[length(tx)+1] <- paste("  Chisq = ", .fmt(ch$statistic,3), ", df = ",
          ch$parameter, ", p-value = ", .fmt(ch$p.value,3), sep="")
      if (any(ch$expected < 5)) 
        tx[length(tx)+1] <- paste(">>> Low cell expected frequencies,",
            "so chi-squared approximation may not be accurate", "\n")
      txchi <- tx

      return(list(n.dim=n.dim, title=txttl, counts=txcnt, 
                  chi=txchi, freq=x, prop=xp))
    }
  }  # one variable

}
