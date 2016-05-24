.ss.numeric <-
function(x, by=NULL, digits.d=NULL, brief, ...) {

  # get variable labels if exist
  gl <- .getlabels()
  x.name <- gl$xn; x.lbl <- gl$xl;
  y.name <- gl$yn; y.lbl <- gl$yl

  max.char <- 0
  if (is.null(by))
    n.lines <- 1
  else {
    bu <- as.factor(unique(na.omit(by)))
    n.lines <- nlevels(bu)
    for (i in 1:nlevels(bu))  # largest level name
      if (nchar(levels(bu)[i]) > max.char) max.char <- nchar(levels(bu)[i])
    vectors <- split(x, by)
  }

  if (is.null(digits.d)) {
    dig.dec <- .max.dd(x) + 1
    if (dig.dec == 1) dig.dec <- 2
  }
  else
    dig.dec <- digits.d
  options(digits.d=dig.dec)
  if (dig.dec > 10  && is.null(digits.d)) {
    cat("\nThese data values contain ", dig.dec, " decimal digits. To enhance\n",
        "the readability of the output, only 4 decimal digits are\n",
        "displayed.  To customize this setting, use the digits.d  parameter.\n",
        "Example for Variables Y and X:  > ss(Y, by=X, digits.d=3)\n",
        sep="")
    dig.dec <- 4
  }


  # get maximum chars in 1st three columns
  max.lv <- 0; max.n <- 0; max.nm <- 0
  if (n.lines > 1) {
    for (i in 1:n.lines) {
      nch.lv <- nchar(as.character(levels(as.factor(by))[i]))
      nch.n <- nchar(as.character(sum(!is.na(vectors[[i]]))))
      nch.nm <- nchar(as.character(sum(is.na(vectors[[i]]))))
      if (nch.lv > max.lv) max.lv <- nch.lv
      if (nch.n > max.n) max.n <- nch.n
      if (nch.nm > max.nm) max.nm <- nch.nm
    }
  }
  else {
    if (max.lv < 3) max.lv <- 3
    max.n <- nchar(as.character(sum(!is.na(x))))
    max.nm <- nchar(as.character(sum(is.na(x))))
  }

  max.ln <- 0  # get max.ln, maximum length of the individual fields
  for (i in 1:n.lines) {
    if (n.lines == 1)
      xx <- x
    else {
      lv <- levels(as.factor(by))[i]
      xx <- vectors[[i]]
    } 
    m <- mean(xx, na.rm=TRUE)
    s <- sd(xx, na.rm=TRUE)
    if (is.na(s) || is.null(s))
      n.ln <- nchar(as.character(round(m, dig.dec))) + dig.dec
    else
      n.ln <- max(nchar(as.character(round(m, dig.dec))), 
                  nchar(as.character(round(s, dig.dec)))) + dig.dec
    if (n.ln > max.ln) max.ln <- n.ln
  }
  if (max.ln < 5) max.ln <- max.ln + 1
  if (max.ln < 10) max.ln <- max.ln + 1


  tx <- character(length = 0)

  # ------
  # output

  # first the title with any variable labels
  txlbl <- .title2(x.name, y.name, x.lbl, y.lbl, is.null(by))
  txlbl <- gsub("\n", "", txlbl)
  for (i in 1:length(txlbl)) tx[length(tx)+1] <- txlbl[i] 
  tx[length(tx)+1] <- ""

  # --------------------------------
  # the stats loop
  for (i in 1:n.lines) {
    if (n.lines == 1) {
      xx <- x
      p.lv <- ""
    }
    else {
      lv <- levels(as.factor(by))[i]
      xx <- vectors[[i]]
      p.lv <- format(lv, width=max.char+1)
    }

    m <- NA; s <- NA; mn <- NA; mx <- NA; md=NA
    sk <- NA; kt <- NA; q1 <- NA; q3 <- NA; qr <- NA
    # get the descriptive statistics
    n <- sum(!is.na(xx))
    n.miss <- sum(is.na(xx))
    xx <- na.omit(xx)
    if (n>0) m <- mean(xx)
    if (n>1) s <- sd(xx)
    # skewness:  adjusted Fisher-Pearson standardized moment coefficient
    if (n>2) {
      sk.coef <- n / ((n-1)*(n-2))
      sk.sum <- 0
      for (j in 1:n) sk.sum <- sk.sum + (( (xx[j]-m) / s)^3) 
      sk <- sk.coef * sk.sum
    }
    # kurtosis
    if (n>3) {
      kt.coef1 <- (n*(n+1)) / ((n-1)*(n-2)*(n-3))
      kt.coef2 <- 3 * ( ((n-1)^2) / ((n-2)*(n-3)) )
      kt.sum <- 0
      for (j in 1:n) kt.sum <- kt.sum + ( (xx[j]-m)^4 )
      kt <- ( kt.coef1 * (kt.sum/(var(xx)^2)) ) - kt.coef2
    }
    # order stats
    if (n > 0) {
      mn <- min(xx)
      q1 <- quantile(xx, probs=0.25)
      md <- median(xx)
      q3 <- quantile(xx, probs=0.75)
      mx <- max(xx)
      qr <- IQR(xx)
    }


    if (i == 1) { # heading labels 

      if (max.ln < 4) max.ln <- max.ln + 2
      if (max.ln < 8) max.ln <- max.ln + 1
      nbuf <- ifelse (n.lines == 1, 2, 4)

      n.lbl <- .fmtc("n", nchar(as.character(n))+nbuf+max.lv)
      miss.lbl <- .fmtc("miss", nchar(as.character(n.miss))+5)
      m.lbl <- .fmtc("mean", max.ln)
      s.lbl <- .fmtc("sd", max.ln)
      mn.lbl <- .fmtc("min", max.ln)
      md.lbl <- .fmtc("mdn", max.ln)
      mx.lbl <- .fmtc("max", max.ln)

      if (brief)
        tx[length(tx)+1] <- paste(n.lbl, miss.lbl, m.lbl, s.lbl, mn.lbl,
                                 md.lbl, mx.lbl) 
      else {
        sk.lbl <- .fmtc("skew", max.ln)
        kt.lbl <- .fmtc("krts", max.ln)
        q1.lbl <- .fmtc("qrt1", max.ln)
        q3.lbl <- .fmtc("qrt3", max.ln)
        qr.lbl <- .fmtc("IQR", max.ln)
        tx[length(tx)+1] <- paste(n.lbl, miss.lbl, m.lbl, s.lbl, sk.lbl,
              kt.lbl, mn.lbl, q1.lbl, md.lbl, q3.lbl, mx.lbl, qr.lbl)
      }
    }  # end first line, labels

    # write values
    lvl <- .fmtc(p.lv, max.lv)
    n.c <- .fmti(n, max.n+1)
    miss.c <- .fmti(n.miss, max.nm+5)

    if (n == 0) 
      tx[length(tx)+1] <- paste(lvl, n.c, miss.c)

    else if (n == 1) {
        m.c <- .fmt(m, dig.dec, max.ln)
        tx[length(tx)+1] <- paste(lvl, n.c, miss.c, m.c)
    }
    else if (n > 1) {
      m.c <- .fmt(m, dig.dec, max.ln)
      s.c <- .fmt(s, dig.dec, max.ln)
      mn.c <- .fmt(mn, dig.dec, max.ln)
      md.c <- .fmt(md, dig.dec, max.ln)
      mx.c <- .fmt(mx, dig.dec, max.ln)

      if (brief)
        tx[length(tx)+1] <- paste(lvl, n.c, miss.c, m.c, s.c, mn.c, md.c, mx.c)

      else {
        sk.c <- .fmt(sk, dig.dec, max.ln)
        kt.c <- .fmt(kt, dig.dec, max.ln)
        q1.c <- .fmt(q1, dig.dec, max.ln)
        q3.c <- .fmt(q3, dig.dec, max.ln)
        qr.c <- .fmt(qr, dig.dec, max.ln)
        tx[length(tx)+1] <- paste(lvl, n.c, miss.c, m.c, s.c, sk.c, kt.c,
             mn.c, q1.c, md.c, q3.c, mx.c, qr.c)
      }
    }

  }  # for each line

     if (n.lines == 1)
       return(list(txlbl=txlbl, tx=tx, n=n, n.miss=n.miss, m=m, s=s, sk=sk,
                   kt=kt, mn=mn, q1=q1, md=md, q3=q3, mx=mx, qr=qr))
     else
       return(list(tx=tx))  # contains title with var labels and stats

}
