corCFA <- 
function(mimm=NULL, x=mycor, data=mydata, fac.names=NULL, 

         Rmd=NULL, explain=getOption("explain"),
         interpret=getOption("interpret"), results=getOption("results"),

         labels=c("include", "exclude", "only"),

         min.cor=.10, min.res=.05, iter=50, grid=TRUE, 

         resid=TRUE, item.cor=TRUE, sort=TRUE,

         main=NULL, heat.map=TRUE, bottom=3, right=3, 

         pdf.file=NULL, pdf.width=5, pdf.height=5,

         F1=NULL, F2=NULL, F3=NULL, F4=NULL, F5=NULL,
         F6=NULL, F7=NULL, F8=NULL, F9=NULL, F10=NULL,
         F11=NULL, F12=NULL,

         fun.call=NULL) {

  if (exists(deparse(substitute(data)), where=.GlobalEnv)) 
    dname <- deparse(substitute(data))
  else
    dname <- NULL

  if (is.null(fun.call)) fun.call <- match.call()

  labels <- match.arg(labels)

  if (!is.null(pdf.file))
    if (!grepl(".pdf", pdf.file)) pdf.file <- paste(pdf.file, ".pdf", sep="")

  max.fname <- 0
  if (!is.null(fac.names)) {
    for (i in 1:length(fac.names))
      if (nchar(fac.names[i]) > max.fname) max.fname <- nchar(fac.names[i])
  }

  if (is.null(dname)  &&  !is.null(Rmd)) {
    cat("\n"); stop(call.=FALSE, "\n","------\n",
      "Need to read from the data table (frame) to generate a Rmd. \n\n")
  }

  if (labels!="only") {
    cor.nm <- deparse(substitute(x))
    .cor.exists(cor.nm)  # see if matrix exists in one of the 3 locations
    if (class(x) == "out_all")
      x <- eval(parse(text=paste(cor.nm, "$cors", sep="")))  # go to $cors 
  }
  else  # if only labels, then need the data 
    if (is.null(dname)) {
      cat("\n"); stop(call.=FALSE, "\n","------\n",
        "No data table (frame) exists from which to get the labels\n\n")
    }


  df.name <- deparse(substitute(data))
  options(dname = df.name)

  NFmax <- 12


  # translate variable names into column positions
  if (labels!="only") {
    NVOld <- as.integer(nrow(x))
    vars.all <- as.list(seq_along(as.data.frame(x)))
    names(vars.all) <- names(as.data.frame(x))
    nm <- dimnames(x)[[1]]
  }
  else {
    NVOld <- nrow(data)
    vars.all <- as.list(seq_along(data))
    names(vars.all) <- names(data)
    nm <- names(data)
  }


  if (!is.null(mimm)) {
    nm.mimm <- deparse(substitute(mimm))
    
    c <- ""  # strip blanks and tildes
    for (i in 1:nchar(mimm)) {
      s <- substr(mimm,i,i) 
      if (s!=" " && s!="~" ) c <- paste(c,s, sep="") 
    }

    c <- strsplit(c, "\n")[[1]]  # split in factor = strings
    d <- c[which(nchar(c)>0)]  # remove empty strings

    # split into sets of factor name and vars
    NF <- length(d)
    fac.names <- character(length=NFmax)
    vars <- character(length=NFmax)
    for (i in 1:NF) {
      f <- strsplit(d[i], "=")[[1]]
      fac.names[i] <- f[1]
      vars[i] <- gsub("+", ",", f[2], fixed=TRUE)
      vars[i] <- paste("c(", vars[i], ")", sep="")
    }

    F1n <- eval(parse(text=vars[1]), vars.all, parent.frame())
    F2n <- eval(parse(text=vars[2]), vars.all, parent.frame())
    F3n <- eval(parse(text=vars[3]), vars.all, parent.frame())
    F4n <- eval(parse(text=vars[4]), vars.all, parent.frame())
    F5n <- eval(parse(text=vars[5]), vars.all, parent.frame())
    F6n <- eval(parse(text=vars[6]), vars.all, parent.frame())
    F7n <- eval(parse(text=vars[7]), vars.all, parent.frame())
    F8n <- eval(parse(text=vars[8]), vars.all, parent.frame())
    F9n <- eval(parse(text=vars[9]), vars.all, parent.frame())
    F10n <- eval(parse(text=vars[10]), vars.all, parent.frame())
    F11n <- eval(parse(text=vars[11]), vars.all, parent.frame())
    F12n <- eval(parse(text=vars[12]), vars.all, parent.frame())

    for (i in 1:NFmax) {
      fnum <- eval(parse(text=paste("F", toString(i), "n", sep="")))
      if (nchar(vars[i] == 0)) fnum <- NULL 
    }
  }

  else {
    F1n <- eval(substitute(F1), vars.all, parent.frame())
    F2n <- eval(substitute(F2), vars.all, parent.frame())
    F3n <- eval(substitute(F3), vars.all, parent.frame())
    F4n <- eval(substitute(F4), vars.all, parent.frame())
    F5n <- eval(substitute(F5), vars.all, parent.frame())
    F6n <- eval(substitute(F6), vars.all, parent.frame())
    F7n <- eval(substitute(F7), vars.all, parent.frame())
    F8n <- eval(substitute(F8), vars.all, parent.frame())
    F9n <- eval(substitute(F9), vars.all, parent.frame())
    F10n <- eval(substitute(F10), vars.all, parent.frame())
    F11n <- eval(substitute(F11), vars.all, parent.frame())
    F12n <- eval(substitute(F12), vars.all, parent.frame())

    # get NF, number of factors
    NF <- 0
    for (i in 1:NFmax) {
      fnum <- eval(parse(text=paste("F", toString(i), "n", sep="")))
      if (!is.null(fnum)) NF <- NF + 1
    }
  }

  Label <- c(F1n,F2n,F3n,F4n,F5n,F6n,F7n,F8n,F9n,F10n,F11n,F12n)

  if (NF == 0) {
    cat("\n"); stop(call.=FALSE, "\n","------\n",
      "Number of Factors: ", NF, "\n",
      "Need to specify some factors.", "\n\n",
      "For example, F1=c(...), F2=c(...), etc.\n\n")
  }

  if (!is.null(fac.names)) if (length(fac.names) < NF) {
    cat("\n"); stop(call.=FALSE, "\n","------\n",
      "Only ", length(fac.names), " factor names entered for ", NF, " factors\n\n")
  }

  # get the ordinal position of the first and last vars in Group i
  # get NItems
  LblCut <- matrix(nrow=NF, ncol=2)
  NItems <- 0
  for (i in 1:NF) {
    LblCut[i,1] <- NItems + 1
    cFac <- eval(parse(text=paste("F", toString(i), "n", sep="")))
    if (length(cFac) == 0) {
      cat("\n"); stop(call.=FALSE, "\n","------\n",
          "Factor Number ", i, " has no items.\n",
          "Each factor must have at least one item.\n\n")
    }
    NItems <- NItems + length(cFac)
    LblCut[i,2] <- NItems
  }

  for (i in 1:NItems) {
    if (Label[i] > NVOld) {
      cat("\n"); stop(call.=FALSE, "\n","------\n",
        "Number of items in correlation matrix: ", NVOld, "\n",
        "Item number in Factor specification: ", Label[i], "\n\n",
        "Specified item does not exist in this correlation matrix.\n\n")
    }
  }

  # display labels by factor if labels=="only"
  txlbl <- ""
  if (labels == "only") {
    tx <- character(length = 0)
    tx[length(tx)+1] <- ""

    for (i in 1:NF) {
      tx[length(tx)+1] <- paste("F", toString(i), sep="")
      if (max.fname > 0)
        tx[length(tx)] <- paste(tx[length(tx)], " - ", fac.names[i], sep="")
      tx[length(tx)+1] <- .dash2(30)
      for (j in LblCut[i,1]:LblCut[i,2]) {
        options(xname = nm[Label[j]])
        tx[length(tx)+1] <- paste(nm[Label[j]], ": ", .getlabels()$xl, sep="")
      }
      tx[length(tx)+1] <- ""
    }

    txlbl <- tx
    class(txlbl) <- "out_piece"
    output <- list(out_labels=txlbl)
  }  # end labels only


  else { # proceed with the analysis

  # --------------------------------------------------------
  # re-order R matrix

  outR <- x[Label,Label]
  nm.new <- colnames(outR)

  # get width of largest variable (item) name
  cc <- as.character(dimnames(outR)[[1]])
  max.chr <- 0
  for (i in 1:NItems)  
    if (nchar(cc[i]) > max.chr) max.chr <- nchar(cc[i])
  if (max.chr < 4) max.chr <- 4


  # --------------------------------------------------------
  # MIMM CFA

  # expand R matrix to include rows/cols for factors
  rr <- matrix(rep(0, NF*NItems), nrow=NF)
  cc <- matrix(rep(0, NF*(NItems+NF)), nrow=(NItems+NF))
  outR <- cbind(rbind(outR,rr),cc)

  alpha <- numeric(length=NF)
  omega <- numeric(length=NF)

  out <- .mimm(outR, LblCut, NItems, NF, iter)

  nmF <- character(length=NF)
  for (i in 1:NF) nmF[i] <- paste("F", toString(i), sep="")

  NVTot <- NItems + NF

  # assign names
  nm  <- character(length=NVTot)
  nm <- c(nm.new, nmF)
  dimnames(out$R) <- list(nm, nm)


  # --------------------------------------------------------
  # Sort within each group by the group factor loading

  if (sort) {

    # get new ordering, factor by factor
    pt <- numeric(length=NItems)
    newLabel <- numeric(length=NItems)
    for (ifac in 1:NF) {
      n1 <- LblCut[ifac,1]
      n2 <- LblCut[ifac,2]
      irow <- NItems + ifac
      for (j in n1:n2) pt[j] <- out$R[irow, j]
      o <- order(pt[n1:n2], decreasing=TRUE)
      for (i in 1:(n2-n1+1)) newLabel[n1-1+i] <- Label[n1-1+o[i]]
    }
    Label <- newLabel

    outR <- x[Label,Label]

    nm.new <- colnames(outR)

    # expand R matrix to include rows/cols for factors
    rr <- matrix(rep(0, NF*NItems), nrow=NF)
    cc <- matrix(rep(0, NF*(NItems+NF)), nrow=(NItems+NF))
    outR <- cbind(rbind(outR,rr),cc)

    # MIMM CFA
    alpha <- numeric(length=NF)
    omega <- numeric(length=NF)

    out <- .mimm(outR, LblCut, NItems, NF, iter)

    nmF <- character(length=NF)
    if (is.null(fac.names))
      for (i in 1:NF) nmF[i] <- paste("F", toString(i), sep="")
    else
     for (i in 1:NF) nmF[i] <- fac.names[i] 

    # assign names
    nm  <- character(length=NVTot)
    nm <- c(nm.new, nmF)
    dimnames(out$R) <- list(nm, nm)
  }

  # --------------------------------------------------------
  if (heat.map) {
    if (is.null(main)) main <- "Item Correlations/Communalities"
   .corcolors(out$R, NItems, main, bottom, right, diag=NULL,
              pdf.file, pdf.width, pdf.height)
  }



  title_scales <- "  FACTOR / SCALE COMPOSITION"

  txlbl <- ""
  tx <- character(length = 0)

  anyLabels <- FALSE
  for (i in 1:NItems) {
    options(xname = nm.new[i])
    if (!is.null(.getlabels()$xl)) anyLabels <- TRUE
  }

  for (i in 1:NF) {
    tx[length(tx)+1] <- paste("F", toString(i), sep="")
    if (max.fname > 0)
      tx[length(tx)] <- paste(tx[length(tx)], " - ", fac.names[i], sep="")
    if (!anyLabels) { # horizontal
      for (j in LblCut[i,1]:LblCut[i,2])
      tx[length(tx)] <- paste(tx[length(tx)], " ", nm.new[j])
    }
    else {  # vertical
      tx[length(tx)+1] <- .dash2(30)
      for (j in LblCut[i,1]:LblCut[i,2]) {   
        options(xname = nm.new[j])
        tx[length(tx)+1] <- paste(nm.new[j], ": ", xW(.getlabels()$xl), sep="")
      }
    }
    if (i < NF) tx[length(tx)+1] <- ""
  }

  txlbl <- tx



  title_rel <- "  RELIABILITY ANALYSIS"

  txrel <- ""
  tx <- character(length = 0)

  if (iter > 0) {
    buf <- ifelse(max.fname > 0, max.fname+6, 3)
    tx[length(tx)+1] <- paste(" Scale", paste(rep(" ", buf), collapse=""),
      "Alpha    Omega", sep="")
    tx[length(tx)+1] <- .dash2(23)
    if (max.fname > 0) tx[length(tx)] <- paste(tx[length(tx)],
      .dash2(max.fname+3), sep="")
  }
  else {
    tx[length(tx)+1] <- " Scale   Alpha"
    tx[length(tx)+1] <- .dash2(14)
  }
  for (i in 1:NF) {
    Fnm <- paste("F", as.character(i), sep="")
    tx[length(tx)+1] <- paste("  ", Fnm) 
    if (max.fname > 0)
      tx[length(tx)] <- paste(tx[length(tx)], " - ",
        .fmtc(fac.names[i], w=max.fname, j="left"), sep="")
    tx[length(tx)] <- paste(tx[length(tx)], " ", .fmt(out$Alpha[i],3, w=6)) 
    if (iter > 0)
      tx[length(tx)] <- paste(tx[length(tx)], " ", .fmt(out$Omega[i],3, w=6))
    else
      out$Omega <- NULL
  }
  txrel <- tx



  title_sol <- "  SOLUTION"

  txind <- ""
  tx <- character(length = 0)

  if (iter > 0 ) {
    MaxLbl <- NItems

    buf <- max.chr - 4
    if (buf < 0) buf <- 0


    if (is.null(options()$knitr.in.progress))
       tx[length(tx)+1] <- paste('Indicator Analysis\n')
    tx[length(tx)+1] <- paste('Fac', ' Indi', .fmtc(" ",buf+1), 'Pat', '   Unique',
       ' Factors with which an indicator correlates too')
    tx[length(tx)+1] <- paste('tor', ' cator', .fmtc("",buf), 'tern', '  ness',
        '   highly, and other indicator diagnostics.')
    tx[length(tx)+1] <- .dash2(75)

    for (IFac in 1:NF) {
      tx[length(tx)+1] <- ""
      Fnm <- paste("F", as.character(IFac), sep="")

      for (Item in LblCut[IFac,1]:LblCut[IFac,2]) {
        Bad <- integer(length=0)
        Lam <- out$R[NItems+IFac,Item]
        Unique <- 1 - Lam**2

        if (Lam <= 0)
          unq <- "   xxxx"
        else
          unq <- .fmt(Unique,3,7)

        tx[length(tx)+1] <- paste(.fmtc(Fnm,3), .fmtc(nm.new[Item],max.chr),
              .fmt(Lam,3,7), unq, "   ")

        if (Lam>0 && Unique>0) {
          for (I in 1:NF)
            if (abs(out$R[NItems+I,Item]) > Lam) Bad[length(Bad)+1] <- I
          if (length(Bad) > 0) for (IBad in 1:length(Bad))
            tx[length(tx)] <- paste(tx[length(tx)], paste("F", Bad[IBad], " ", sep=""))
        }
        else if (Lam <= 0)
          tx[length(tx)] <- paste(tx[length(tx)], '** Negative Loading on Own Factor **')
        else if (Unique <= 0) {
          if (LblCut[IFac,2]-LblCut[IFac,1] > 0)
            tx[length(tx)] <- paste(tx[length(tx)], '** Improper Loading **')
          else
            tx[length(tx)] <- paste(tx[length(tx)], '** Factor Defined by Only One Item **')
        }

        Bad <- rep(0, NItems)
      }  # each item within a factor

    }  # each factor

    txind <- tx
  }


  # --------------------------------------------------------
  # Solution
  txsol <- ""
  tx <- character(length = 0)

  if (iter > 0)
    if (is.null(options()$knitr.in.progress))
      tx[length(tx)+1] <- "Factor / Item Correlations \n"
  else
    if (is.null(options()$knitr.in.progress))
      tx[length(tx)+1] <- "Item-Scale and Scale-Scale Correlations\n"

  # print the solution
  if(grid) boundary <- LblCut[,2] else boundary <-  NULL
  if (item.cor)
    txcrs <- .prntbl(out$R, 2, cut=min.cor, cc=NULL, cors=TRUE, bnd=boundary)
  else
    txcrs <- .prntbl(out$R[1:NVTot,(NItems+1):NVTot], 2, cut=min.cor,
                     cc=NULL, cors=TRUE, bnd=boundary)
  for (i in 1:length(txcrs)) tx[length(tx)+1] <- txcrs[i]

  txsol <- tx




  title_res <- "  RESIDUALS"

  txres <- ""
  txrst <- ""
  if (resid) {
    tx <- character(length = 0)

    phi <- out$R[(NItems+1):(NItems+NF), (NItems+1):(NItems+NF)]
    lambda <- matrix(0, nrow=NItems, ncol=NF)
    iFac <- 1
    for (i in 1:NItems) {
      if (i > LblCut[iFac,2]) iFac <- iFac + 1 
      lambda[i, iFac] <- out$R[i, NItems+iFac]
    }
    est <- lambda %*% phi %*% t(lambda)
    diag(est) <- 1.0
    est <- round(est, 5)
    colnames(est) <- row.names(out$R[1:NItems,1:NItems])
    rownames(est) <- colnames(est)

    res <- out$R[1:NItems,1:NItems] - est
    diag(res) <- 0
    res <- round(res, 5)

    if (is.null(options()$knitr.in.progress))
      tx[length(tx)+1] <- "Item residuals\n"

    txcrs <- .prntbl(res, 2, cut=min.res, cc=NULL, cors=TRUE, bnd=boundary)
    for (i in 1:length(txcrs)) tx[length(tx)+1] <- txcrs[i]

    txres <- tx


    tx <- character(length = 0)

    # sum of squares, sum of abs
    if (is.null(options()$knitr.in.progress))
      tx[length(tx)+1] <- "Residual summaries\n"

    tx[length(tx)+1] <- paste(.fmtc(" ", max.chr+2), "Sum of    Average", "\n",
        .fmtc(" ", max.chr+2), "Squares   Abs Value", "\n",
        .fmtc(" ", max.chr+2), "-------   ---------", sep="")

    cc <- as.character(dimnames(res)[[1]])
    res.avg <- double(length=NItems)

    ssq.tot <- 0
    abv.tot <- 0
    abv.all <- 0
    for (i in 1:NItems) {
      ssq <- 0
      abv <- 0
      for (j in 1:NItems) {
        ssq <- ssq + res[i,j]^2
        abv <- abv + abs(res[i,j])
        abv.all <- abv.all + abs(res[i,j])
      }
      ssq.tot <- ssq.tot + ssq
      res.avg[i] <- abv / (NItems - 1)
      tx[length(tx)+1] <- paste(.fmtc(cc[i],max.chr), "  ", .fmt(ssq,3), "  ", .fmt(res.avg[i],3))
    }
    abv.all.tot <- abv.all / (NItems^2 - NItems)
    tx[length(tx)+1] <- ""
    tx[length(tx)+1] <- paste("Total sum of squares for all items:", .fmt(ssq.tot,3), "\n")
    tx[length(tx)+1] <- paste("Average absolute residual w/o the diagonal:", .fmt(abv.all.tot,3), "\n\n")

    txrst <- tx
  }
  else { # no residuals
    res <- NULL
    est <- NULL
  }





  title_lvn <- "  LAVAAN SPECIFICATION"
  txlvn <- ""

  if (iter > 0) {
    tx <- character(length = 0)

    if (is.null(mimm)) nm.mimm <- "MeasModel"

    tx[length(tx)+1] <- paste(nm.mimm, " <-")

    tx[length(tx)+1] <- paste("\"")
    for (i in 1:NF) {
      if (max.fname > 0) nmF[i] <- fac.names[i]
      tx[length(tx)+1] <- paste("   ", nmF[i], " =~", sep="")

      for (j in LblCut[i,1]:LblCut[i,2]) {
        if (j == LblCut[i,1])
          tx[length(tx)] <- paste(tx[length(tx)], " ", nm.new[j], sep="")
        else  
          tx[length(tx)] <- paste(tx[length(tx)], "+", nm.new[j])
      }
    }
    tx[length(tx)+1] <- paste("\"\n")

    tx[length(tx)+1] <- paste("library(lavaan)")
    tx[length(tx)+1] <- paste("fit <- lavaan::cfa(", nm.mimm, ", data=mydata,",
      " std.ov=TRUE, std.lv=TRUE)", sep="")
    tx[length(tx)+1] <- "summary(fit, fit.measures=TRUE)"
    tx[length(tx)+1] <- ""

    tx[length(tx)+1] <- "--------"
    tx[length(tx)+1] <- paste(">>> The preceding code fits the model from",
      "data frame:  mydata")
    tx[length(tx)+1] <- paste(">>> To access the correlation matrix",
      "directly without the data")
    tx[length(tx)+1] <- paste(">>>   use the following fit statement instead.\n")
    tx[length(tx)+1] <- paste("fit <- lavaan::cfa(", nm.mimm, 
     ", sample.cov=mycor$cors,", " sample.nobs=nnn, std.lv=TRUE)\n", sep="")
    tx[length(tx)+1] <- ">>>   mycor: name of correlation matrix"
    tx[length(tx)+1] <- ">>>   nnn: numeric, number of observations"

    txlvn <- tx
  }



  # R Markdown
  txkfl <- ""
  if (!is.null(Rmd)) {
    if (!grepl(".Rmd", Rmd)) Rmd <- paste(Rmd, ".Rmd", sep="")
    txknt <- .corfa.Rmd(mimm, nm.mimm, dname, fun.call, NItems, NF,
                 iter, item.cor, explain, interpret, results)
    cat(txknt, file=Rmd, sep="\n")
    txkfl <- .showfile2(Rmd, "R Markdown instructions")
  }


  class(title_scales) <- "out_piece"
  class(txlbl) <- "out_piece"
  class(title_rel) <- "out_piece"
  class(txrel) <- "out_piece"
  class(title_sol) <- "out_piece"
  class(txind) <- "out_piece"
  class(txsol) <- "out_piece"
  class(title_res) <- "out_piece"
  class(txres) <- "out_piece"
  class(txrst) <- "out_piece"
  class(title_lvn) <- "out_piece"
  class(txlvn) <- "out_piece"
  class(txkfl) <- "out_piece"

  output <- list(
    call=fun.call,

    out_title_scales=title_scales, out_labels=txlbl,
    out_title_rel=title_rel, out_reliability=txrel,
    out_title_solution=title_sol, out_indicators=txind, out_solution=txsol,
    out_title_residuals=title_res, out_residuals=txres, out_res_stats=txrst,
    out_title_lavaan=title_lvn, out_lavaan=txlvn,  out_Rmd=txkfl,

    ff.cor=out$R[(NItems+1):NVTot,(NItems+1):NVTot],
    if.cor=out$R[1:NItems,(NItems+1):NVTot],
    diag.cor=diag(out$R[1:NItems,1:NItems]),
    alpha=out$Alpha,
    omega=out$Omega,
    pred=est,
    resid=res
    )

  }  # proceed with analysis

  # --------------------------------------------------------
  # return

  class(output) <- "out_all"

  return(output)
}

