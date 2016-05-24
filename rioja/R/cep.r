read.CEP <- function(fName, mValue=-99.9, impZero=0) {
#  if (.Machine$sizeof.pointer == 8) {
#     print("This function is not yet working on 64-bit OSs.\nTry read.cep in the vegan package as an alternative.")
#     return (NULL)
#  }
  if (!is.character(fName))
    stop("Expecting filename as first argument")
  if (!is.numeric(mValue))
    stop("Argument mValue must be numeric")
  if (!is.numeric(impZero) && !is.na(impZero))
    stop("Argument impZero must be numeric or NA")
  ret <- .Call("ReadCornellFile", fName, as.double(mValue), as.double(impZero), PACKAGE="rioja")
  if (!is.null(ret[[1]])) {
    out <- as.data.frame(ret[[1]])
    rownames(out) <- ret[[2]]
  }
  else {
    errMsg <- ret[[3]]
    if (substring(errMsg, 1, 18) == "Error reading data") {
       errMsg <- paste(errMsg, "\nIf the file came from Windows this problem may be caused by the wrong end-of-line characters for your operating system.\nOpen the file in a text editor and convert if necessary")
    }
    stop(errMsg)
  }
  if (any(ret$Summary > 0)) {
    eMess <- paste("File", fName, "contains")
    if (ret$Summary[1] > 0)
        eMess <- paste(eMess, ret$Summary[1], "missing value(s)")
    if (ret$Summary[1] > 0 && ret$Summary[2] > 0)
        eMess <- paste(eMess, "and")
    if (ret$Summary[2] > 0)
        eMess <- paste(eMess, ret$Summary[2], "column(s) with no data (which have been removed)")
    eMess <- paste(eMess, ".", sep="")
    warning(eMess)
  }
  return(out)
}

write.CEP <- function(x, fName, fTitle=fName, type = "condensed", nSig=5, nCouplets=4, mValue=-99.9) {
  x <- as.matrix(x)
  if (!is.numeric(x))
    stop("output matric must be numeric")
  if (!is.character(fName))
    stop("Expecting filename as first argument")
  itype <- pmatch(type, c("condensed", "full"), nomatch=0) - 1
  if (itype < 0)
    stop("Argument type must be \"condensed\" or \"full\" (partial match allowed).")
  if (!is.numeric(nSig) && !is.numeric(nCouplets))
    stop("Arguments nSig and nCouplets must be numeric")
  if (is.null(colnames(x)))
     colnames(x) <- 1:ncol(x)
  if (is.null(rownames(x)))
     rownames(x) <- 1:nrow(x)
  rN <- rownames(x)
  cN <- colnames(x)
  if (any(nchar(rN)>8)) {
     rownames(x) <- make.cepnames(rN)
     warning("Some sample names were longer than 8 characters and have been abbreviated - see returned list.")
  }
  if (any(nchar(cN)>8)) {
    colnames(x) <- make.cepnames(cN)
    warning("Some species names were longer than 8 characters and have been abbreviated - see returned list.")
  }
  ret <- .Call("WriteCornellFile", x, fName, fTitle, as.integer(itype), as.integer(nSig), as.integer(nCouplets), as.integer(nSig), as.integer(nCouplets), as.double(mValue), PACKAGE="rioja")
  if (nchar(ret)>0) {
    if (substr(ret, 1, 6) == "Cannot")
      stop (ret)
    else
      warning(ret)
  }
  rN <- cbind(original=rN, written=rownames(x))
  cN <- cbind(original=cN, written=colnames(x))
  invisible (list(site.names=rN, sp.names=cN))
}

read.Tilia <- function(fName) {
  if (!is.character(fName))
    stop("Expecting filename as first argument")
  ret <- .Call("ReadTiliaFile", fName, PACKAGE="rioja")
  if (!is.null(ret[[1]])) {
    df <- as.data.frame(ret[[1]])
    rownames(df) <- make.unique(ret[[2]])
    vars=data.frame(CodeNames=ret[[3]], CodeNums=ret[[4]], FullNames=ret[[5]], Sums=ret[[6]])
    levels = data.frame(Names=make.unique(ret[[2]]), Depths = ret[[7]])
    mt <- na.omit(grep(toupper("Chron*"), toupper(colnames(df))))
    if (length(mt) > 0) {
       ch <- df[, mt]
       df <- df[, -mt]
       levels <- cbind(levels, ch)
    }
    out <- list(data=df, vars=vars, levels=levels)
  }
  else {
    out <- NULL
    stop(ret[[8]])
  }
  return(out)
}

#read.C2Model <- function(fName) {
#    if (requireNamespace("RODBC", quietly=TRUE)==FALSE) {
#       stop("This function requires package RODBC")
#    }
#    channel <- RODBC::odbcConnectExcel(fName)
#    tabs <- RODBC::sqlTables(channel, errors=TRUE)
#    tabs <- tabs[tabs$TABLE_TYPE == "TABLE", ]
#    ntabs <- nrow(tabs)
#    X <- vector("list", length=ntabs)
#    names(X) <- tabs$TABLE_NAME
#    for (i in 1:ntabs) {
#       X1 <- RODBC::sqlFetch(channel, tabs$TABLE_NAME[i], colnames=FALSE, rownames=FALSE)
#       if (tabs$TABLE_NAME[i] == "Summary") {
#          X[[i]] <- as.character(X1[!is.na(X1), 1])
#       } else {
#          rownames(X1) <- X1[, 2]
#          X[[i]] <- X1[, -c(1:3)]
#       }
#    }
#    RODBC::odbcCloseAll()
#    class(X) <- "C2"
#    X
#}

#print.C2 <- function(x, ...) {
#   n <- grep(":", x$Summary)
#   cat(x$Summary[n], sep="\n")
#   cat("\nModel contains the following components:\n\n")
#   cat(names(x), sep="\n")
#}

#summary.C2 <- function(x, ...) {
#   cat(as.character(x$Summary), sep="\n")
#}

