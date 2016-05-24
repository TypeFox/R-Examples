### write.fwf.R
###------------------------------------------------------------------------
### What: Write fixed width format - code
### $Id: write.fwf.R 1959 2015-04-25 07:56:23Z warnes $
### Time-stamp: <2008-08-05 12:11:27 ggorjan>
###------------------------------------------------------------------------

write.fwf <- function(x,
                      file="",
                      append=FALSE,
                      quote=FALSE,
                      sep=" ",
                      na="",
                      rownames=FALSE,
                      colnames=TRUE,
                      rowCol=NULL,
                      justify="left",
                      formatInfo=FALSE,
                      quoteInfo=TRUE,
                      width=NULL,
                      eol="\n",
                      qmethod=c("escape", "double"),
                      scientific=TRUE,
                      ...)
{
  ## --- Setup ---

  dapply <- function(x, FUN, ..., simplify=TRUE)
      {
          if(is.data.frame(x))
              return(sapply(x, FUN, ..., simplify=simplify))
          else if(is.matrix(x))
              return(apply(x, 2, FUN, ...))
          else
              stop("x must be a data.frame or a matrix")
      }

  if(!(is.data.frame(x) || is.matrix(x)))
    stop("'x' must be a data.frame or matrix")
  if(length(na) > 1)
    stop("only single value can be defined for 'na'")

  if(!scientific)
      {
          option.scipen <- getOption("scipen")
          on.exit( function() options("scipen"=option.scipen) )
          options("scipen"=100)
      }


  if(rownames) {
    x <- as.data.frame(x)
    x <- cbind(rownames(x), x)
    rowColVal <- ifelse(!is.null(rowCol), rowCol, "row")
    colnames(x)[1] <- rowColVal
  }
  colnamesMy <- colnames(x)
  if(length(colnamesMy)==0)
      colnamesMy <- paste( "V", 1:ncol(x), sep="")

  nRow <- nrow(x)
  nCol <- length(colnamesMy)

  widthNULL <- is.null(width)
  if(!widthNULL && length(width) != nCol) {
    warning("recycling 'width'")
    widthOld <- width
    width <- integer(length=nCol)
    width[] <- widthOld
  }

  ## --- Format info ---

  retFormat <- data.frame(colname=colnamesMy,
                          nlevels=0,
                          position=0,
                          width=0,
                          digits=0,
                          exp=0,
                          stringsAsFactors=FALSE)

  ## Which columns are numeric like
  isNum <- dapply(x, is.numeric)
  ## is.numeric picks also Date and POSIXt
  isNum <- isNum & !(dapply(x, inherits, what="Date") |
                     dapply(x, inherits, what="POSIXt"))

  ## Which columns are factors --> convert them to character
  isFac <- dapply(x, is.factor)
  if(any(isFac))
      ## This conditional is necessary because if x is a matrix, even if
      ## all(isFAC==FALSE), this assignment will coerce it to mode
      ## character.  This isn't a problem for dataframes.
      x[, isFac] <- sapply(x[, isFac, drop=FALSE], as.character)

  ## Collect information about how format() will format columns.
  ## We need to get this info now, since format will turn all columns to character
  tmp <- dapply(x, format.info, ..., simplify=FALSE)
  if(is.matrix(x)) tmp <- as.data.frame(tmp)
  tmp1 <- sapply(tmp, length)
  tmp <- t(as.data.frame(tmp))
  retFormat$width <- tmp[, 1]
  ## Collect other details for numeric columns
  if(any(isNum)) {
    ## Numeric columns with digits
    test <- tmp1 > 1
    if(any(test)) {
      retFormat[test, c("digits", "exp")] <- tmp[test, c(2, 3)]
      ## Numeric columns with scientific notation
      test2 <- tmp[test, 3] > 0
      if(any(test2)) ## adding +1; see ?format.info
        retFormat[test, ][test2, "exp"] <- retFormat[test, ][test2, "exp"] + 1
    }
  }

  ## --- Format ---

  ## store original object in 'y'
  y <- x

  ## Formatting (to character)
  for(i in 1:nCol) {
    if(widthNULL) {
      tmp <- NULL
    } else {
      tmp <- width[i]
    }
    ## Due to na.encode bug in format() in 2.7.1; na.encode=TRUE should
    ## return NA values and not "NA", but even then we rely on the
    ## following test to "fiddle" with the value in 'na' argument since -
    ## NA should not increase the width of column with width 1, while wider
    ## value for 'na' should increase the width
    test <- is.na(y[, i])
    ## Make a copy to make sure we get character after first format() - Date class caused problems
    x2 <- character(length=nRow)
    ## Add formatted values
    x2[!test] <- format(y[!test, i], justify=justify, width=tmp, ...)
    ## Add 'na' value
    x2[test] <- na
    ## Replace the original
    x[, i] <- x2
    ## Collect width (again)
    tmp2 <- format.info(x2, ...)[1]
    ## Reformat if 'na' value change the width of the column
    if(tmp2 != retFormat[i, "width"]) {
      retFormat[i, "width"] <- tmp2
      ## ifelse() makes sure that numeric columns are justified to right
      x[, i] <- format(x[, i], justify=ifelse(isNum[i], "right", justify),
                       width=tmp, ...)
    }
    ## Reformat 'na' value if it is narrower than the width of the column
    if(nchar(na) < retFormat[i, "width"]) {
      x[test, i] <- format(na, justify=ifelse(isNum[i], "right", justify),
                           width=retFormat[i, "width"], ...)
    }
  }

  ## Number of levels for "non-numeric"" columns
  if(any(!isNum)) {
    retFormat[!isNum, "nlevels"] <- dapply(x[, !isNum, drop=FALSE],
                                           function(z) length(unique(z)))
  }

  ## Check that width was not to small
  if(!widthNULL) {
    test <- retFormat$width > width
    if(any(test)) {
      tmpCol <- paste(colnamesMy[test], collapse=", ")
      tmpWidth <- paste(width[test], collapse=", ")
      tmpNeed <- paste(retFormat$width[test], collapse=", ")
      stop(paste("'width' (", tmpWidth, ") was too small for columns: ",
                 tmpCol, "\n 'width' should be at least (", tmpNeed, ")",
                 sep=""))
    }
  }

  ## --- Write ---

  if(colnames) {
    if(rownames && is.null(rowCol)) colnamesMy <- colnamesMy[-1]
    write.table(t(as.matrix(colnamesMy)),
                file=file,
                append=append,
                quote=quote,
                sep=sep,
                eol=eol,
                na=na,
                row.names=FALSE,
                col.names=FALSE,
                qmethod=qmethod)
 }

  write.table(x=x,
              file=file,
              append=(colnames || append),
              quote=quote,
              sep=sep,
              eol=eol,
              na=na,
              row.names=FALSE,
              col.names=FALSE,
              qmethod=qmethod)

  ## --- Return format and fixed width information ---

  if(formatInfo) {
    ## be carefull with these ifelse constructs
    retFormat$position[1] <- ifelse(quote, ifelse(quoteInfo, 1, 2), 1)
    if(ifelse(quote, quoteInfo, FALSE)) retFormat$width <- retFormat$width + 2
    N <- nrow(retFormat)
    if(N > 1) {
      for(i in 2:N) {
        retFormat$position[i] <- retFormat$position[i - 1] +
          retFormat$width[i - 1] + nchar(x=sep, type="chars") +
            ifelse(quote, ifelse(quoteInfo, 0, 1), 0)
      }
    }
    if(rownames && is.null(rowCol)) {
      retFormat <- retFormat[-1,]
      rownames(retFormat) <- 1:(N-1)
    }
    return(retFormat)
  }
}

###------------------------------------------------------------------------
### write.fwf.R ends here
