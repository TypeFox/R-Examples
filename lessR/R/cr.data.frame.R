.cr.data.frame <-
function(x, miss, show.n, digits.d,
         graphics, main, bottom, right, 
         pdf, pdf.width, pdf.height, ...)  {

  if (!is.null(digits.d) && digits.d<1) {
    cat("\n"); stop(call.=FALSE, "\n","------\n",      
        "\n>>> digits d is ", digits.d,  " and must be at least 1.\n")
  }

  if (miss == "pairwise") miss.type <- "pairwise.complete.obs"
  else if (miss == "listwise") miss.type <- "complete.obs"
  else if (miss == "everything") miss.type <- "everything"

  max.digits <- 0
  dig.warning <- FALSE

  # check for valid numeric variables and get largest decimal digits
  not.num <- integer(length=ncol(x))
  i.not <- 0
  for (i in 1:ncol(x)) {

    x.name <- names(x)[i]
    options(xname = x.name)

    #nu <- length(unique(na.omit(x[,i])))
    #if (.is.num.cat(x[,i], n.cat)) {
      #i.not <- i.not + 1
      #not.num[i.not] <- i
      #cat("\n")
      #cat("\n>>> Note:", x.name,  "is technically numeric, but only has ", 
          #nu, "<= n.cat =", n.cat, " levels,\n",
         #"      so treat as a categorical variable.\n",
         #"    To obtain the correlations decrease  n.cat  to specify a",
         #"lower number\n",
         #"      of unique values, such as with the function: theme\n",
         #"    Perhaps make this variable a factor with the R factor function.\n")
    #}

    if (!is.numeric(x[,i])) {
      i.not <- i.not + 1
      not.num[i.not] <- i
    }

    if (is.null(digits.d)) {
      dig.dec <- .max.dd(x[,i]) + 1 
      if (dig.dec < 2) dig.dec <- 2
      if (dig.dec > max.digits  &&  !dig.warning) max.digits <- dig.dec
      if (dig.dec > 10  &&  !dig.warning) {
        cat("\nThese data values contain ", dig.dec, " decimal digits. To enhance\n",
            "the readability of the output, only 4 decimal digits are\n",
            "displayed.  To customize this setting, use the digits.d  parameter.\n",
            "Example for Variables Y and X:  > cr(Y, by=X, digits.d=3)\n\n",
            sep="")
        dig.warning <- TRUE
        max.digits <- 4
      }
    }
  }

  # background
  tx <- character(length = 0)

  if (i.not > 0) {
    tx[length(tx)+1] <- paste("The following non-numeric variables are deleted",
        "from the analysis")
    for (i in 1:i.not) tx[length(tx)+1] <- paste(i, ". ", names(x)[not.num[i]],
      sep="")
    x <- x[, -not.num[1:i.not]]
    if (is.null(dim(x))) { 
      cat("\n"); stop(call.=FALSE, "\n","------\n",
        "A correlation matrix requires at least 2 variables.\n\n")
    }
  }
  
  if (!is.null(digits.d)) max.digits <- digits.d
  crs <- round(cor(x, use=miss.type, ...), digits=max.digits)  # cor matrix

  n.vars <- nrow(crs)
  tx[length(tx)+1] <- ""
  tx[length(tx)+1] <- paste("The correlation matrix contains",
    n.vars, "variables")

  txb <- tx


  # missing
  tx <- character(length = 0)

  tx[length(tx)+1] <- paste("Missing data deletion: ", miss)

  if (miss == "listwise") 
    tx[length(tx)+1] <- paste("   Sample size after deleted rows:",
      sum(complete.cases(x)), "\n")

  else if (miss == "pairwise") {
    n <- as.integer(crs)
    n <- matrix(n, nrow=n.vars)
    dimnames(n) <- dimnames(crs)
    for (i in 1:n.vars) {
      for (j in 1:n.vars) {
        n[i,j] <- sum(!is.na(x[i] - x[j]))  # non-missing values
      }
    }
    options(xname = "Missing Data Analysis")
    tot.miss <- sum(is.na(x))
    if (tot.miss == 0) 
      tx[length(tx)+1] <- paste("\n>>> No missing data")
    else {
      .ss.factor(as.vector(n), brief=TRUE)
      if (n.vars <= 15  ||  show.n) {
        txn <- .prntbl(n, 2, cc=" ")
        for (i in 1:length(txn)) tx[length(tx)+1] <- txn[i]
      }
      else
        tx[length(tx)+1] <- paste("To view the sample size for each ",
          "correlation, re-run the\n",
          "analysis with show.n=TRUE\n", sep="")
    }
  }  # end pairwise

  txm <- tx


  # cor matrix
  tx <- character(length = 0)

  if (n.vars <= 15) {
    tx[length(tx)+1] <- paste("Correlation Matrix")
    txcrs <- .prntbl(crs, 2, cc=" ")
    for (i in 1:length(txcrs)) tx[length(tx)+1] <- txcrs[i]
  }
  else {
    tx[length(tx)+1] <- "Variables in the correlation matrix:"
    tx[length(tx)+1] <- toString(names(x))
    tx[length(tx)+1] <- paste("\nTo view the correlation matrix, ",
     "enter the name of the returned object\n", 
     "followed by  $cors  such as  mycor$cors\n", sep="")
  }

  txc <- tx


  if (graphics  ||  pdf) {

    # set up graphics system for 2 windows
    #if (!pdf) {
      #.graphwin(2)
      #dev.set(which=3)
    #}
    #else { 
      #pdf.file <- "Cor_SPmatrix.pdf"
      #pdf(file=pdf.file, width=pdf.width, height=pdf.height)
    #}
  
    # scatter plot matrix
      #panel2.smooth <- function (x, y, pch=par("pch"), cex=.9,
        #col.pt=getOption("col.stroke.pt"), col.smooth=getOption("col.stroke.bar"),
        #span=2/3, iter=3, ...) 
      #{
          #points(x, y, pch=pch, col=col.pt, cex=cex)
          #ok <- is.finite(x) & is.finite(y)
          #if (any(ok)) 
            #lines(lowess(x[ok], y[ok], f=span, iter=iter), col=col.smooth, ...)
      #}

    #pairs(x, panel=panel2.smooth)

    #if (pdf) {  # terminate pdf graphics
      #dev.off()
      #.showfile(pdf.file, "scatter plot matrix")
    #}

    # heat map
    if (!pdf) {
      .graphwin(1)
      dev.set(which=3) 
    }
    else { 
      pdf.file <- "Cor_HeatMap.pdf"
      pdf(file=pdf.file, width=pdf.width, height=pdf.height)
    }

    if (is.null(main)) main <- "Correlations"

    for (i in 1:ncol(crs)) crs[i,i] <- 0
    cat("\nNote: To provide more color separation for off-diagonal\n",
        "      elements, the diagonal elements of the matrix for\n",
        "      computing the heat map are set to 0.\n", sep="")

    max.color <- getOption("col.heat")
    hmcols <- colorRampPalette(c("white",max.color))(256)

    heatmap(crs[1:ncol(crs),1:ncol(crs)], Rowv=NA, Colv="Rowv", symm=TRUE,
      col=hmcols, margins=c(bottom, right), main=main)

    if (pdf) {  # terminate pdf graphics
      dev.off()
      .showfile(pdf.file, "heat map")
      cat("\n\n")
    }
  }

  else
    if (is.null(options()$knitr.in.progress))
      cat("\n>>> To view comments, enter the name of the saved object, e.g., mycor\n")
      cat("\n>>> To view a heat map, set:  graphics=TRUE\n\n")

  return(list(txb=txb, txm=txm, txc=txc, cors=crs))

}
