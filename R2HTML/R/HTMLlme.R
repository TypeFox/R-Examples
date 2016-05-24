# $Id: HTMLlme.R 47 2008-05-23 17:29:31Z mentus $
# R2HTML.lme.r
# Formattted HTML printout of lme/nlme results
# Remove R2HTML:::  when in namespace
# Dieter Menne, (Dr. Menne Biomed Software Tübingen)
# dieter.menne@menne-biomed.de
#
# Revision 0.9, 15.2.2005 # first version forwarded to Eric

#% rm(list=ls(all=TRUE)) # test code, make sure we start fresh
#% library(nlme)
#% library(R2HTML)

# Assumes the following entry in R2HTML.CSS
#TD.cellinsideLeft {
#  padding: 5 5;
#  background: #FFFFFF;
#  text-align:left
#}
# In addition, for easier reading I use
# p {text-align:left}
# which overrides the default align=center.
#

# -------------------- HTML.summary.lme --------------------------------
# After print.summary.lme
HTML.summary.lme = function (x, file=HTMLGetFile(),
  digits = max(3,getOption("digits")-3),
  use.cormat = TRUE, # no symbolic.cor used when true
  symbolic.cor = p>4,
  signif.stars = getOption("show.signif.stars"),
  append=TRUE,verbose = FALSE, ...)
{
    verbose <- verbose || attr(x, "verbose")
    method <- x$method
    header = as.character()
    if (inherits(x, "nlme")) {
      header["Model"] <- "<font class='call'>Nonlinear</font>"
    }
    else {
      header["Model"] <- "<font class='call'>Linear</font>"
    }
    header["Method"]=ifelse(method == "REML", "REML","Max likelihood")
    header["Data"]= deparse(x$call$data)
    if (!is.null(x$call$subset)) {
      header["Subset"] <-
         deparse(asOneSidedFormula(x$call$subset)[[2]])
    }
    header["Observations"] <- x$dims[["N"]]
    f <- "Number of Groups: "
    dd <- x$dims
    Ngrps <- dd$ngrps[1:dd$Q]
    if ((lNgrps <- length(Ngrps)) == 1) {
      header["Groups"] <- Ngrps
    }
    else {
      sNgrps <- 1:lNgrps
      aux <- rep(names(Ngrps), sNgrps)
      aux <- split(aux, array(rep(sNgrps, lNgrps), c(lNgrps,
          lNgrps))[!lower.tri(diag(lNgrps))])
      names(Ngrps) <- unlist(lapply(aux, paste, collapse = " in "))
      header <- c(header,t(Ngrps)[1,])
    }
    if (verbose && !is.null(x$numItr)) {
      header["Iterations"] <- x$numIter
    }
    HTML(as.matrix(header),classcellinside="cellinsideLeft",file=file)
    HTML(data.frame(AIC = x$AIC, BIC = x$BIC, logLik = c(x$logLik),
        row.names = " "),file=file)
    # Random effects
    lapply(summary(x$modelStruct),HTML.reStruct, sigma = x$sigma,
       reEstimates = x$coef$random, verbose = verbose,file=file,
       digits=digits,use.cormat=use.cormat)
    # Fixed effects
    HTML(as.title("Fixed effects"),HR=3,file=file)
    fixF <- x$call$fixed
    if (inherits(fixF, "formula") || is.call(fixF)) {
      f <- deparse(x$call$fixed)
    }
    else {
      f <- deparse(lapply(fixF, function(el) as.name(deparse(el))))
    }
    HTMLli(paste("Formula:",f),file=file)
    xtTab <- as.data.frame(x$tTable)
    # Note: lme has a "DF" column, which should be printed as an int
    # According to Eric, this will be corrected for in a later
    # version, when individual formatting is possible
    # We force has.Pvalue to TRUE, because the column name is p-value, not Pr(..
    # If this should lead to problems, I could also rename the column
    HTML.coefmat(xtTab,signif.stars=TRUE,file=file,digits=digits,
      has.Pvalue=TRUE,...)
    if (nrow(x$tTable) > 1) {
      corr <- x$corFixed
      p <- NCOL(corr)
      if (p > 1) {
        HTMLli("Correlation of Coefficients:\n", file = file,append = TRUE, ...)
        # use of cormat has priority over symbolic
        if (use.cormat) {
            HTML.cormat(corr, file = file, ...)
        } else # conventional display
        if (symbolic.cor)
          HTML(symnum(corr)[-1, -p], file = file,...)
        else {
          HTML(corr[-1, -p, drop = FALSE], file = file, ...)
        }
      }
    }
    HTML(as.title("Standardized Within-Group Residuals"),HR=3,file=file)
    HTML(x$residuals,file=file)
}


# -------------------- HTML.reStruct ----------------------------------
# Modelled after print.reStruct in nlme
# This function used internally
HTML.reStruct = function (x, sigma = 1, reEstimates, verbose = FALSE,
  file=HTMLGetFile(),digits=max(3,getOption("digits")-3),use.cormat=FALSE,...)
{
  if (nlme::isInitialized(x)) {
    nobj <- length(x)
    if (is.null(namx <- names(x)))
       names(x) <- nobj:1
    aux <- t(array(rep(names(x), nobj), c(nobj, nobj)))
    aux[lower.tri(aux)] <- ""
    x[] <- rev(x)
    names(x) <- rev(apply(aux, 1, function(x) paste(x[x !=
            ""], collapse = " %in% ")))
    HTML(as.title("Random effects"),HR=3,file=file)
    for (i in seq(along = x)) {
      sm <-summary(x[[i]])
      p <- NCOL(sm)
      Level <- names(x)[i]
      resid  <-  p==length(x)
      HTML.summary.pdDiag(sm,sigma,Level=Level,resid=resid,file=file)
      if (p > 1) { # for p==1 we let HTML.summary.pdDiag do the job
        # use of cormat has priority over symbolic
        if (use.cormat) {
           HTML.cormat(sm, file = file, ...)
        } else # conventional display
            HTML.matrix(sm, file = file, ...)
      }
      if (verbose) {
        HTMLli("Random effects estimates",file=file)
        HTML(reEstimates[[i]],file=file,digits=digits)
      }
      HTML("\n",file=file)
     }
    }
    else {
      HTMLli("Uninitialized random effects structure",file=file)
    }
}

# -------------------- HTML.summary.pdDiag ---------------------------------
# Modelled after print.summary.pdDiag in nlme
HTML.summary.pdDiag <-
function (x, sigma = 1, rdig = 3, Level = NULL, resid = FALSE, ...)
{
  if (!(is.null(form <- attr(x, "formula")))) {
    f <- "Formula: "
    if (inherits(form, "formula")) {
      f <- paste(f,deparse(form))
      if (!is.null(Level)) {
        f <- paste(f," |", Level)
      }
    }
    else {
    if (length(form) == 1) {
      paste(f,deparse(form[[1]]))
      if (!is.null(Level)) {
        f <- paste(f," |", Level)
      }
    }
    else {
      f <- paste(f,deparse(lapply(form, function(el) as.name(deparse(el)))))
      f <- paste(f,"\n Level:", Level)
    }
  }
  HTMLli(f)
  }
  if (ncol(x) == 1) {
    if (resid) {
      HTML(array(sigma * c(attr(x, "stdDev"), 1),
        c(1, 2), list("StdDev:", c(names(attr(x, "stdDev")),
        "Residual"))), ...)
    }
    else {
      HTML(array(sigma * attr(x, "stdDev"), c(1, 1),
            list("StdDev:", names(attr(x, "stdDev")))),...)
    }
  }
  else {
    HTMLli(paste(" Structure: ", attr(x, "structName"),"\n", sep = ""))
    if (attr(x, "noCorrelation") | (1 >= (p <- dim(x)[2]))) {
      if (resid) {
        HTML(array(sigma * c(attr(x, "stdDev"), 1),
           c(1, p + 1), list("StdDev:", c(names(attr(x,
           "stdDev")), "Residual"))), ...)
      }
      else {
        HTML(array(sigma * attr(x, "stdDev"), c(1,p),
        list("StdDev:", names(attr(x, "stdDev")))),...)
      }
    }
    else {
      ll <- lower.tri(x)
      stdDev <- attr(x, "stdDev")
      x[ll] <- format(round(x[ll], digits = rdig),...)
      x[!ll] <- ""
      xx <- array("", dim(x), list(names(attr(x, "stdDev")),
             c("StdDev", "Corr", rep("", p - 2))))
      xx[, 1] <- format(sigma * attr(x, "stdDev"))
      xx[-1, -1] <- x[-1, -p]
      if (!is.null(colnames(x))) {
        xx[1, -1] <- abbreviate(colnames(x)[-p], minlength = rdig + 3)
      }
      if (resid) {
        x <- array("", dim(xx) + c(1, 0), list(c(rownames(xx),
             "Residual"), colnames(xx)))
        x[1:p, ] <- xx
        x[, 1] <- format(sigma * c(stdDev, 1))
        xx <- x
      }
      HTML(xx, ..., quote = FALSE)
    }
  }
  invisible(x)
}

