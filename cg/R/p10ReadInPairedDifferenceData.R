## $Id: $
## Paired Groups /  Difference Case

setClass("cgPairedDifferenceData",
         representation(dfr="data.frame",
                        dfru="data.frame",
                        dfr.gcfmt="data.frame",
                        settings="list"))

prepareCGPairedDifferenceData <- function(dfr, format="listed",
                                          analysisname="",
                                          endptname="",
                                          endptunits="",
                                          logscale=TRUE,
                                          zeroscore=NULL,
                                          addconstant=NULL,
                                          digits=NULL,
                                          expunitname="",
                                          refgrp=NULL,
                                          stamps=FALSE) {
  ##
  ## PURPOSE: Read-in the input data frame and convert it to
  ## "univariate style" if needed for the usual S modeling conventions.
  ##
  ## Input argument handling
  dataformat <-  validArgMatch(format, c("listed","groupcolumns"))
  validBoolean(logscale)
  validBoolean(stamps)
  validCharacter(expunitname)
  digits <- validArgDigits(digits)

  if(!is.data.frame(dfr)) {
    stop(cgMessage("The input data needs to be of class data.frame"))
  }
  if(!is.null(zeroscore) && !is.null(addconstant)) {
    stop(cgMessage("The zeroscore and addconstant arguments",
                   "cannot be specified at the same time.",
                   seeHelpFile("prepareCGPairedDifferenceData")))
  }
  if(logscale==FALSE && (!is.null(zeroscore) || !is.null(addconstant))) {
    stop(cgMessage("Both the zeroscore and addconstant arguments",
                   "must be NULL when logscale=FALSE.",
                   seeHelpFile("prepareCGPairedDifferenceData")))
  }

  if(expunitname=="" && dataformat=="groupcolumns" && ncol(dfr)==2) {
    expunitname <- "Experimental Unit"
  }
  else if((expunitname==""  && dataformat=="groupcolumns" && ncol(dfr)==3) || dataformat=="listed") {
    expunitname <- names(dfr)[1]
  }
  ## else keep specified expunitname argument

  ## End input argument handling

  if(dataformat=="groupcolumns") {
    validCGPairedDiffGroupColDfr(dfr)

    numcols <- ncol(dfr)
    numrows <- nrow(dfr)

    ## If the dataframe has 3 columns then the first column must be labels for the
    ## paired observation experimental units, such as subject identifiers.
    ## Code uses the term "expunit"
    
    ## If only 2 columns make it 3 and preserve a version of it
    ## that can later be
    ## used for calculating paired differences
    if(numcols==2) {
      dfr <- cbind(expunit=as.character(1:numrows), dfr)
    ##  expunitname <- "expunit"
    }
    else { ## numcols==3
      ## expunitname <- names(dfr)[1]
      names(dfr)[1] <- "expunit"
    }

    ## Preserve a version that can later be
    ## used for accessing paired differences
    dfr.gcfmt <- dfr ## "gc" acronymn for "groupcolumns"
    
    names(dfr.gcfmt)[2:3] <- c("grp1","grp2")
    
    ## Transformation to univariate form, called dfru
    dfrnames <- names(dfr)
    dfrendpt <- dfr[, 2:3]
    grpnames <- dfrnames[2:3]

    grp <- rep(grpnames, each=numrows)
    expunit <- rep(dfr[, 1], times=2)
    endpt <- unlist(dfrendpt)

    grpf <- factor(grp, levels=grpnames)
    expunitf <- factorInSeq(expunit)
    dfru <- data.frame(expunitf=expunitf, grpf=grpf,
                       endpt=I(endpt), row.names=NULL)
    dfru <- validCGPairedDiffListedDfr(dfru)

  }

  else if(dataformat=="listed") {
    dfru <- validCGPairedDiffListedDfr(dfr)
    
    ## Ensure ordering by group
    numcols <- ncol(dfr)
    numpairs <- nrow(dfr)/2

    ## Assume the number of columns equals 3
    grp <- dfr[, 2]

    grpnames <- unique(grp)
    grpf <- factor(grp, levels=grpnames)
    dfru <- dfr[order(grpf), ]

    expunit <- dfr[, 1]
    expunitf <- factorInSeq(expunit)
    dfru <- dfru[order(grpf, expunitf), ]
    row.names(dfru) <- NULL
    names(dfru) <- c("expunitf","grpf","endpt")

    ## Preserve a version that can later be
    ## used for calculating paired differences
    ## "gc" acronymn for "groupcolumns"
    dfr.gcfmt <- with(dfru,
                      data.frame(expunit=as.character(exptunitf[1:numpairs]),
                                 grp1=endpt[1:numpairs],
                                 grp2=endpt[-(1:numpairs)]))
  }
  
 
  digits <- as.integer(ifelse(is.null(digits),
                              getNumDigits(stripmiss(dfru[,3])), digits))

  if(!is.expression(endptname) && endptname=="") { endptname <- "Endpoint" }
  ## if(expunitname=="") { expunitname <- "Experimental Unit" }
  grpnames <- levels(dfru$grpf)

  ## Do we need to replace 0's with a non-zero score?
  if(logscale && !is.null(zeroscore)) {

    endpt <- dfru$endpt
    
    if(all(endpt > 0)) {
      stop(cgMessage("There do not appear to be any zero",
                     "values in the data, so zeroscore should",
                     "be set to NULL.",
                     seeHelpFile("prepareCGPairedDifferenceData")))
    }
    zeroscore <- validZeroScore(zeroscore, endpt)
    
    if(zeroscore >= min(endpt[endpt > 0])) {
      stop(cgMessage("The zeroscore must evaluate to a",
                     "value that is less than all the",
                     "the rest of the endpoint values."),
           seeHelpFile("prepareCGPairedDifferenceData"))
    }
    
    dfru$endpt[dfru$endpt==0] <- zeroscore
    dfr.gcfmt$grp1[dfr.gcfmt$grp1] <- zeroscore
    dfr.gcfmt$grp2[dfr.gcfmt$grp2] <- zeroscore
  }

  ## Or perhaps add a constant to all values so log is defined?
  if(logscale && !is.null(addconstant)) {
    if(!is.numeric(addconstant)) {
      addconstant <- validArgMatch(addconstant, "simple", "addconstant")
    }
    endpt <- dfru$endpt
    if(all(endpt > 0)) {
      stop(cgMessage("There do not appear to be any zero",
                     "values in the data, so addconstant should",
                     "be set to NULL.",
                     seeHelpFile("prepareCGPairedDifferenceData")))
    }
    addconstant <- validAddConstant(addconstant, endpt)
    dfru$endpt <- endpt + addconstant
    dfr.gcfmt$grp1 <- dfr.gcfmt$grp1 + addconstant  
    dfr.gcfmt$grp2 <- dfr.gcfmt$grp2 + addconstant  
  }

  ## If logscale, are all endpoint values greater than 0?
  if(logscale && any(dfru$endpt <= 0)) {
    stop(cgMessage("There seems to be at least one zero or",
                   "negative value in the endpoint data prepared",
                   "for analysis, so a",
                   "log scale analysis cannot be used, as all",
                   "values need to be greater than zero. ",
                   seeHelpFile("prepareCGPairedDifferenceData")))
  }


  ## Define reference group
  if(!is.null(refgrp)) {
    if(!is.element(refgrp, grpnames)) {
      stop(cgMessage("The refgrp argument does not appear",
                     "to match one of the group labels.",
                     seeHelpFile("prepareCGPairedDifferenceData")))
    }
    ## else keep specified refgrp value
  }
  else {
    refgrp <- grpnames[1]
  }

  ## Derive differences
  if(refgrp==grpnames[1]) {
    dfr.gcfmt$diffendpt <- with(dfr.gcfmt, grp2 - grp1)
  }
  else {
    dfr.gcfmt$diffendpt <- with(dfr.gcfmt, grp1 - grp2)
  }
  
  if(logscale) {
    if(refgrp==grpnames[1]) {
      dfr.gcfmt$difflogendpt <- with(dfr.gcfmt, log(grp2) - log(grp1))
    }
    else {
      dfr.gcfmt$difflogendpt <- with(dfr.gcfmt, log(grp1) - log(grp2))
    }
  }
  
  returnObj <- new("cgPairedDifferenceData",
                   dfr=dfr,
                   dfru=dfru,
                   dfr.gcfmt=dfr.gcfmt,
                   settings=list(
                     analysisname=analysisname,
                     endptname=endptname,
                     endptunits=endptunits,
                     endptscale=if(logscale) { "log" } else { "original" },
                     zeroscore=zeroscore,
                     addconstant=addconstant,
                     digits=digits,
                     grpnames=grpnames,
                     expunitname=expunitname,
                     refgrp=refgrp,
                     stamps=stamps))
  
  returnObj  
  
}


validCGPairedDiffGroupColDfr <- function(dfr) {
  ##
  ## PURPOSE:  Error checking on dfr in groupcolumn format, the input data frame
  ##
  numcols <- ncol(dfr)
  if(numcols < 2 || numcols > 3) {
    stop(cgMessage("The groupcolumn input data format needs to be 2 or 3",
                   "columns. If two columns, then",
                   "the first column needs to be the group identifier",
                   "and the second needs to be the endpoint.",
                   "If three columns, then",
                   "the first column needs to be the subject identifier",
                   "the second column needs to be the group identifier",
                   "and the third needs to be the endpoint.",
                   seeHelpFile("prepareCGPairedDiffData")))
  }

  ## No missing (NA) or non-numeric values in endpoint are allowed
  ## 2 Column Case
  if(numcols==2 && ( any(is.na((dfr[, 1:2]))) || !all(sapply(dfr, is.numeric)))) {
    stop(cgMessage("There seem to be missing or otherwise non-numeric",
                   "values in the endpoint data portion. These are",
                   "not allowed.",
                   seeHelpFile("prepareCGPairedDiffData")))
  }
  ## 3 Column Case
  if(numcols==3 && ( any(is.na(dfr[, 2:3])) || !all(sapply(dfr[, 2:3], is.numeric)))) {
    stop(cgMessage("There seem to be missing or otherwise non-numeric",
                   "values in the endpoint data portion. These are",
                   "not allowed.",
                   seeHelpFile("prepareCGPairedDiffData")))
  }
  ## Detect if the experimental unit identifier has duplicate values
  else if(numcols==3 && (length(unique(as.character(dfr[,1]))) !=  nrow(dfr))) {
    stop(cgMessage("The first column experimental unit identifier",
                   "needs to have all unique values.",
                   seeHelpFile("prepareCGPairedDiffData")))
  }
  
  else return(TRUE)
}


validCGPairedDiffListedDfr <- function(dfr) {
  ##
  ## PURPOSE:  Error checking on dfr in listed format.
  numcols <- ncol(dfr)
  if(numcols != 3) {
    stop(cgMessage("The listed input data format needs to be 3",
                   "columns.", 
                   "The first column needs to be the experimental unit identifier,",
                   "the second column needs to be the group identifier,",
                   "and the third is the endpoints.",
                   seeHelpFile("prepareCGPairedDiffData")))
  }
  
  ## 3 Column Case: Ensure the first and second columns are character
  ## and the first has all unique values, and the second has only 2 unique
  ## values
  if(numcols==3 &&
     ( (length(unique(as.character(dfr[, 1]))) != nrow(dfr)/2) ||
      (length(unique(as.character(dfr[, 2]))) != 2) )) {
    stop(cgMessage("The first column of the listed input data format",
                   "needs to have two sets of distinct values since",
                   "it is the experimental unit identifier.",
                   "The second column of the listed input data format",
                   "needs to have exactly 2 distinct values since",
                   "it is the group identifier.",
                   seeHelpFile("prepareCGPairedDiffData")))
  }

  ## In both 2 and 3 column cases the third column needs to be numeric
  if(any(!is.numeric(dfr[, numcols]))) {
    stop(cgMessage("The last column of the data frame needs to be",
                   "numeric since it is the endpoint",
                   seeHelpFile("prepareCGPairedDiffData")))
  }

  ## No missing values allowed
  if(any(is.na(dfr[, numcols]))) {
    stop(cgMessage("There seem to be missing",
                   "values in the endpoint data portion. These are",
                   "not allowed.",
                   seeHelpFile("prepareCGPairedDiffData")))
  }

  ## Number of rows needs to be even
  if(nrow(dfr) %% 2 != 0) {
    stop(cgMessage("There seem to be at least one experimental",
                   "unit where the pair of observations is not complete.",
                   seeHelpFile("prepareCGPairedDiffData")))
  }
  
  return(dfr)
}





