# Functions to process variables for binary classification

# the call to globalVariables removes 'Note' messages during R CMD CHECK about global variables in 
# columns used for ggplot
globalVariables(c('bin.names','bin.names.short','prob','bin.width','x','bin.count','iv'))

#
# public functions
#

BinaryPredictor <- function(iv, dv, ...){
  # generic method for binary predictors
  UseMethod("BinaryPredictor", iv)
}


BinaryPredictor.default <- function(iv, dv, ...){
  result <- DefaultResults(iv, dv) # get default values
  result$keep <- FALSE
  result$reason <- "Not numeric or factor"
  return(result)
}


BinaryPredictor.factor <- function(iv, dv, min.power=0.01, min.robustness=0.5, max.missing=0.99, max.levels=20, civ=NULL, copy.data=FALSE, name=NULL, ...){
  # execute basic checks -- all missing, unique values, etc.
  result <- BasicPredictorChecks(iv, civ, max.missing)
  # set the variable name to the deparsed iv, or use the name argument if supplied
  if (is.null(name)){
    result$name <- deparse(substitute(iv))
  } else {
    result$name <- name
  }
  
  if (!result$keep){
    # this variable failed basic checks, return
    return(result)
  }
  
  # check if number of levels exceed threshold
  if (length(levels(iv)) > max.levels){
    iv <- MergeLevels(iv, max.levels=max.levels)
  }
  
  # compute predictive power and WOE
  result$predictivePower <- PredictivePowerCv(iv, dv, ...)
  result$woe <- Woe(iv=iv, dv=dv, civ=civ)
  # remove the woe argument since it's a full-length vector of WOE values, keep only bin values in woe.levels
  result$woe$woe <- NULL
  
  # execute more checks using power and woe
  if (is.na(result$predictivePower$mean) | is.infinite(result$predictivePower$mean) | 
      is.nan(result$predictivePower$mean) | (result$predictivePower$mean < min.power)){
    result$keep <- FALSE
    result$reason <- "Low power"
    return(result)
  } else if (result$predictivePower$robustness < min.robustness) {
    result$keep <- FALSE
    result$reason <- "Low robustness"
    return(result)
  }
  # copy data to result object if the copy.data option is TRUE and if this is a continuous variable
  if (copy.data){
    if (!is.null(civ)){
      result$iv <- civ # copy the continuous data
      result$dv <- dv
    }
  }
  return(result)
}


BinaryPredictor.numeric <- function(iv, dv, min.power=0.01, min.robustness=0.5, max.missing=0.99, copy.data=FALSE, name=NULL, ...){
  # execute basic checks
  result.check <- BasicPredictorChecks(iv, iv, max.missing)
    # set the variable name to the deparsed iv, or use the name argument if supplied
  if (is.null(name)){
    result.check$name <- deparse(substitute(iv))
  } else {
    result.check$name <- name
  }
  
  if (!result.check$keep){
    # basic checks failed, return results
    return(result.check)
  }
  
  # extra arguments in ... may go to BinaryCut or PredictivePowerCv, sort them out
  args.BinaryCut <- names(formals(BinaryCut))
  args.PredictivePowerCv <- names(formals(PredictivePowerCv))
  dots <- list(...)
  
  # convert from numeric to factor and pass to BinaryPredictor.factor
  # use do.call to send only the relevant arguments to BinaryCut and BinaryPredictor.factor
  fiv <- do.call('BinaryCut', c( 
    list(iv=iv, dv=dv),
    dots[names(dots) %in% args.BinaryCut]))
  # include the continuous data civ so that linearity is computed in Woe
  result <- do.call('BinaryPredictor', c(
    list(iv=fiv, dv=dv, min.power=min.power, min.robustness=min.robustness, 
    max.missing=max.missing, civ=iv, 
    copy.data=copy.data, name=result.check$name), 
    dots[names(dots) %in% args.PredictivePowerCv]) )
  return(result)
}


BinaryPredictor.data.frame <- function(iv, dv, min.power=0.01, min.robustness=0.5, max.missing=0.99, verbose=FALSE, 
                                       copy.data=FALSE, ...){  
  # first argument is a data frame, loop over each element
  bplist <- list()
  nl <- 1
  if (verbose){cat("\nProcessing columns:\n")}
  # loop for each column in data frame
  for (col in names(iv)){
    if (verbose){cat("  ", col, "\n", sep="")}
    # get info for this predictor and save in list
    bplist[[nl]] <- BinaryPredictor(iv=iv[[col]], dv=dv, 
      min.power=min.power, min.robustness=min.robustness, max.missing=max.missing, copy.data=copy.data, ...)
    # set the name
    bplist[[nl]]$name <- col
    # increment counter
    nl <- nl + 1
  }
  # set names of list elements to match data frame columns
  names(bplist) <- names(iv)
  class(bplist) <- "BinaryPredictorList"
  attr(bplist, "nrow") <- nrow(iv)
  attr(bplist, "ntrue") <- table(dv)[2] # the second value is considered true / positive by convention
  return(bplist)
}


plot.BinaryPredictor <- function(x, y=NULL, type="bin", plot.missing=TRUE, ...){
  if (x$class == "factor"){
    p <- PlotFactor(x, ...)
  } else if ((x$class %in% c("numeric","integer")) & (type == "glm")){
    p <- PlotNumericGlm(x, plot.missing)
  } else if ((x$class %in% c("numeric","integer")) & (type == "bin")){
    p <- PlotNumericBin(x, plot.missing)
  } else {
    stop("Cannot print BinaryPredictor with class ", x$class)
  }
  return(p)
}


print.BinaryPredictorList <- function(x, file=NULL, silent=FALSE, ...){
  # this function generates summary info for a binary predictor list
  # generate a data frame summary of results
  # initialize vectors
  N <- length(x)
  keep <- rep(F,N)
  reason <- rep("x",N)
  name <- names(x)
  classname <- rep("x",N)
  missing <- rep(0,N)
  predictivePower <- rep(0,N)
  robustness <- rep(0,N)
  linearity <- rep(0,N)
  
  # loop to fill in values
  for (i in 1:N){
    keep[i] <- x[[i]]$keep
    reason[i] <- x[[i]]$reason
    classname[i] <- x[[i]]$class
    missing[i] <- x[[i]]$missing
    # if the predictive power was calculated then extract the values
    if (length(x[[i]]$predictivePower) > 1) {
      # power was calculated, fill in values
      predictivePower[i] <- x[[i]]$predictivePower$mean
      robustness[i] <- x[[i]]$predictivePower$robustness
      linearity[i] <- x[[i]]$woe$linearity
    } else {
      # power was not calculated, set to missing
      predictivePower[i] <- NA
      robustness[i] <- NA
      linearity[i] <- NA
    }
  }
  
  # create data frame with summary information
  df <- data.frame(
    name = name,
    predictivePower = predictivePower,
    robustness = robustness,
    linearity = linearity,
    missing = missing,
    class = classname,
    keep = keep,
    reason = reason
    )
  
  # data frame of valid variables sorted by power
  if (sum(df$keep) > 0){
    dfValid <- df[df$keep,]
    dfValid <- dfValid[ order(dfValid$predictivePower, decreasing=TRUE), ]
    dfValid$keep <- NULL
    dfValid$reason <- NULL
    dfValid$name <- PrependRowNumber(dfValid$name)
  } else {
    dfValid <- data.frame()
  }
  
  # data frame of invalid variables
  if (sum(!df$keep) > 0){
    dfInvalid <- df[!df$keep,]
    dfInvalid$keep <- NULL
    dfInvalid <- dfInvalid[ order(dfInvalid$predictivePower, dfInvalid$name, decreasing=TRUE), ]
    dfInvalid$name <- PrependRowNumber(dfInvalid$name)
  } else {
    dfInvalid <- data.frame()
  }
  
  # data frame of variables to linearize, rank by 1-linearity * power
  if (sum(df$keep & !is.na(df$linearity)) > 0) {
    dfLin <- df[ df$keep & !is.na(df$linearity), ]
    dfLin <- dfLin[ order( (1 - abs(dfLin$linearity)) * dfLin$predictivePower, decreasing=TRUE), ]
    dfLin$keep<- NULL
    dfLin$reason <- NULL
    dfLin$class <- NULL
    dfLin$missing <- NULL
    dfLin$name <- PrependRowNumber(dfLin$name)
  } else {
    dfLin <- data.frame()
  }
  
  # open an output file
  if (!is.null(file)) {
    # a file was provided, write output to the file
    sink(file=file)
  }
  
  # set width so wide records aren't wrapped
  old.width <- getOption("width")
  options(width=500)
  
  # write summary information
  nrow_ <- attr(x, "nrow")
  ntrue_ <- attr(x, "ntrue")
  nfalse_ <- nrow_ - ntrue_
  if (!silent){
    cat("\nRecords: ", nrow_)
    cat("\nDependent variable true rate:", format(ntrue_/nrow_, digits=4), " True/False values:",  ntrue_, "/", nfalse_, "\n")
    sp <- "  " # indentation spaces
    cat("\n\nValid variables ranked by predictive power:\n")
    print(dfValid, row.names=rep(sp, nrow(dfValid)), right=FALSE)
    cat("\n\nLinearization candidates ranked by (1-linearity) * predictivePower:\n")
    print(dfLin, row.names=rep(sp, nrow(dfLin)), right=FALSE)
    cat("\n\nInvalid variables:\n")
    print(dfInvalid, row.names=rep(sp, nrow(dfInvalid)), right=FALSE)
  }
  if (!is.null(file)) {
    # a file was provided, turn off sink
    sink()
  }
  
  # return width to original setting
  options(width=old.width)
  
  # return data frame of summary info
  return(df)
}


#
# private functions
#

PlotData <- function(obj){
  minwidth <- 0.03
  # generate dataframe for plot
  d <- data.frame(
    true.count = obj$woe$true.count,
    bin.count = obj$woe$bin.count,
    bin.names = as.factor(names(obj$woe$woe.levels)),
    logit = obj$woe$woe.levels )
  
  # set bin width
  d$bin.width <- d$bin.count / max(d$bin.count)
  
  # probability of true value within bin
  d$prob <- d$true.count / d$bin.count 
  
  # ensure bins are at least minimum width
  d$bin.width[ d$bin.width < minwidth ] <- minwidth
  # order bins by the names so bars are ordered, too
  levels(d$bin.names) <- levels(d$bin.names)[ order(levels(d$bin.names)) ]
  return(d)
}


PlotFactor <- function(obj, ...){
  # get a data frame representing bins
  df <- PlotData(obj)
  # plot a factor variable
  titleStr <- paste(obj$name, "\n", 
    "Power: ",  format(obj$predictivePower$mean, digits=3), 
    "  Robustness: ", format(obj$predictivePower$robustness, digits=3), "\n", 
    "  Bin count min: ", min(df$bin.count), 
    "  max: ", max(df$bin.count), sep="")
  df$bin.names.short <- ShortenStrings(df$bin.names, ...)
  p <- ggplot(df, aes(x=bin.names.short, y=prob, width=bin.width)) + geom_bar(stat="identity", position="identity")
  p <- p + coord_flip()
  p <- p + ylab("Probability dependent variable is true") + xlab(paste("Bins of", obj$name))
  p <- p + labs(title=titleStr)
  return(p)
}


PlotNumericGlm <- function(obj, plot.missing){
  if (class(obj$dv)=="factor"){
    dv <- as.numeric(obj$dv) - 1
  }
  df <- data.frame(iv=obj$iv, dv=dv)
  p <- ggplot(df, aes(x=iv, y=dv))
  p <- p + geom_point(position=position_jitter(height=0.1))
  p <- p + stat_smooth(method="glm", family="binomial")
  return(p)
}


PlotNumericBin <- function(obj, plot.missing){
  # get a data frame representing bins
  df <- PlotData(obj)
  # copy the data frame and exclude missing values if present in the last row
  if (obj$missing > 0){
    # missing values present
    idx <- 1 : nrow(df)-1
    missingBinCount <- obj$woe$bin.count[nrow(df)]
    missingLogit <- df$logit[nrow(df)]
  } else {
    # no missing values
    idx <- 1 : nrow(df)
    missingBinCount <- 0
    missingLogit <- NA
  }
  df <- df[idx, ] # remove row for missing bin if present
  
  # x values for plot are from mean of values in each bin
  df$x <- obj$woe$civ.bin.mean[idx]
    
  # create a string for plot title
  titleStr <- paste(obj$name, "\n", 
    "Power: ",  format(obj$predictivePower$mean, digits=3), 
    "  Robustness: ", format(obj$predictivePower$robustness, digits=3), 
    "  Linearity: ", format(obj$woe$linearity, digits=3), "\n", 
    "  Bin count min: ", min(df$bin.count), 
    "  max: ", max(df$bin.count), 
    "  Missing: ", missingBinCount, "/", format(obj$missing, digits=3), " Logit:", format(missingLogit, digits=4) )
  p <- ggplot(df, aes(x=x, y=logit))
  p <- p + ylab("Logit of dependent variable")
  # add a line representing missing values
  if (obj$missing>0 & plot.missing){
    p <- p + geom_hline(yintercept=missingLogit, color="red", linetype=2, alpha=0.7)
  }
  # plot points, edit titles
  p <- p + geom_point(aes(size=bin.count)) + labs(title=titleStr) + xlab(paste("Bins of", obj$name))
  # add a line computed with regression weighted by bin counts, if there are only 2 points then don't smooth
  if (length(df$bin.count)>2){
    p <- p + stat_smooth(method="lm", aes(weight=bin.count))
  } else {
    p <- p + geom_line(color="blue")
  }
  return(p)
}


DefaultResults <- function(iv, civ){
  # set a list of default values that will be modified in later steps
  result <- list(
    name = "", # TODO: get unevaluated name
    keep = TRUE,
    reason = "",
    missing = sum(is.na(iv)) / length(iv), # fraction of missing values
    class = class(iv)[1], # use only the first element, class() can return two elements for POSIXct, POSIXt
    predictivePower = NA,
    woe = NA
  )
  # handle the case where a continuous variable is provided, use its class rather than iv
  if (!is.null(civ)){
    result$class <- class(civ)[1]
    result$missing <- sum(is.na(civ)) / length(civ)
  }
  class(result) <- "BinaryPredictor"
  return(result)
}


BasicPredictorChecks <- function(iv, civ, max.missing){
  # execute basic checks on the variable
  result <- DefaultResults(iv, civ)
  
  # basic checks
  if (length(unique(iv)) == 1){
    result$keep <- FALSE
    result$reason <- "Single value"
    return(result)
  } else if (result$missing > max.missing){
    result$keep <- FALSE
    result$reason <- "Missing values"
    return(result)
  }
  return(result)
}


PrependRowNumber <- function(s){
  # given a vector of strings "a","b", appends row number " 1 a", " 2 b"
  nums <- sprintf("%-5d", 1:length(s))
  return(paste(nums, s, sep=""))
}
