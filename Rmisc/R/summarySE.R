#' Summarizes data
#' 
#' Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
#' 
#' @param data a data frame
#' @param measurevar the name of a column that contains the variable to be summariezed
#' @param groupvars a vector containing names of columns that contain grouping variables
#' @param na.rm a boolean that indicates whether to ignore NA's
#' @param conf.interval the percent range of the confidence interval (default is 95%)
#' @param .drop should combinations of variables that do not appear in the input data be preserved (FALSE) or dropped (TRUE, default)
#' 
#' @return a data frame with count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
#' 
#' @references http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)
#' 
#' @importFrom plyr ddply rename
#' 
#' @export
#' 
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This is does the summary; it's not easy to understand...
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun= function(xx, col, na.rm) {
                   c( N    = length2(xx[,col], na.rm=na.rm),
                      mean = mean   (xx[,col], na.rm=na.rm),
                      sd   = sd     (xx[,col], na.rm=na.rm)
                   )
                 },
                 measurevar,
                 na.rm
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean"=measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

#' Normalize within-group data
#' 
#' Norms the data within specified groups in a data frame; it normalizes each
#' subject (identified by idvar) so that they have the same mean, within each group
#' specified by betweenvars.
#' 
#' @param data a data frame.
#' @param idvar the name of a column that identifies each subject (or matched subjects)
#' @param measurevar the name of a column that contains the variable to be summariezed
#' @param betweenvars a vector containing names of columns that are between-subjects variables
#' @param na.rm a boolean that indicates whether to ignore NA's
#' @param .drop should combinations of variables that do not appear in the input data be preserved (FALSE) or dropped (TRUE, default)
#' 
#' @return a data frame with normalized data
#' 
#' @references http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)
#' 
#' @importFrom plyr ddply
#' 
#' @export
#' 
normDataWithin <- function(data=NULL, idvar, measurevar, betweenvars=NULL,
                           na.rm=FALSE, .drop=TRUE) {
  # Measure var on left, idvar + between vars on right of formula.
  data.subjMean <- ddply(data, c(idvar, betweenvars), .drop=.drop,
                         .fun = function(xx, col, na.rm) {
                           c(subjMean = mean(xx[,col], na.rm=na.rm))
                         },
                         measurevar,
                         na.rm
  )
  
  # Put the subject means with original data
  data <- merge(data, data.subjMean)
  
  # Get the normalized data in a new column
  measureNormedVar <- paste(measurevar, "Normed", sep="")
  data[,measureNormedVar] <- data[,measurevar] - data[,"subjMean"] +
    mean(data[,measurevar], na.rm=na.rm)
  
  # Remove this subject mean column
  data$subjMean <- NULL
  
  return(data)
}

#' Summarize within-subjects data
#' 
#' Summarizes data, handling within-subjects variables by removing inter-subject variability.
#' It will still work if there are no within-S variables. Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
#' If there are within-subject variables, calculate adjusted values using method from Morey (2008).
#' 
#' @param data a data frame
#' @param measurevar the name of a column that contains the variable to be summariezed
#' @param betweenvars a vector containing names of columns that are between-subjects variables
#' @param withinvars a vector containing names of columns that are within-subjects variables
#' @param idvar the name of a column that identifies each subject (or matched subjects)
#' @param na.rm a boolean that indicates whether to ignore NA's
#' @param conf.interval the percent range of the confidence interval (default is 95%)
#' @param .drop should combinations of variables that do not appear in the input data be preserved (FALSE) or dropped (TRUE, default)
#' 
#' @return a data frame with count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
#' 
#' @references http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)
#' 
#' @autoImports
#' 
#' @export
#' 
summarySEwithin <- function(data=NULL, measurevar, betweenvars=NULL, withinvars=NULL,
                            idvar=NULL, na.rm=FALSE, conf.interval=.95, .drop=TRUE) {
  
  # Ensure that the betweenvars and withinvars are factors
  factorvars <- sapply(data[, c(betweenvars, withinvars), drop=FALSE], FUN=is.factor)
  if (!all(factorvars)) {
    nonfactorvars <- names(factorvars)[!factorvars]
    message("Automatically converting the following non-factors to factors: ",
            paste(nonfactorvars, collapse = ", "))
    data[nonfactorvars] <- lapply(data[nonfactorvars], factor)
  }
  
  # Norm each subject's data    
  data <- normDataWithin(data, idvar, measurevar, betweenvars, na.rm, .drop=.drop)
  
  # This is the name of the new column
  measureNormedVar <- paste(measurevar, "Normed", sep="")
  
  # Replace the original data column with the normed one
  data[,measurevar] <- data[,measureNormedVar]
  
  # Collapse the normed data - now we can treat between and within vars the same
  datac <- summarySE(data, measurevar, groupvars=c(betweenvars, withinvars), na.rm=na.rm,
                     conf.interval=conf.interval, .drop=.drop)
  
  # Apply correction from Morey (2008) to the standard error and confidence interval
  #  Get the product of the number of conditions of within-S variables
  nWithinGroups    <- prod(sapply(datac[,withinvars, drop=FALSE], FUN=nlevels))
  correctionFactor <- sqrt( nWithinGroups / (nWithinGroups-1) )
  
  # Apply the correction factor
  datac$sd <- datac$sd * correctionFactor
  datac$se <- datac$se * correctionFactor
  datac$ci <- datac$ci * correctionFactor
  
  return(datac)
}