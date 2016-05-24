setClass(Class="DFPdata",
    representation=representation(Subject="ANY", 
                                  Condition="ANY",  
                                  RT="numeric",  
                                  Correct="ANY",  
                                  Salience="matrix" ),
    )  

setValidity(Class="DFPdata", method=function(object) {

        if (!(length(object@Subject)==length(object@Condition) &
              length(object@Subject)==length(object@RT) &
              length(object@Subject)==length(object@Correct) &
              length(object@Subject)==dim(object@Salience)[1] )) {
            return("Input data not of equal length.")
        }

        if (!is.numeric(object@Salience)){return("Invalid data type for Salience." ) }
        if (any(object@RT<0)) { return("There are negative response times.") }
        if ( !is.logical(object@Correct) ) {
            if (! (object@Correct == 1 | object@Correct == 0 ) ) {
                return("Invalid levels for Correct.")
            }
        }
        return(TRUE)
    } )

setMethod("initialize", "DFPdata", function(.Object, ...) {
    value <- callNextMethod()
    validObject(value)
    value 
    } )

setMethod("show", "DFPdata", function(object) {
    cat("*** Double Factorial Paradigm Data ***\n")
    dat <- data.frame(object@Subject, object@Condition, object@RT, object@Correct)
    nchannels <- dim(object@Salience)[2]
    for ( i in 1:nchannels ) {
        dat <- cbind(dat, object@Salience[,i])
    }
    channels <- paste(rep("Channel", nchannels), 1:nchannels)
    names(dat) <- c("Subject", "Condition", "RT", "Correct", channels)
    show(dat)
  } )

setMethod("print", "DFPdata", function(x, ...) {
    cat("*** Double Factorial Paradigm Data ***\n")
    dat <- data.frame(x@Subject, x@Condition, x@RT, x@Correct)
    nchannels <- dim(x@Salience)[2]
    for ( i in 1:nchannels ) {
        dat <- cbind(dat, x@Salience[,i])
    }
    channels <- paste(rep("Channel", nchannels), 1:nchannels)
    names(dat) <- c("Subject", "Condition", "RT", "Correct", channels)
    print(dat)
  } )


setMethod("summary", "DFPdata", function(object) {
    cat("*** Double Factorial Paradigm Data ***\n")

    dat <- data.frame(object@Subject, object@Condition, object@RT, object@Correct)
    nchannels <- dim(object@Salience)[2]
    for ( i in 1:nchannels ) {
        dat <- cbind(dat, object@Salience[,i])
    }
    channels <- paste(rep("Channel", nchannels), 1:nchannels, sep="_")
    names(dat) <- c("Subject", "Condition", "RT", "Correct", channels)

    if( max(object@Salience) <= 2 ) { 
        salience <- c("-", "L", "H") 
    } else { salience <- sort(unique(object@salience)) }
    for ( cond in sort(unique(object@Condition)) ) {
      cat("Condition:  ", cond, "\n======================================\n")
      for ( subj in sort(unique(dat$Subject[dat$Condition==cond])) ) {
         cat("Subject:  ", subj, "\n--------------------------------------\n")
         sdat <- subset(dat, dat$Subject==subj & dat$Condition==cond)

         cat("\tMean RT\t\tAccuracy\n")
         for (l1 in sort(unique(sdat$Channel_1))) {
           for (l2 in sort(unique(sdat$Channel_2))) {
              salience2 <- paste(salience[l1+1],salience[l2+1], ":", sep="")
              cat(salience2, "\t")
              cat(signif(mean(sdat$RT[sdat$Channel_1==l1 & sdat$Channel_2==l2]),5),"  \t")
              cat(mean(sdat$Correct[sdat$Channel_1==l1 & sdat$Channel_2==l2] ))
                cat("\n")
            }
         }
        cat("\n")
      }
      cat("\n")
    }
  } )


#DFPdata <- function(x, ...) UseMethod("DFPdata")
#
#DFPdata.default <- function(Subject, Condition, RT, Correct, Channel1, Channel2, ...) {
#    if (!(length(Subject)==length(Condition) &
#          length(Subject)==length(RT) &
#          length(Subject)==length(Correct) &
#          length(Subject)==length(Channel1) &
#          length(Subject)==length(Channel2))) {stop("Input data not of equal length.")}
#
#
#    Subject <- as.factor(Subject)
#    Condition <- as.factor(Condition)
#
#    if ( !is.numeric(RT) ) { stop("Invalid response time data type (must be numeric).") }
#    if ( any(RT < 0) ) { stop("There are negative response times.") }
#    if ( !is.logical(Correct) ) {
#        if (! (Correct == 1 | Correct == 0 ) ) {
#            stop("Invalid levels for Correct.")
#        }
#        Correct <- as.logical(Correct)
#    }
#
#    if( any(Channel1%%1 > 0)) {stop("Invalid salience level on Channel 1") }
#    if( any(Channel2%%1 > 0)) {stop("Invalid salience level on Channel 2") }
#
#    Channel1 <- as.ordered(Channel1)
#    Channel2 <- as.ordered(Channel2)
#
#    dat <- data.frame(Subject=Subject, Condition=Condition, RT=RT,
#                        Correct=Correct, Channel1=Channel1, Channel2=Channel2) 
#    #class(dat) <- "DFPdata"
#    dat
#}
