###################################################
### chunk number 1:
###################################################


# 'algo.quality' calculates quality values
# like specifity, sensitivity for a surveillance method
#
# Parameters:
#      survResObj: object of class survRes, which includes the state chain and
#                               the computed alarm chain

######################################################################
## Hot fix function fixing two issues in the algo.quality function.
##
## Author: Michael Hoehle <http://www.math.su.se/~hoehle>
## Date:   2015-11-24
##
## 1) The function does not work if state or alarms are coded as TRUE/FALSE
##    instead of 0/1.
## 2) algo.quality doesn't work for sts objects.
##
## The function now branches on the appropriate thing to do depending on
## what class the argument is. This is not necessarily very good object
## oriented programming, but it works for now.
######################################################################

algo.quality <- function (sts, penalty = 20) {
  if (class(sts) == "survRes") {
    state <- sts$disProgObj$state[sts$control$range] * 1
    alarm <- sts$alarm * 1
  } else {
    if (class(sts) == "sts") {
      if (ncol(sts) > 1) { stop("Function only works for univariate objects.") }
      state <- sts@state*1
      alarm <- alarms(sts)*1
    } else {
      stop(paste0("Class ",class(sts)," not supported!"))
    }
  }

  state <- factor(state, levels = c(0, 1))
  alarm <- factor(alarm, levels = c(0, 1))
  confusionTable <- table(state, alarm)
  sens = confusionTable[2, 2]/(confusionTable[2, 2] + confusionTable[2,
      1])
  spec = confusionTable[1, 1]/(confusionTable[1, 2] + confusionTable[1,
      1])
  TP = confusionTable[2, 2]
  FN = confusionTable[2, 1]
  TN = confusionTable[1, 1]
  FP = confusionTable[1, 2]
  dist = sqrt(((1 - spec) - 0)^2 + (sens - 1)^2)
  if (!(is.element(1, state))) {
    lag = 0
  }
  else {
    lag <- c()
    outbegins <- c()
    varA <- which(state == 1)
    outbegins <- c(outbegins, varA[1])
    if (length(varA) > 1) {
      varB <- diff(varA)
      outbegins <- c(outbegins, varA[which(varB != 1) +
                                     1])
    }
    count <- 1
    for (i in outbegins) {
      if (count < length(outbegins)) {
        pos <- match(1, alarm[i:min(i + penalty, (outbegins[count +
                                                            1] - 1))])
        if (is.na(pos)) {
          lag <- c(lag, penalty)
        }
        else {
          lag <- c(lag, pos - 1)
        }
      }
      else {
        pos <- match(1, alarm[i:min(i + penalty, length(alarm))])
        if (is.na(pos)) {
          lag <- c(lag, penalty)
        }
        else {
          lag <- c(lag, pos - 1)
        }
      }
      count <- count + 1
    }
    lag <- mean(lag)
  }
  result <- list(TP = TP, FP = FP, TN = TN, FN = FN, sens = sens,
                 spec = spec, dist = dist, mlag = lag)
  class(result) <- "algoQV"
  return(result)
}



###################################################
### chunk number 2:
###################################################

print.algoQV <- function(x,...) {
  qualityValues <- c("TP", "FP", "TN", "FN", "Sens", "Spec", "dist", "mlag" )
  class(x) <- "list"
  result <- t(as.matrix(x))
  #Give the result matrix names
  dimnames(result)[[2]] <- qualityValues
  #Print to screen
  print(result)
  invisible()
}



###################################################
### chunk number 3:
###################################################
xtable.algoQV <- function(x, caption = NULL, label = NULL, align = NULL,
    digits = NULL, display = NULL, ...)  {
  n <- names(x)
  x <- matrix(x,nrow=1)
  dimnames(x)[[2]] <- n
  xtable(x,caption, label, align, digits, display, ...)
}


###################################################
### chunk number 4:
###################################################

# 'algo.call' calls the defined surveillance algorithms for
# a specified observed vector.
#
# Parameter
#       disProgObj: object of class survRes, which includes the state chain, the observed
#       control: specifies which surveillance systems should be used with their parameters.
#                The parameter funcName and range must be specified where funcName must be
#                the apropriate function (without 'algo.')
#       range (in control): positions in observed which should be computed


algo.call <- function(disProgObj, control = list( list(funcName = "rki1", range = range),
                                                   list(funcName = "rki", range = range, b = 2, w = 4, actY = TRUE),
                                                   list(funcName = "rki", range = range, b = 2, w = 5, actY = TRUE) ) ) {
  #Function to apply one algorithm to the disProgObj
  onecall <- function(i) {
    do.call(paste("algo.",control[[i]]$funcName, sep=""),
            list(disProgObj = disProgObj, control = control[[i]]))
  }

  #Apply each algorithm in the control list to the disProgObj
  survResults <- lapply(1:length(control),onecall)

  #Create some fancy naming..
  names(survResults) <- lapply(survResults,function(survObj) {survObj$control$name})

  #Done
  return(survResults)
}


###################################################
### chunk number 5:
###################################################

algo.compare <- function(survResList){
  return(t(sapply(survResList,algo.quality)))
}




###################################################
### chunk number 6:
###################################################

algo.summary <- function(compMatrices){

  # check if the input is large enough for summing
  if(length(compMatrices) < 1){
        stop("It's an empty list !")
  }
  if(length(compMatrices) == 1){
        return(compMatrices[[1]])
  }

 #Stupid conversion...
  compMatrices <- lapply(compMatrices,function(one) {
    n <- dimnames(one)
    one <- matrix(as.numeric(one),nrow=dim(one)[[1]])
    dimnames(one) <- n
    return(one)
  })

  # Compute the whole result
  wholeResult = compMatrices[[1]]
  lag = matrix(0,length(compMatrices),length(wholeResult[,1]))
  lag[1,] = wholeResult[,8]

  for(i in 2:length(compMatrices)){
    wholeResult = wholeResult + compMatrices[[i]]
    lag[i,] = compMatrices[[i]][,8]
  }

  # Sens (TP)
  wholeResult[,5] = wholeResult[,1]/(wholeResult[,1]+wholeResult[,4])
  # Spec (TN/(TN+FP))
  wholeResult[,6] = wholeResult[,3]/(wholeResult[,2]+wholeResult[,3])
  # dist
  wholeResult[,7] = sqrt((wholeResult[,6]-1)^2 + (wholeResult[,5]-1)^2)
  # median(lag)
  for(i in 1:length(wholeResult[,1])){
    wholeResult[i,8] = mean(lag[,i])
  }

  #class(wholeResult) <- "compMatrix" # comparison matrix
  return(wholeResult)
}




