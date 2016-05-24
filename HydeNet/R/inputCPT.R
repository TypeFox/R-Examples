#' @rdname cpt
#' @export
#' @importFrom stats terms
#' 
inputCPT <- function(x, factorLevels, reduce=TRUE, ...) UseMethod("inputCPT")

#' @rdname cpt
#' @export
inputCPT.formula <- function(formula, factorLevels, reduce=TRUE, ...)
{
  variables       <- as.character(attr(stats::terms(formula), "variables"))[-1]
  dependentVar    <- variables[1]
  independentVars <- variables[-1]
  inputCPT_workhorse(variables, dependentVar, independentVars, 
                     factorLevels, reduce, ...)
}

#' @rdname cpt
#' @export

inputCPT.list <- function(x, factorLevels, reduce=TRUE, ...)
{
  Check <- ArgumentCheck::newArgCheck()
  
  if (!all(c("y","x") %in% names(x)))
  ArgumentCheck::addError(paste0("List object 'x' must contain character vectors ",
                                 "'y' and 'x'. See help('cpt')."),
                          Check)

  if (!all(unlist(lapply(x,is.character))))
  ArgumentCheck::addError(paste0("List object 'x' must contain character vectors ",
                                 "'y' and 'x'. See help('cpt')."),
                          Check)
  
  if (length(x[["y"]]) != 1)
  ArgumentCheck::addError(paste0("Element 'y' of list object 'x' must be a character ",
                                 "vector of length 1. See help('cpt')."),
                          Check)
  
  ArgumentCheck::finishArgCheck(Check)
  
  variables       <- c(x[["y"]], x[["x"]])
  dependentVar    <- x[["y"]]
  independentVars <- x[["x"]]
  
  inputCPT_workhorse(variables, dependentVar, independentVars, 
                     factorLevels, reduce, ...)
}


#******** UNEXPORTED FUNCTION

inputCPT_workhorse <- function(variables, dependentVar, independentVars,
                               factorLevels, reduce=TRUE, ...)
{
  hbar <- paste(paste(rep("-",80),collapse=""),"\n",sep="")
  factorEntryCommand <- function(variableName){
    cat(hbar, "Enter Factor Levels for node '", variableName,"':\n\n",
        "If this is a binary variable, enter '<yn>' as a shortcut.\n",
        "When finished, enter '<z>'.\n",
        "To repeat entry of the last inputted factor level, enter '<b>'.\n",
        "To start over entirely, enter '<s>'.\n",
        "To quit, enter <q>.", hbar, sep="")
  }
  
  if(missing(factorLevels)){  # solicit the names of factor levels from the console
    factorLevels <- vector(mode="list")
    for(i in seq_along(variables)){
      escapeFlag <- 0
      levelIndex <- 1
      tmp <- vector("character")
      factorEntryCommand(variables[i])
      while(!escapeFlag){
        IO <- readline(paste0("Level ",levelIndex," of '",variables[i],"':   "))
        if(IO == "<yn>"){
          if(levelIndex == 1) {
            tmp <- c("No","Yes")
            escapeFlag <- 1
          } else cat("(NOTE: <yn> only works when entering the FIRST factor level.)\n")
        } else if(IO == "<s>") {
          levelIndex <- 1
          tmp <- vector("character")
          factorEntryCommand(variables[i])
        } else if(IO == "<b>"){
          levelIndex <- max(c(1,levelIndex-1))
          if(levelIndex == 1) tmp <- vector("character") else tmp <- tmp[1:(levelIndex-1)]
        } else if(IO == "<z>"){
          escapeFlag <- 1
        } else if(IO == "<q>"){
          stop("User requested termination.")
        }
        else {
          tmp <- c(tmp, IO)
          levelIndex <- levelIndex + 1
        }
      }
      factorLevels[[variables[i]]] <- tmp
    }
  } else{
    if(!all(variables %in% names(factorLevels))){
      stop(paste("Variables",paste(variables,collapse=", "),
                 "not all in parameter 'factorLevels'."))
    }
    factorLevels <- factorLevels[variables]
    if(!all(unlist(lapply(factorLevels, is.character)))) {
      stop("Incompatible 'factorLevels' argument. See help('inputCPT()').")
    }
  }
  facValWidths <- unlist(lapply(factorLevels, function(x) max(nchar(x))))
  
  # input the conditional probabilities
  data <- expand.grid(factorLevels)
  
  
  if(reduce){
    cat(hbar,
        "NOTE: parameter 'reduce' is set to TRUE in inputCPT().\n",
        "      Conditional probabilities Pr(",dependentVar,"=",
        factorLevels[[dependentVar]][1]," | ", paste(independentVars,collapse=", "),
        ") will be calculated\n",
        "      as the complement of the inputted probabilities Pr(", dependentVar,
        " != ",factorLevels[[dependentVar]][1]," | ",
        paste(independentVars,collapse=", "), ").\n", hbar,sep="")
    data <- data[data[,dependentVar] %in% levels(data[,dependentVar])[-1],]
    cat("Enter the following conditional probabilities:\n")
  } else {
    cat(hbar, "Enter the following conditional probabilities, or positive\n",
        "numbers proportional to them (e.g., counts):\n")
  }
  cat("Use '<q>' to halt execution.\n",
      "To go back one step and re-enter, enter '<b>'.\n", hbar, sep="")
  
  formattedDepVarLvls <- format(as.character(data[,dependentVar]),
                                width = facValWidths[dependentVar])
  
  noNegativeProbs <- FALSE
  i <- 1
  optWarn <- options()$warn
  options(warn = -1)
  while(!noNegativeProbs){
    while(i <= nrow(data)){
      valid.IO <- FALSE;
      while(!valid.IO){
        formattedIndepVarLvls <- data[i, independentVars]
        formattedIndepVarLvls <- format(unlist(formattedIndepVarLvls),
                                        width=facValWidths[-1])
        prompt <- paste("Pr(",dependentVar,"=", formattedDepVarLvls[i], " | ",
                        paste(apply(cbind(names(data[i,independentVars]),
                                          formattedIndepVarLvls),
                                    1, paste, collapse="="),
                              collapse=", "),
                        "):   ", sep="")
        IO <- readline(prompt)
        if(IO == "<q>") stop("User requested termination.") else if(IO != "<b>"){
          IO.n <- as.numeric(IO)
          if(is.na(IO.n)) cat("Invalid numeric data entry. Try again:\n") else {
            if(reduce & (IO.n<0 | IO.n>1)){
              cat("Invalid probability given. Enter a number in [0,1]:\n")  
            } else if(IO.n<0){
              cat("Invalid count/probability given. Enter a non-negative number:\n")
            } else{
              valid.IO <- TRUE
              data[i,"wt"] <- IO.n
              i <- i + 1
            }
          }
        } else i <- max(i -1 , 1)
      }
    }
    options(warn = optWarn)
    
    if(reduce){
      # Add complement rows to the conditional probability data frame
      # if reduce=TRUE was used; check for errors involving sum of entered
      # conditional probabilities greater than 1
      complementProbs <- plyr::ddply(data, independentVars,
                                     function(data) c("wt"=1-sum(data$wt)))
      complementProbs[,dependentVar] <- levels(data[,dependentVar])[1]
      data <- rbind(data,complementProbs)
      if(min(data$wt)>=0) noNegativeProbs <- TRUE else{
        cat(hbar,"Invalid set of conditional probabilities given. There exists\n",
            "some combination of conditioning variables such that\n",
            "the sum of Pr(",dependentVar," != ",factorLevels[[dependentVar]][1]," | ",
            paste(independentVars,collapse=", "), ") is greater than 1.\n",
            "Please re-enter the conditional probabilities.\n",
            hbar, sep="")
      }
    } else noNegativeProbs <- TRUE
  } #end while(!noNegativeProbs) loop

  return(cpt(x = list(y = dependentVar, x = independentVars), 
             data = data, wt = data$wt))
} #end function inputCPT()

