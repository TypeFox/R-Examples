#' Probabilities Discrete Variables
#' 
#' Computes the probabilities of discrete variables.
#' 
#' @param stateNames A \code{"character"} array indicating the states of the variable.
#' @param Variable A \code{"numeric"} array containing the records of the variable.
#' @return A list of  \code{"numeric"} arrays:
#' \item{coeff}{Contains the probabilities.}
#' \item{sizeDataLeaf}{Number of records in each leaf of the discrete tree.}
#' @seealso \link{discreteVariablesStates}
#' @export
#' @examples
#' ## Discrete Variable
#' data <- data.frame(X=rep(c("yes", "no", "maybe"), 500))
#' data <- discreteVariables_as.character(data, "X")
#' n <- nrow(data)
#' 
#' ## Probabilities
#' s <- discreteVariablesStates(namevariables="X", discreteData=data)
#' states <- s[[1]]$states
#' p <- probDiscreteVariable(stateNames=states, Variable=data$X)
#' p
#'
probDiscreteVariable <- function(stateNames, Variable)
{
  coeff <- c(); sizeDataLeaf <- c()
  for(j in 1:length(stateNames)){
    sizeDataLeaf <- c(sizeDataLeaf, length(which(Variable==stateNames[j])))
    probability <- (length(which(Variable==stateNames[j]))+1)/(length(Variable)+length(stateNames))###Corrección de laplace
    coeff <- c(coeff, probability) 
  }
  return(list(coeff=coeff,sizeDataLeaf=sizeDataLeaf))
}

#' Goodness of discrete probabilities
#' 
#' Get the loglikelihood and the BIC for discrete models, i.e discrete Bayesian Networks.
#' 
#' @name goodnessDiscreteVariables
#' @rdname goodnessDiscreteVariables
#' @param discreteBN A list of multiples lists. Each list contains two entries,
#' the probabilities and the size of the data which is in each leaf of the discrete tree.
#' @param sameData A logical argument; if \code{FALSE} means differents datasets had been used for the learnings.
#' @return The loglikelihood and the BIC of the discrete network.
#' @examples
#' ## 1. EXAMPLE 
#' ## Discrete data
#' X <- rep(c("yes", "no", "maybe"), 500)
#' Y <- rep(c("M", "F"), 750)
#' data <- data.frame(X=X, Y=Y)
#' disVar <- c("X","Y")
#' data <- discreteVariables_as.character(data, discreteVariables=disVar)
#' n <- nrow(data)
#' 
#' ## Probabilities
#' s <- discreteVariablesStates(namevariables=disVar, discreteData=data)
#' p <- lapply(1:length(s), function(i) probDiscreteVariable(stateNames=
#' s[[i]]$states, Variable=data[,i]))
#'
#' ## Log-likelihood
#' getlogLikelihoodDiscreteBN(p)
#' 
#' ## BIC
#' getBICDiscreteBN(p, sameData = TRUE)
#' 
#' ## 2. EXAMPLE 
#' ## Discrete variables
#' X <- rep(c("1", "2", "3"), 500)
#' data <- data.frame(X=as.character(X))
#' s <- discreteVariablesStates(namevariables="X", discreteData=data)
#' p1 <- probDiscreteVariable(stateNames = s[[1]]$states, Variable = data[,1])
#' 
#' Y <- rep(c("YES", "NO"), 100)
#' data <- data.frame(Y = as.character(Y))
#' s <- discreteVariablesStates(namevariables = "Y", discreteData = data)
#' p2 <- probDiscreteVariable(stateNames = s[[1]]$states, Variable = data[,1])

#' ## Probabilities
#' P <- list(p1,p2)
#'
#' ## Log-likelihood
#' getlogLikelihoodDiscreteBN(P)
#' 
#' ## BIC
#' getBICDiscreteBN(P, sameData = TRUE)


#' @rdname goodnessDiscreteVariables
#' @export
getlogLikelihoodDiscreteBN <- function(discreteBN){
  loglike <- 0
  for(i in 1:length(discreteBN)){
    coeff <- discreteBN[[i]]$coeff[discreteBN[[i]]$coeff!=0]
    size <- discreteBN[[i]]$sizeDataLeaf[discreteBN[[i]]$sizeDataLeaf!=0]
    loglike <- loglike + sum(log(coeff)*(size))
  }
  return(loglike)
}


#' @rdname goodnessDiscreteVariables
#' @export
getBICDiscreteBN <- function (discreteBN, sameData = FALSE) 
{
  if(sameData) l <- length(discreteBN) else l <- 1
  sizeData <- sum(sapply(1:length(discreteBN), function(i) sum(discreteBN[[i]]$sizeDataLeaf)))/l
  nlevel <- sapply(1:length(discreteBN), function(i) length(discreteBN[[i]]$coeff))
  dimension <- sum(sapply(1:length(discreteBN), function(i) length(discreteBN[[i]]$coeff[!(discreteBN[[i]]$coeff%in%c(0,1))])*(1-1/nlevel[i])))
  bic <- getlogLikelihoodDiscreteBN(discreteBN) - (dimension*log(sizeData))/2
  return(bic)
}

#' Prints Discrete Learnings
#' 
#' Shows the results of univariate and conditional learning of a discrete BN.
#' 
#' @param BN A discrete learning.
#' @return The results are shown by screen.
#' @export
printDiscreteBN <- function(BN)
{
  cat("Potential(", BN$Child,")\n", sep="")
  if((length(BN$functions)<2)&&(length(BN$functions[[1]])==2)){
    cat(BN$functions[[1]]$coeff, "\n")
  } else{
    for(j in 1:length(BN$functions)){
      if(is.character(BN$functions[[j]]$interval)){
        cat("Parent:", BN$functions[[j]]$parent, "   \t Range =", paste("\"", BN$functions[[j]]$interval,"\"", sep=""),"\n")
        if(!is.null(BN$functions[[j]]$Px$coeff)) cat(BN$functions[[j]]$Px$coeff, "\n")
      } else {
        cat("Parent:", BN$functions[[j]]$parent, "   \t Range:", BN$functions[[j]]$interval[1], "<",BN$functions[[j]]$parent,"<", BN$functions[[j]]$interval[2], "\n")
        if(is.null(BN$functions[[j]]$Px[[1]])) next
        if(is.numeric(BN$functions[[j]]$Px[[1]])) cat(BN$functions[[j]]$Px[[1]], "\n")
        else for(i in 1:length(BN$functions[[j]]$Px)) cat(BN$functions[[j]]$Px[[i]]$coeff, "\n")
      }
    }
  }
  cat("\n")
}

