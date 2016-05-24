#' Learning MoTBFs in a Network
#'
#' Learn mixtures of truncated basis functions in a full hybrid network.
#' 
#' @param graph A network of the class \code{"bn"}, \code{"graphNEL"} or \code{"network"}.
#' @param data A datase of class \code{"data.frame"}; it can contain continuous and discrete variables.
#' @param numIntervals A \code{"numeric"} value indicating the maximum number of intervals in which we 
#' want to split the domain of the parent variables.
#' @param POTENTIAL_TYPE A \code{"character"} string specifying the posibles potential
#' types, must be one of \code{"MOP"} or \code{"MTE"}.
#' @param maxParam A \code{"numeric"} value which indicate the maximum number of coefficients in the function. By default it is \code{NULL}; 
#' if not, the output is the function which gets the best BIC with at most this number of parameters.
#' @param s A \code{"numeric"} coefficient which fixes the confidence of the prior knowledge 
#' we are going to introduce. By default it is \code{NULL}, only we must modify it if we want 
#' to incorporate prior information to the fits.
#' @param priorData Prior dataset with values of the variables we have information apriori about.
#' This dataset must be of \code{"data.frame"} class.
#' @return A list of lists. Each list contains two elements
#' \item{Child}{A \code{"charater"} string which contains the name of the child variable.}
#' \item{functions}{A list with three elements: the name of the parents, a \code{"numeric"} vector
#' with the limits of the interval and the fitted function in this interval.}
#' @details If the variable is discrete then it computes the probabilities and the size of each leaf.
#' @seealso \link{printBN} and \link{ecoli}
#' @export
#' @examples
#' 
#' ## Dataset Ecoli
#' require(MoTBFs)
#' data(ecoli)
#' data <- ecoli[,-c(1)] ## remove variable sequence
#' 
#' ## Directed acyclic graph
#' dag <- LearningHC(data)
#' 
#' ## Learning BN
#' intervals <- 3
#' potential <- "MOP"
#' P1 <- MoTBFs_Learning(graph = dag, data = data, POTENTIAL_TYPE=potential,
#' maxParam = 5)
#' printBN(P1)
#' 
#'  ## Learning BN
#' intervals <- 4
#' potential <- "MTE"
#' P2 <- MoTBFs_Learning(graph = dag, data = data, numIntervals = intervals, POTENTIAL_TYPE=potential,
#' maxParam = 15)
#' printBN(P2)
#' 
#' 
MoTBFs_Learning <- function(graph, data, numIntervals, POTENTIAL_TYPE, maxParam=NULL,s=NULL, priorData=NULL)
{
  childrenAndParents <- getChildParentsFromGraph(graph,colnames(data))
  MoTBFs <- c()
  for(i in 1:length(childrenAndParents)){
    if(length(childrenAndParents[[i]])==1){
      
      ## Single variable
      Child <- childrenAndParents[[i]]
      if(Child%in%colnames(priorData)) priorParent <- priorData[,Child]
      else priorParent <- NULL
      
      if(is.numeric(data[,Child])){
        
        ## Numeric child
        if(is.null(priorParent)) PX <- univMoTBF(data[,Child], POTENTIAL_TYPE)
        else PX <- learnMoTBFpriorInformation(priorParent, data[,Child], s, POTENTIAL_TYPE, maxParam=maxParam)$posteriorFunction 
        UnivMoTBFs <- list(PX)
        information <- list(Child=Child, functions=UnivMoTBFs)
        MoTBFs[[length(MoTBFs)+1]] <- information
      }else{
        
        ##Discrete child
        states <- discreteVariablesStates(Child, data)[[1]]$states
        prob <- probDiscreteVariable(states, data[,Child])
        UnivDisc <- list(prob)
        information <- list(Child=Child, functions=UnivDisc)
        MoTBFs[[length(MoTBFs)+1]] <- information
      }
    } else{
      
      ## Conditional variables
      Child <- childrenAndParents[[i]][1]
      Parents <- childrenAndParents[[i]][2:length(childrenAndParents[[i]])]
      if((Child%in%colnames(priorData))&&(all(Parents%in%colnames(priorData)))){
        priorChild <- priorData[,childrenAndParents[[i]]]
        colnames(priorChild) <- childrenAndParents[[i]]
      } else if (Child%in%colnames(priorData)) {
        priorChild <- as.data.frame(priorData[,Child])
        colnames(priorChild) <- Child
      } else{ 
        priorChild <- NULL
      }
      
      result <- conditionalMethod(data, Parents, Child, numIntervals, POTENTIAL_TYPE, maxParam=maxParam, s, priorChild)
      information <- list(Child=Child, functions=result)
      MoTBFs[[length(MoTBFs)+1]] <- information
    }
  }
  return(MoTBFs) 
}

#' Prints BN Results
#' 
#' Prints the results of a hybrid Bayesian network
#' 
#' @param MoTBF.BN The output of the method \code{MoTBFs_Learning()}.
#' @return The results of the fitted functions in the full network.
#' @seealso \link{MoTBFs_Learning}
#' @export
#' @examples
#' 
#' ## Dataset Ecoli
#' require(MoTBFs)
#' data(ecoli)
#' data <- ecoli[,-c(1)] ## remove variable sequence
#' 
#' ## Directed acyclic graph
#' dag <- LearningHC(data)
#' 
#' ## Learning BN
#' intervals <- 3
#' potential <- "MOP"
#' P <- MoTBFs_Learning(graph = dag, data = data, POTENTIAL_TYPE=potential,
#' maxParam = 15)
#' printBN(P)
#' 
printBN <- function(MoTBF.BN)
{
  for(i in 1:length(MoTBF.BN)){
    if(is.motbf(MoTBF.BN[[i]]$functions[[length(MoTBF.BN[[i]]$functions)]])||is.motbf(MoTBF.BN[[i]]$functions[[length(MoTBF.BN[[i]]$functions)]][[length(MoTBF.BN[[i]]$functions[[length(MoTBF.BN[[i]]$functions)]])]])){
      if((length(MoTBF.BN[[i]]$functions)==1)&&(is.motbf(MoTBF.BN[[i]]$functions[[1]]))){
        cat("Potential(", MoTBF.BN[[i]]$Child, ")\n", sep="")
        print(              MoTBF.BN[[i]]$functions[[1]])
        cat("\n")
      } else{
        cat("Potential(", MoTBF.BN[[i]]$Child,")\n", sep="")
        printConditional(MoTBF.BN[[i]]$functions)
        cat("\n")
      } 
    } else {
      printDiscreteBN(MoTBF.BN[[i]])
      cat("\n")
    }
  }
}

#' BIC of an MoTBF BN
#' 
#' Compute the Bic score and the loglikelihood from the fitted MoTBFs functions 
#' of a hybrid Bayesian network.
#' 
#' @name goodnessMoTBFBN
#' @rdname goodnessMoTBFBN
#' @param MoTBF.BN The output of the 'MoTBF_Learning' method.
#' @param data The dataset of class "data frame".
#' @return A numeric value giving the log-likelihood of the BN.
#' @seealso \link{MoTBFs_Learning}
#' @examples
#' 
#' ## Dataset Ecoli
#' require(MoTBFs)
#' data(ecoli)
#' data <- ecoli[,-c(1)] ## remove variable sequence
#' 
#' ## Directed acyclic graph
#' dag <- LearningHC(data)
#' 
#' ## Learning BN
#' intervals <- 3
#' potential <- "MOP"
#' P1 <- MoTBFs_Learning(graph = dag, data = data, POTENTIAL_TYPE=potential,
#' maxParam = 5)
#' logLikelihood.MoTBFBN(P1, data) ##BIC$LogLikelihood
#' BIC <- BiC.MoTBFBN(P1, data)
#' BIC$BIC
#'
#' ## Learning BN
#' intervals <- 2
#' potential <- "MTE"
#' P2 <- MoTBFs_Learning(graph = dag, data = data, POTENTIAL_TYPE=potential,
#' maxParam = 10)
#' logLikelihood.MoTBFBN(P2, data) ##BIC$LogLikelihood
#' BIC <- BiC.MoTBFBN(P2, data)
#' BIC$BIC 


#' @rdname goodnessMoTBFBN
#' @export
logLikelihood.MoTBFBN <- function(MoTBF.BN, data)
{
  originalData <- data; valuesT <- c()
  for(i in 1:length(MoTBF.BN)){
    
    ## Single Variable
    if(((length(MoTBF.BN[[i]]$functions)==1)&&(is.motbf(MoTBF.BN[[i]]$functions[[1]])))||((length(MoTBF.BN[[i]]$functions)==1)&&(length(MoTBF.BN[[i]]$functions[[1]])==2))){
      data <- originalData; Y <- data[, MoTBF.BN[[i]]$Child]
      if(is.numeric(Y)){
        
        ## Continuous child
        values <- as.function(MoTBF.BN[[i]]$functions[[1]])(Y)
        valuesT <- c(valuesT,values)
      } else{
        
        ## Discrete Child
        st <- sum(MoTBF.BN[[i]]$functions[[1]]$sizeDataLeaf)
        for(j in 1:length(MoTBF.BN[[i]]$functions[[1]]$coeff)){
          values=rep(MoTBF.BN[[i]]$functions[[1]]$coeff[j],MoTBF.BN[[i]]$functions[[1]]$sizeDataLeaf[j]*length(Y)/st)
          valuesT=c(valuesT,values)
        }
      }
    } else{
      
      ##Conditional Variables
      N <- length(valuesT); Parents <- c()
      for(t in 1:length(MoTBF.BN[[i]]$functions)){
        parent <- MoTBF.BN[[i]]$functions[[t]]$parent
        Parents <- c(Parents, parent)
      }
      fullData <- c(); enter <- c(); data <- originalData
      disc <- names(which(sapply(data, is.character)==T))
      if(length(disc)!=0) states <- discreteVariablesStates(disc, data)
      for(j in 1:length(MoTBF.BN[[i]]$functions)){
        if(length(MoTBF.BN[[i]]$functions[[j]])==2){
          if(is.numeric(data[, MoTBF.BN[[i]]$functions[[j]]$parent])) data <- splitdata(data, MoTBF.BN[[i]]$functions[[j]]$parent , MoTBF.BN[[i]]$functions[[j]]$interval[1], MoTBF.BN[[i]]$functions[[j]]$interval[2])
          else data <- subset(data,(data[,MoTBF.BN[[i]]$functions[[j]]$parent]==MoTBF.BN[[i]]$functions[[j]]$interval))
          
          fullData[[length(fullData)+1]] <- data
          enter <- c(enter, "NO")
        } else{
          if(is.numeric(data[, MoTBF.BN[[i]]$functions[[j]]$parent])) data <- splitdata(data, MoTBF.BN[[i]]$functions[[j]]$parent , MoTBF.BN[[i]]$functions[[j]]$interval[1], round(MoTBF.BN[[i]]$functions[[j]]$interval[2], digits=2)+2^-1)
          else data <- subset(data,(data[,MoTBF.BN[[i]]$functions[[j]]$parent]==MoTBF.BN[[i]]$functions[[j]]$interval))         
          
          fullData[[length(fullData)+1]] <- data
          enter <- c(enter, "YES")
          
          if(is.character(data[,MoTBF.BN[[i]]$Child])){
            prob <- MoTBF.BN[[i]]$functions[[j]]$Px
            if(is.null(prob$coeff)) prob <- prob[[1]]
            values <- c(); st <- sum(prob$sizeDataLeaf)
            for(m in 1:length(prob$coeff)){
              val <- rep(prob$coeff[m], prob$sizeDataLeaf[m]*length(data[,MoTBF.BN[[i]]$Child])/st)
              values <- c(values, val)
            }
          } else {
            Y <- data[,MoTBF.BN[[i]]$Child]
            if(!is.motbf(MoTBF.BN[[i]]$functions[[j]]$Px)){
              values <- c()
            } else{
              values <- as.function(MoTBF.BN[[i]]$functions[[j]]$Px)(Y)
            }
          }
          
          valuesT <- c(valuesT,values)
          if(j!=length(MoTBF.BN[[i]]$functions)){
            if(MoTBF.BN[[i]]$functions[[j+1]]$parent==MoTBF.BN[[i]]$functions[[1]]$parent){
              data <- originalData
            }else{
              parents <- Parents[1:(j+1)]; pos <- 0
              first <- max(which(parents==MoTBF.BN[[i]]$functions[[1]]$parent))
              if(MoTBF.BN[[i]]$functions[[j]]$parent==MoTBF.BN[[i]]$functions[[j+1]]$parent){
                for(h in j:1){
                  if(MoTBF.BN[[i]]$functions[[h+1]]$parent==MoTBF.BN[[i]]$functions[[h]]$parent) pos <- h
                  else break
                }
                if(is.character(MoTBF.BN[[i]]$functions[[pos]]$interval)){
                  s <- which(names(which(sapply(data, is.character)==T))==MoTBF.BN[[i]]$functions[[pos]]$parent)
                  sta <- states[[s]]$states
                  estate <- which(sta==MoTBF.BN[[i]]$functions[[pos]]$interval)
                  if(estate!=1){
                    t <- -1; pos <- discretePos(pos, MoTBF.BN[[i]]$functions, sta, estate,first, t)
                  }
                } else {
                  t <- c(); pos <- fun1(pos, MoTBF.BN[[i]]$functions,first,j, t)
                }
                data <- fullData[[pos-1]]
                
              }else {
                
                if(is.character(MoTBF.BN[[i]]$functions[[j+1]]$interval)){
                  for(h in j:1){
                    if(MoTBF.BN[[i]]$functions[[j+1]]$parent==MoTBF.BN[[i]]$functions[[h]]$parent){
                      pos <- h; break
                    } 
                  }
                  s <- which(names(which(sapply(data, is.character)==T))==MoTBF.BN[[i]]$functions[[pos]]$parent)
                  sta <- states[[s]]$states
                  estate <- which(sta==MoTBF.BN[[i]]$functions[[pos]]$interval)
                  if(estate!=1){
                    t <- -1; pos <- discretePos(pos, MoTBF.BN[[i]]$functions, sta, estate,first, t)
                  }
                } else {
                  if(length(disc)==0) estates <- 0
                  else estates <- discreteVariablesStates(disc, fullData[[first]])
                  pos <- fun(parents, MoTBF.BN[[i]]$functions, first, j, pos, estates, fullData[[first]])
                }
                
                if(pos==1||pos==0) data <- originalData
                else data <- fullData[[pos-1]]
              }
            }
          }
        }
      }
      if((length(valuesT)-N)!=nrow(originalData)) valuesT <- duplicatedValues(fullData, valuesT, N, enter)
    }
  }
  loglike <- sum(log(valuesT))
  return(loglike)
}

#' @rdname goodnessMoTBFBN
#' @export
BiC.MoTBFBN <- function(MoTBF.BN, data)
{
  t <- logLikelihood.MoTBFBN(MoTBF.BN, data)
  param <- c()
  for(i in 1:length(MoTBF.BN)){
    if(length(MoTBF.BN[[i]]$functions)==1){
      if(is.motbf(MoTBF.BN[[i]]$functions[[1]])) param <- c(param,length(coef(MoTBF.BN[[i]]$functions[[1]]))+1)
      if(is.list(MoTBF.BN[[i]]$functions[[1]])){
        if(is.null(MoTBF.BN[[i]]$functions[[1]]$Px)) param <- c(param,length(MoTBF.BN[[i]]$functions[[1]]$coeff))
        else param <- c(param,length(coef(MoTBF.BN[[i]]$functions[[1]]$Px))+1)
      }
    }else{
      for(j in 1:length(MoTBF.BN[[i]]$functions)){
        if(is.motbf(MoTBF.BN[[i]]$functions[[j]]$Px)) param <- c(param,length(coef(MoTBF.BN[[i]]$functions[[j]]$Px))+1)
        else param <- c(param,length(MoTBF.BN[[i]]$functions[[j]]$Px$coeff))
      }
    }
  }
  
  BiC <- t - 1/2*(sum(param))*log(nrow(data)) 
  return(list(LogLikelihood=t, BIC=BiC))
}


fun <- function(parents, interval, first, j, pos, states, originalData)
{
  pos <- which(parents[(first+1):(j+1)]%in%interval[[first+1]]$parent)
  
  if((length(pos)!=1)&&(!is.null(pos))){
    no <- c()
    for(t in 1:(length(pos)-1)){
      if(is.character(interval[first+pos[t]][[1]]$interval)){
        s <- which(names(which(sapply(data, is.character)==T))==interval[first+pos[t+1]][[1]]$parent)
        sta <- states[[s]]$states
        estate1 <- which(sta==interval[first+pos[t]][[1]]$interval)
        estate2 <- which(sta==interval[first+pos[t+1]][[1]]$interval)
        if((estate1+1)==estate2) no <- c(no,pos[t], pos[t+1])
      }else {
        if(interval[first+pos[t]][[1]]$interval[2]==interval[first+pos[t+1]][[1]]$interval[1]) no <- c(no,pos[t], pos[t+1])
      }
    }
    if(any(duplicated(no))) no <- no[-which(duplicated(no))]
    if(is.null(no)) f1 <- first+1
    else f1 <- first+no[length(no)]
  } else {
    f1  <- first+1
  }
  
  if(f1==(j+1)) return(first + no[1])
  if(pos[1]==0) return(0)
  pos <- fun(parents, interval, f1, j, pos, states, originalData)
}

fun1 <- function(pos, interval, first,j, t)
{
  if(is.null(interval[[pos-1]]$Px)){
    return(pos)
  }else{
    for(h in (pos-1):first){
      if(is.null(interval[[h]]$Px)){
        pos <- h; break
      } 
    }
    if(pos==1) return(pos+1)
    if(interval[[pos]]$parent!=interval[[j+1]]$parent) return(pos+1)
    pos <- fun1(pos, interval,first,j, pos)
  }
}

discretePos <- function(pos, interval, sta, estate,first, t)
{
  for(h in (pos-1):first){
    if(interval[[pos]]$parent==interval[[h]]$parent) {
      if(interval[[h]]$interval!=sta[estate-1]) break
      else pos <- h
    }
  }
  if(estate==1) return(pos)
  if(sta[estate-1]==sta[1]) return(pos)
  if(h==1) return(pos)
  if(t==pos) {
    if(is.null(interval[[pos-1]]$Px)){
      return(pos)
    } else {
      sta <- sta[-estate]; pos <- discretePos(pos, interval, sta, estate-1, first, pos)
    }
  }
  
  pos <- discretePos(pos, interval, sta, estate-1,first, pos)
}

duplicatedValues <- function(fullData, valuesT, N, enter)
{
  coln=c()
  for(i in 1:length(fullData)){
    if(enter[i]=="NO") next
    coln <- c(coln,as.numeric(rownames(fullData[[i]])))
  }
  ind <- which(duplicated(coln, fromLast=T))
  if((length(ind)!=0)&&(all(ind>=0))) valuesT <- valuesT[-(N+ind)]
  return(valuesT)
}



