#BASED ON <<A (Not So) Short Introduction to S4>>
setClass (
  Class = "DecisionTable",
  representation = representation(
    decisionTable = "matrix"
  ),
  validity = function(object){
    #cat("~~~ DecisionTable: inspector ~~~ \n")
    if(nrow(object@decisionTable) < 2 || ncol(object@decisionTable) < 2){
      stop ("[DecisionTable: validation] the minimum number of row and columns are 2")
    }else{}
    return(TRUE)
  }
)
#*******************************************************
#CONSTRUCTOR
setMethod (
  f="initialize",
  signature="DecisionTable",
  definition=function(.Object,decisionTable){
    #cat ("~~~~~ DecisionTable: initializator ~~~~~ \n")
    
    if(!missing(decisionTable)){
      colnames(decisionTable) <- paste("C",1:ncol(decisionTable),sep="")
      dtColNames <- colnames(decisionTable)
      dtColNames[length(dtColNames)] <- "D"
      colnames(decisionTable) <- dtColNames
      rownames(decisionTable) <- paste("R",1:nrow(decisionTable),sep="")
      #THIS DOESN'T WORK ON TESTS
      #if(length(decisionTable)>0){      
      #decisionTable <- changeDtNames(decisionTable)
      #}
      .Object@decisionTable <- decisionTable
      validObject(.Object)# call of the inspector
    }else{
      #decisionTable <- matrix(nrow=0,ncol=0)
      .Object@decisionTable <- matrix(nrow=0,ncol=0)
    }
    return(.Object)
  }
)

#CONSTRUCTOR (USER FRIENDLY)
decisionTable <- function(theDecisionTable){
  #cat ("~~~~~ DecisionTable: constructor ~~~~~ \n")
  new (Class="DecisionTable",decisionTable=theDecisionTable)
}


#*******************************************************
#ACCESSORS

#Returns the complete decision table as numeric matrix
setGeneric("getDecisionTable",function(object){standardGeneric ("getDecisionTable")})
setMethod("getDecisionTable","DecisionTable",
          function(object){
            return(object@decisionTable)
          }
)


#Returns the conditions of the decision table as numeric matrix
setGeneric("getCondition",function(object){standardGeneric ("getCondition")})
setMethod("getCondition","DecisionTable",
          function(object){
            decisionTable <- object@decisionTable
            condition <- decisionTable[,1:(ncol(decisionTable)-1)]
            return(condition)
          }
)

#Returns the decision of the decision table as vector
setGeneric("getDecision",function(object){standardGeneric ("getDecision")})
setMethod("getDecision","DecisionTable",
          function(object){
            decisionTable <- object@decisionTable
            decision <- decisionTable[,ncol(decisionTable)]
            return(decision)
          }
)

#Returns a subset of rules as numeric matrix
setGeneric("getRule",function(object,ruleIndex){standardGeneric ("getRule")})
setMethod("getRule","DecisionTable",
          function(object,ruleIndex){
            return(object@decisionTable[ruleIndex,])
          }
)


#*******************************************************
#GENERIC METODS

#Summary

setMethod ("print","DecisionTable",
           function(x,...){
             cat("*** Class DecisionTable, method Print *** \n")
             if(length(x@decisionTable) != 0){
               print(formatC(x@decisionTable),quote=FALSE)
             }else{}
             cat("******* End Print (DecisionTable) ******* \n")
           }
)

setMethod("show","DecisionTable",
          function(object){
            cat("*** Class DecisionTable, method Show *** \n")
            cat("* decisionTable (limited to a matrix 10x10) = \n")
            if(length(object@decisionTable) != 0){
              nrowShow <- min(10,nrow(object@decisionTable))
              ncolShow <- min(10,ncol(object@decisionTable))
              print(formatC(object@decisionTable[1:nrowShow,1:ncolShow]),quote=FALSE)
            }else{}            
            cat("******* End Show (DecisionTable) ******* \n")
          }
)

#*******************************************************
#METHODS

### checkConsistency takes
### - A DecisionTable
###
### checkConsistency returns a boolean vector indicating which rules are inconsistent or contradictory in the decision table given. Its a summary of the consistency matrix
### The consistency vector is the row result of the consistency matrix
###
setGeneric (name = "checkConsistency",def = function(object){standardGeneric("checkConsistency")})
setMethod(
  f = "checkConsistency",
  signature = "DecisionTable",
  definition = function(object){
    
    consistencyMatrix <- .computeConsistencyMatrix(object@decisionTable)
    expectedSum <- c(1:nrow(consistencyMatrix))#If the rule is consistent its sum will be equal to its row number
    obtainedSum <- apply(consistencyMatrix,1,function(theRule){sum(theRule,na.rm = TRUE)})
    consistentRules <- (obtainedSum == expectedSum)
    return(consistentRules)
  }
)


### computeConsistencyMatrix takes
### - A DecisionTable
###
### computeConsistencyMatrix returns a boolean diagonal matrix indicating inconsistency between rules. It must be interpreted by columns
###
setGeneric (name = "computeConsistencyMatrix",def = function(object){standardGeneric("computeConsistencyMatrix")})
setMethod(
  f = "computeConsistencyMatrix",
  signature = "DecisionTable",
  definition = function(object){
    
    consistencyMatrix <- .computeConsistencyMatrix(object@decisionTable)
    return(consistencyMatrix)
  }
)


### computeDiscernibilityMatrix takes
### - A DecisionTable
###
### computeDiscernibilityMatrix returns an object of the class DiscernibilityMatrix
###
setGeneric (name = "computeDiscernibilityMatrix",def = function(object){standardGeneric("computeDiscernibilityMatrix")})
setMethod(
  f = "computeDiscernibilityMatrix",
  signature = "DecisionTable",
  definition = function(object){
    
    discernibilityMatrix <- .computeDiscernibilityMatrix(object@decisionTable)#Duplicated rules could be removed before to speed performance
    dt <- new(Class="DiscernibilityMatrix",discernibilityMatrix = discernibilityMatrix)
    return(dt)
  }
)


### removeDuplicatedRules takes
### - A DecisionTable
###
### removeDuplicatedRulesDT returns a new decision table object without duplicated rules
###
setGeneric (name = "removeDuplicatedRulesDT",def = function(object){standardGeneric("removeDuplicatedRulesDT")})
setMethod(
  f = "removeDuplicatedRulesDT",
  signature = "DecisionTable",
  definition = function(object){
    
    dtMat <- object@decisionTable
    newDecisionTable <- unique(dtMat)
    newDT <- new(Class="DecisionTable",decisionTable = newDecisionTable)
    return(newDT)
  }
)


### findFirstConditionReduct takes
### - A DecisionTable
###
### findFirstConditionReduct returns one condition reduct object with the least number of conditions. This reduct belongs to the family of the smallest reducts of the decision table. It maybe the only one
###
setGeneric (name = "findFirstConditionReduct",def = function(object){standardGeneric("findFirstConditionReduct")})
setMethod(
  f = "findFirstConditionReduct",
  signature = "DecisionTable",
  definition = function(object){
    
    dtMat <- object@decisionTable
    dm <- computeDiscernibilityMatrix(object)
    coreNoDecision <- computeCore(dm)
    firstReduct <- .findFirstReductsFromCore(dtMat,coreNoDecision)
    firstReductWithDecision <- append(firstReduct,ncol(dtMat))#add the decision column id
    cr <- new(Class="ConditionReduct",decisionTable = object, columnIds = firstReductWithDecision)
    return(cr)
  }
)


### findSmallestReductFamilyFromCore returns a list of ConditionReduct objects representing the smallest reducts, all of them with the same number of conditions
###
setGeneric (name = "findSmallestReductFamilyFromCore",def = function(object){standardGeneric("findSmallestReductFamilyFromCore")})
setMethod(
  f = "findSmallestReductFamilyFromCore",
  signature = "DecisionTable",
  definition = function(object){
    lres <- list()
    dtMat <- object@decisionTable
    dm <- computeDiscernibilityMatrix(object)
    coreNoDecision <- computeCore(dm)
    reductFamilyList <- .findSmallestReductFamilyFromCore(dtMat,coreNoDecision)
    listLen <- length(reductFamilyList)
    for(i in 1:listLen){
      reductVector <- reductFamilyList[[i]]
      reductVector <- append(reductVector,ncol(dtMat))#add the decision column id
      cr <- new(Class="ConditionReduct",decisionTable = object,columnIds = reductVector)
      lres[[i]] <- cr
    }
    return(lres)
  }
)


### findAllReductsFromCore returns a list of ConditionReduct objects representing all the reducts found in the decision table
###
setGeneric (name = "findAllReductsFromCore",def = function(object){standardGeneric("findAllReductsFromCore")})
setMethod(
  f = "findAllReductsFromCore",
  signature = "DecisionTable",
  definition = function(object){
    lres <- list()
    dtMat <- object@decisionTable
    dm <- computeDiscernibilityMatrix(object)
    coreNoDecision <- computeCore(dm)
    reductFamilyList <- .findAllReductsFromCore(dtMat,coreNoDecision)
    listLen <- length(reductFamilyList)
    for(i in 1:listLen){
      reductVector <- reductFamilyList[[i]]
      reductVector <- append(reductVector,ncol(dtMat))#add the decision column id
      cr <- new(Class="ConditionReduct",decisionTable = object,columnIds = reductVector)
      lres[[i]] <- cr
    }
    return(lres)
  }
)


### simplifyDecisionTable returns a value reduct
###
setGeneric (name = "simplifyDecisionTable",def = function(object){standardGeneric("simplifyDecisionTable")})
setMethod(
  f = "simplifyDecisionTable",
  signature = "DecisionTable",
  definition = function(object){
    #REDUCTION STEPS - BOOK PG 71
    #1 - Computation of reducts of condition attributes which is equivalent to elimination of some column from the decision table
    #2 - Elimination of duplicated rows
    #3 - Elimination of superflous values of attributes
    
    
    #ALBER IMPLEMENTATION FOR MINIMUM/OPTIMAL REDUCT
    #1 - Elimination of duplicated rows
    #2 - Compute discerniility matrix and obtain CORE
    #3a- Check if CORE is a reduct (the CORE is not always a reduct)
    #3b- Find the first smallest reduct
    #3c- Find all the smallest reducts
    #3d- Find all reducts
    #4 - Elimination of superflous values of attributes of reduct(s)(find value reduct)
    
    #Elimination of duplicated rows
    noDuplicatedDT <- removeDuplicatedRulesDT(object)
    #Compute discerniility matrix and obtain CORE
    #Check if CORE is a reduct (the CORE is not always a reduct)
    #Find the first smallest reduct
    cr <- findFirstConditionReduct(noDuplicatedDT)
    noDuplicatedCR <- removeDuplicatedRulesCR(cr)
    #Find value reduct
    vr <- computeValueReduct(noDuplicatedCR)
    NoDuplicatedVr <- removeDuplicatedRulesVR(vr)
    return(NoDuplicatedVr)
  }
)

#*******************************************************
#UTIL


### .inflateValueReduct
### - theValueReductWithDecision is a numeric matrix of the value reduct including an additional column representing the decision's table decision
### - theDecisionTableConditionDecisionIds is a numeric vector with the condition and decision (column) indexes (positions) of the reduct in the decition table
###
### .inflateValueReduct returns a numeric matrix representing th reduced knowge of the original decision table calculated from the given value reduct
###
.inflateValueReduct <- function(theValueReductWithDecision, theReductConditionIds,theDecisionTableConditionDecisionIds){
  resRowCounter <- nrow(theValueReductWithDecision)#Value reduct number of rules is different from the decision table number of rules
  resColCounter <- length(theDecisionTableConditionDecisionIds)# An additional column for the decision
  res <- matrix(data = NA, nrow = resRowCounter, ncol = resColCounter)
  for(i in 1:resRowCounter){
    valueReductRule <- theValueReductWithDecision[i,]
    copyCounter <- 0
    for(j in 1:resColCounter){
      if(is.element(j, theReductConditionIds)){
        copyCounter <- (copyCounter + 1)
        res[i,j] <- valueReductRule[copyCounter]
      }else{}
      if(j == resColCounter){
        res[i,j] <- valueReductRule[length(valueReductRule)]#Always copies the decision
      }else{}
    }
  }
  return(res)
}


### .findFirstReductsFromCore
### - theDecisionTable is a matrix where each row is a rule. A rule is composed of condition(all elements except the last one) and decision(the last vector element)
### - theCoreNoDecision is a numeric vector  indicating the decision table columns id which are part of the core without the decision column id
###
### .findFirstReductsFromCore returns a numeric vector with condition (column) indexes which are a reduct of theDecisionTable
###
.findFirstReductsFromCore <- function(theDecisionTable,theCoreNoDecision){
  #make the combinations in case a reduct is not obtained just by adding 1 condition
  
  theCore <- theCoreNoDecision
  if(.isReduct(theDecisionTable, theCore)){#Tests if the core is a reduct itself
    res <- theCore
  }else{
    conditionCount <- (ncol(theDecisionTable) - 1)
    notCoreConditionsId <- setdiff(c(1:conditionCount),theCore)# COnditions which are not part of the core
    res <- vector(mode = "numeric", length = 0)
    
    for(i in 1:(length(notCoreConditionsId))){
      if(length(notCoreConditionsId) > 1){
        combinations <- combn(notCoreConditionsId,i)#All the cOmbinations of not core conditions  
        combinationsCount <- ncol(combinations)
      }else if(length(notCoreConditionsId) == 1){
        combinations <- notCoreConditionsId
        combinationsCount <- 1
      }else{}
      
      for(j in 1:combinationsCount){
        if(is.vector(combinations)){
          combinationNoCore <- combinations[j]#Takes one possible not core condition combination at the time
        }else if(is.matrix(combinations)){
          combinationNoCore <- combinations[,j]#Takes one possible not core condition combination at the time
        }else{}
        combinationCore <- union(theCore, combinationNoCore)#Adds the core to the no core conditions
        combinationCore <- sort(combinationCore)
        if(.isReduct(theDecisionTable, combinationCore)){
          res <- combinationCore
          break #Makes sure it returns just one, the first found solution
        }else{}
      }
      if(length(res) > 0){
        break #Stop testing combinations
      }else{}
    }
  }
  return(res)
}





### .findSmallestReductFamilyFromCore
### - theDecisionTable is a matrix where each row is a rule. A rule is composed of condition(all elements except the last one) and decision(the last vector element)
### - theCoreNoDecision is a numeric vector  indicating the decision table columns id which are part of the core without the decision column id
###
### .findSmallestReductFamilyFromCore returns a list of numeric vectors indicating the reducts with the least number of conditions possible
###
.findSmallestReductFamilyFromCore <- function(theDecisionTable,theCoreNoDecision){
  #make the combinations in case a reduct is not obtained just by adding 1 condition
  
  theCore <- theCoreNoDecision
  lres <- list()
  if(.isReduct(theDecisionTable, theCore)){#Tests if the core is a reduct itself
    lres[[1]] <- theCore
  }else{
    conditionCount <- (ncol(theDecisionTable) - 1)
    notCoreConditionsId <- setdiff(c(1:conditionCount),theCore)# COnditions which are not part of the core
    
    counter <- 0
    for(i in 1:(length(notCoreConditionsId))){
      if(length(notCoreConditionsId) > 1){
        combinations <- combn(notCoreConditionsId,i)#All the cOmbinations of not core conditions  
        combinationsCount <- ncol(combinations)
      }else if(length(notCoreConditionsId) == 1){
        combinations <- notCoreConditionsId
        combinationsCount <- 1
      }else{}
      
      for(j in 1:combinationsCount){
        if(is.vector(combinations)){
          combinationNoCore <- combinations[j]#Takes one possible not core condition combination at the time
        }else if(is.matrix(combinations)){
          combinationNoCore <- combinations[,j]#Takes one possible not core condition combination at the time
        }else{}
        combinationCore <- union(theCore, combinationNoCore)#Adds the core to the no core conditions
        combinationCore <- sort(combinationCore)
        if(.isReduct(theDecisionTable, combinationCore)){
          counter <- (counter + 1)
          lres[[counter]] <- combinationCore
        }else{}
      }
      if(counter != 0){
        break#Makes sure it returns only the smallest family of condition reducts
      }
    }
  }
  return(lres)
}


### .findAllReductsFromCore
### - theDecisionTable is a matrix where each row is a rule. A rule is composed of condition(all elements except the last one) and decision(the last vector element)
### - theCoreNoDecision is a numeric vector  indicating the decision table columns id which are part of the core without the decision column id
###
### .findAllReductsFromCore returns a list of numeric vectors indicating the all the reducts in the decision table
###
.findAllReductsFromCore <- function(theDecisionTable,theCoreNoDecision){
  #make the combinations in case a reduct is not obtained just by adding 1 condition
  
  theCore <- theCoreNoDecision
  lres <- list()
  counter <- 0
  if(.isReduct(theDecisionTable, theCore)){#Tests if the core is a reduct itself
    counter <- (counter + 1)
    lres[[counter]] <- theCore
  }else{}
  conditionCount <- (ncol(theDecisionTable) - 1)
  notCoreConditionsId <- setdiff(c(1:conditionCount),theCore)# COnditions which are not part of the core
  
  for(i in 1:(length(notCoreConditionsId))){
    if(length(notCoreConditionsId) > 1){
      combinations <- combn(notCoreConditionsId,i)#All the cOmbinations of not core conditions  
      combinationsCount <- ncol(combinations)
    }else if(length(notCoreConditionsId) == 1){
      combinations <- notCoreConditionsId
      combinationsCount <- 1
    }else{}
    
    for(j in 1:combinationsCount){
      if(is.vector(combinations)){
        combinationNoCore <- combinations[j]#Takes one possible not core condition combination at the time
      }else if(is.matrix(combinations)){
        combinationNoCore <- combinations[,j]#Takes one possible not core condition combination at the time
      }else{}
      combinationCore <- union(theCore, combinationNoCore)#Adds the core to the no core conditions
      combinationCore <- sort(combinationCore)
      if(.isReduct(theDecisionTable, combinationCore)){
        counter <- (counter + 1)
        lres[[counter]] <- combinationCore
      }else{}
    }
  }
  return(lres)
}















### .computeDiscernibilityMatrix
### - theDecisionTable is a matrix where each row is a rule. A rule is composed of condition(all elements except the last one) and decision(the last vector element)
###
###
### .computeDiscernibilityMatrix returns an array of rules x rules x conditions according to "ROUGH SETS - Theoretical aspects of reasoning about data" section 5.5, page 60
###
.computeDiscernibilityMatrix <- function(theDecisionTable){
  id <- c(1:nrow(theDecisionTable))#Rule id
  conditionCount <- ncol(theDecisionTable) - 1
  ruleCount <- nrow(theDecisionTable)
  
  theDTid <- cbind(theDecisionTable,id)#Attach rule id to decision table for APPLY
  partialResult <- apply(theDTid,1,.computeDiscernibilityVector,theDecisionTable=theDecisionTable)
  
  #Takes the partialResult and builds the expected array
  discernibilityMatrix <- array(NA,c(ruleCount,ruleCount,conditionCount))
  mat <- matrix()
  for(i in 1:ruleCount){
    mat <- matrix(partialResult[,i],ncol=conditionCount,byrow=TRUE)
    for(j in 1:nrow(mat)){
      discernibilityMatrix[i,j,] <- mat[j,]
    }
  }
  
  ruleNames <- paste("R",1:nrow(discernibilityMatrix),sep="")
  attNames <- paste("C",1:conditionCount,sep="")
  dimensionNames <- list(ruleNames,ruleNames,attNames)
  dimnames(discernibilityMatrix) <- dimensionNames
  return(discernibilityMatrix)  
}


### .computeConsistencyMatrix
### - theDecisionTable is a matrix where each row is a rule. A rule is composed of condition(all elements except the last one) and decision(the last vector element)
###
###
### .computeConsistencyMatrix returns a diagonal matrix of booleans indicating the consistency beetwen rules
###
.computeConsistencyMatrix <- function(theDecisionTable){
  id <- c(1:nrow(theDecisionTable))#Rule id
  theDTid <- cbind(theDecisionTable,id)#Attach rule id to decision table for APPLY
  consistencyMatrix <- apply(theDTid,1,.computeConsistencyVector,theDecisionTable=theDecisionTable)
  rownames(consistencyMatrix) <- paste("R",1:nrow(consistencyMatrix),sep="")
  colnames(consistencyMatrix) <- paste("R",1:ncol(consistencyMatrix),sep="")
  return(consistencyMatrix)  
}


### .computeDiscernibilityVector
### - theRule is a vector where the las element is considered to be the rule ID inside theDecisionTable and the remaining vector elements are the rule 
### - theDecisionTable is a matrix where each row is a rule. A rule is composed of condition(all elements except the last one) and decision(the last vector element)
###
### .computeDiscernibilityVector returns a boolean matrix comparing elementwise the rule against the Decision table's rules (TRUE means elements are different)
###
.computeDiscernibilityVector<- function(theRule,theDecisionTable){
  #OPTIMISATION: Just compute discernibility for one rule on each decision category  
  
  ruleId <- theRule[length(theRule)]#The last element is the rule ID, the remaining are the rule
  rule <- theRule[-length(theRule)]#The last element is the rule ID, the remaining are the rule
  ruleCondition <- rule[-length(rule)]#The last column is the decision, the rest conditions
  ruleConditionCount <- length(ruleCondition)
  conditionTable <- theDecisionTable[,-ncol(theDecisionTable)]#The last column are the decisions, the rest conditions
  ruleCount <- nrow(theDecisionTable)
  conditionCount <- length(rule) - 1#The last element is the decision, the remaining are the condition
  
  #theDisTable <- vector(mode = "logical", length = (ruleCount*conditionCount))
  theDisTable <- array(NA, c(ruleCount,conditionCount))
  
  for(i in 1:ruleCount){
    iniPos <- ((i * ruleConditionCount)-(ruleConditionCount-1))
    endPos <-  (i * ruleConditionCount)
    if(i <= ruleId){
      theDisTable[iniPos:endPos] <- NA
    }else if(i > ruleId){
      theDisTable[iniPos:endPos] <- (ruleCondition != conditionTable[i,])
    }else{}
  }
  #m <- matrix(theDisTable,nrow=ruleCount,ncol=conditionCount,byrow = TRUE)#This returns the result as a matrix
  return(theDisTable)
}


### .computeConsistencyVector
### - theRule is a vector where the las element is considered to be the rule ID inside theDecisionTable and the remaining vector elements are the rule 
### - theDecisionTable is a matrix where each row is a rule. A rule is composed of condition(all elements except the last one) and decision(the last vector element)
###
### .computeConsistencyVector returns a boolean vector indicating the consistency of a rule with the following rules in the decision table
###
.computeConsistencyVector <- function(theRule,theDecisionTable){
  ruleId <- theRule[length(theRule)]#The last element is the rule ID, the remaining are the rule
  rule <- theRule[1:(length(theRule)-1)]
  ruleCount <- nrow(theDecisionTable)
  theConVector <- vector(mode = "logical", length = ruleCount)
  theValue <- NA
  for(i in 1:ruleCount){
    if(i < ruleId){
      theValue <- NA
    }else if(i == ruleId){
      theValue <- TRUE# A rule is consistent with itself
    }else if(i > ruleId){
      theValue <- .checkRuleConsistency(rule,theDecisionTable[i,])
    }else{}
    theConVector[i] <- theValue
  }
  #If there is an inconsistency in the vector, then the rule is inconsistent  too
  if(sum(theConVector, na.rm = TRUE) != ((length(theConVector) + 1) - ruleId)){
    theConVector[ruleId] <- FALSE
  }else{}
  return(theConVector)
}


### .checkRuleConsistency
### - rule1 is a vector where the las element is considered to be a decision an the rest conditions
### - rule2 is a vector where the las element is considered to be a decision an the rest conditions
###
### .checkRuleConsistency returns a boolean indicating if the rules are consistent
###
.checkRuleConsistency <- function(rule1,rule2){
  areConsistent <- FALSE
  condR1 <- rule1[-length(rule1)]
  condR2 <- rule2[-length(rule2)]
  decR1 <- rule1[length(rule1)]
  decR2 <- rule2[length(rule2)]
  if(identical(condR1,condR2)){
    if(identical(decR1,decR2)){
      areConsistent <- TRUE# rule1 & rule2 have exactly the same condition and decision
    }else{
      areConsistent <- FALSE# Same condition, different decision
    }
  }else{
    areConsistent <- TRUE# Different condition, the decision doesn't matter
  }
  return(areConsistent)
}


### .isReduct
### - theDecisionTable is a matrix where each row is a rule. A rule is composed of condition(all elements except the last one) and decision(the last vector element)
### - theReductConditionIds is a numeric vector indicating the indexes of the columns (conditions) of the reduct in the decision table
###
### .isReduct returns a boolean indicating if theReductConditionIds are or not a reduct of the given decision table  
###
.isReduct <- function(theDecisionTable, theReductConditionIds){
  reduct <- theDecisionTable[,theReductConditionIds]#Reduct building from decision table
  conditions <- theDecisionTable[,-ncol(theDecisionTable)]#Removes the Decision from the decision table
  indRed <- .computeIndescernibilityFunction(reduct)
  indCon <- .computeIndescernibilityFunction(conditions)
  ans <- identical(indCon,indRed)
  return(ans)
}