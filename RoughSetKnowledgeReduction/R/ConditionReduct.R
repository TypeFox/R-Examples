setClass (
  Class = "ConditionReduct",
  representation = representation(
    decisionTable = "DecisionTable",
    columnIds = "numeric"
  ),
  validity = function(object){
    #cat("~~~ ConditionReduct: inspector ~~~ \n")
    if(length(object@columnIds) < 2){
      stop ("[ConditionReduct: validation] the minimum number of columnIds Id is 2 (1 condition and the decision)")
    }else{}
    
    len <- length(object@columnIds)
    dt <- object@decisionTable
    dtMat <- getDecisionTable(dt)
    columnCount <- ncol(dtMat)
    
    if(len > 0){
      decisionFlag <- FALSE
      for(i in 1:len){
        el <- object@columnIds[i]
        if(el > columnCount){
          stop ("[ConditionReduct: validation] a column Id is out of the decision table")        
        }else{}
        if(el == columnCount){
          decisionFlag <- TRUE
        }else{}
      }
      if(decisionFlag == FALSE){#Check that the decision Id is included
        stop ("[ConditionReduct: validation] a column Id pointing to the decision is needed")
      }
    }
    return(TRUE)
  }
)

#*******************************************************
#CONSTRUCTOR
setMethod (
  f="initialize",
  signature="ConditionReduct",
  definition=function(.Object,decisionTable,columnIds){
    #cat ("~~~~~ ConditionReduct: initializator ~~~~~ \n")
    if(!missing(decisionTable)){
      if(!missing(columnIds)){      
        .Object@decisionTable <- decisionTable
        .Object@columnIds <- columnIds
        validObject(.Object)# call of the inspector
      }else{
        .Object@columnIds <- vector(mode = "numeric", length = 0)
      }  
    }else{
      .Object@decisionTable <- new(Class="DecisionTable")
    }
    return(.Object)
  }
)


#CONSTRUCTOR (USER FRIENDLY)
conditionReduct <- function(theDecisionTable,theColumnIds){
  #cat ("~~~~~ ConditionReduct: constructor ~~~~~ \n")
  new (Class="ConditionReduct", decisionTable = theDecisionTable, columnIds = theColumnIds)
}

#*******************************************************
#ACCESSORS

#Returns the Decision Table object
setGeneric("getConditionReductDecisionTable",function(object){standardGeneric ("getConditionReductDecisionTable")})
setMethod("getConditionReductDecisionTable","ConditionReduct",
          function(object){
            return(object@decisionTable)
          }
)

#Returns the numeric vector column Ids
setGeneric("getColumnIds",function(object){standardGeneric ("getColumnIds")})
setMethod("getColumnIds","ConditionReduct",
          function(object){
            return(object@columnIds)
          }
)

#Returns the condition reduct inflated as a numeric matrix
setGeneric("getConditionReduct",function(object){standardGeneric ("getConditionReduct")})
setMethod("getConditionReduct","ConditionReduct",
          function(object){
            cids <- object@columnIds
            dt <- object@decisionTable
            dtMat <- getDecisionTable(dt)
            res <- .buildConditionReductMatrix(dtMat,cids)
            return(res)
          }
)

#*******************************************************
#GENERIC METODS

#Summary

setMethod ("print","ConditionReduct",
           function(x,...){
             cat("*** Class ConditionReduct, method Print *** \n")
             cids <- x@columnIds
             dt <- x@decisionTable
             dtMat <- getDecisionTable(dt)
             if(length(dtMat) != 0){
               printMat <- .buildConditionReductMatrix(dtMat,cids)
               print(formatC(printMat),quote=FALSE)
             }else{}
             cat("******* End Print (ConditionReduct) ******* \n")
           }
)


setMethod("show","ConditionReduct",
          function(object){
            cat("*** Class ConditionReduct, method Show *** \n")
            cat("* ConditionReduct (limited to a matrix 10x10) = \n")
            cids <- object@columnIds
            dt <- object@decisionTable
            dtMat <- getDecisionTable(dt)
            if(length(dtMat) != 0){
              printMat <- .buildConditionReductMatrix(dtMat,cids)
              nrowShow <- min(10,nrow(printMat))
              ncolShow <- min(10,ncol(printMat))
              print(formatC(printMat[1:nrowShow,1:ncolShow]),quote=FALSE)
            }else{}
            cat("******* End Show (ConditionReduct) ******* \n")
          }             
)


#*******************************************************
#METHODS


### isConditionReduct returns a boolean indicating if the condition reduct is a condition reduct of the decision table
setGeneric (name = "isConditionReduct",def = function(object){standardGeneric("isConditionReduct")})
setMethod(
  f = "isConditionReduct",
  signature = "ConditionReduct",
  definition = function(object){
    
    cids <- object@columnIds
    dt <- object@decisionTable
    dtMat <- getDecisionTable(dt)
    res <- .isConditionReduct(dtMat,cids)
    return(res)
  }
)



### computeValueReduct returns an object of type ValueReduct
setGeneric (name = "computeValueReduct",def = function(object){standardGeneric("computeValueReduct")})
setMethod(
  f = "computeValueReduct",
  signature = "ConditionReduct",
  definition = function(object){
    
    cids <- object@columnIds
    dt <- object@decisionTable
    dtMat <- getDecisionTable(dt)
    decisionTableColumnCounter <- ncol(dtMat)
    reductWithDecision <- dtMat[,cids]
    valueReductRaw <- .computeValueReducts(reductWithDecision)
    valueReductOriginRuleId <- valueReductRaw[,ncol(valueReductRaw)]#Gets the column with the origin rule Id
    valueReduct <- valueReductRaw[,-ncol(valueReductRaw)]#Removes the column with the origin rule Id
    valueReductInflated <- .inflateValueReductMatrix(valueReduct, cids, decisionTableColumnCounter)
    #ok
    cn <- colnames(dtMat)
    colnames(valueReductInflated) <- cn
    rownames(valueReductInflated) <- paste("R",valueReductOriginRuleId,sep="")
    vr <- new(Class="ValueReduct", conditionReduct = object, valueReduct = valueReductInflated)
    return(vr)
  }
)



### removeDuplicatedRulesCR returns a new condition reduct object without duplicated rules in the decision table from the condition reduct perspective
setGeneric (name = "removeDuplicatedRulesCR",def = function(object){standardGeneric("removeDuplicatedRulesCR")})
setMethod(
  f = "removeDuplicatedRulesCR",
  signature = "ConditionReduct",
  definition = function(object){
    
    dt <- object@decisionTable
    cids <- object@columnIds
    dtMat <- getDecisionTable(dt)
    
    crMat <- dtMat[,cids]
    dup <- duplicated(crMat)
    posVec <- .boolean2positionVector(dup)
    rowCount <- nrow(dtMat)
    posVecDif <- setdiff(c(1:rowCount),posVec)
    newDtMat <- dtMat[posVecDif,]
    
    newDt <- new(Class="DecisionTable",decisionTable = newDtMat)
    newCr <- new(Class="ConditionReduct",decisionTable = newDt,columnIds = cids)
    return(newCr)
  }
)


#*******************************************************
#UTIL


### .buildConditionReductMatrix
### - theDecisionTableMatrix is a numeric matrix representing a decision table
### - theColumnIds is a numeric vector representing the column Ids which conform the reduct. The decision Id is needed
###
###
### .buildConditionReductMatrix returns a formated numeric matrix where the columns not part of theColumnIds are replaced by NA
###
.buildConditionReductMatrix <- function(theDecisionTableMatrix,theColumnIds){  
  
  dtMat <- theDecisionTableMatrix
  colCount <- ncol(dtMat)
  rowCount <- nrow(dtMat)
  for(i in 1:colCount){
    if(!is.element(i,theColumnIds)){
      dtMat[,i] <- NA
    }
  }
  return(dtMat)
}


### .isConditionReduct
### - theDecisionTable is a matrix where each row is a rule. A rule is composed of condition(all elements except the last one) and decision(the last vector element)
### - theReductConditions is a numeric vector indicating the indexes of the columns (conditions) of the reduct in the decision table
###
### .isConditionReduct returns a boolean indicating if theReductConditions are or not a reduct of the given decision table  
###
.isConditionReduct <- function(theDecisionTable, theReductConditions){
  
  cids <- theReductConditions[-length(theReductConditions)]#Removes the decision from the columns Ids
  reduct <- theDecisionTable[,cids]#Reduct taken from decision table
  conditions <- theDecisionTable[,-ncol(theDecisionTable)]#Removes the Decision from the decision table
  indRed <- .computeIndescernibilityFunction(reduct)
  indCon <- .computeIndescernibilityFunction(conditions)
  ans <- identical(indCon,indRed)
  return(ans)
}


### .computeIndescernibilityFunction
### - theRows is a numeric vector/matrix where each row is a rule or a rule's condition or decision
###
### .computeIndescernibilityFunction returns a numeric vector indicating the Indescernibility Function; that's is the Row Id where a row is duplicated evaluated top-down
###
.computeIndescernibilityFunction <- function(theRows){
  
  if(is.vector(theRows)){
    theRows <- matrix(theRows,ncol=1)
  }else{}
  
  rowCount <- nrow(theRows)  
  result <- vector(mode = "numeric", length = rowCount)
  
  for(i in 1:rowCount){
    theCategory <- result[i]
    if(theCategory == 0){#The category hasn?t been assigned yet
      theCategory <- i#The category is the number of the row
      referenceRow <- theRows[i,]
      for(j in i:rowCount){
        testedRow <- theRows[j,]
        if(identical(referenceRow,testedRow)){
          result[j] <- i
        }else{}
      }
    }else{}
  }
  return(result)
}


### .computeValueReducts
### - theReductWithDecision is a reduct (condition reduct) of a decision table
###
### .computeValueReducts returns a matrix with rules derived from the calculation of attribute reducts of theReductWithDecision. The last column indicates the rule from which the reduced rule was generated
###
.computeValueReducts <- function(theReductWithDecision){
  
  ruleCount <- nrow(theReductWithDecision)
  ruleIds <- as.matrix(c(1:ruleCount))
  columnCount <- ncol(theReductWithDecision)
  
  #********************************************************************************************************************************************  
  #PARALELL COMPUTING - SNOWFALL
  #Get the value reducts of each rule in theReductWithDecision
  #********************************************************************************************************************************************
  valueReducts <- apply(ruleIds,1,.computeValueReduct,theReductWithDecision)
  #********************************************************************************************************************************************
  
  reducedRulevector <- vector(mode = "numeric",length = 0)
  #Process the results and build rules  
  for (i in 1:length(valueReducts)){
    valueReductMatrix <- valueReducts[[i]][[1]]#The 1st list element is teh matrix which columns are conditions indexes which are value reducts
    valueReductCounter <- ncol(valueReductMatrix)
    rule <- theReductWithDecision[i,]#Since the returned reducts are the smallest, the number of elements on the list is the same of rules in theReductWithDecision
    for(j in 1:valueReductCounter){
      valueReductConditionIndex <- valueReductMatrix[,j]
      reducedRule <- NA
      for(k in 1:(columnCount + 1)){
        if(is.element(k,valueReductConditionIndex)){
          reducedRule[k] <- rule[k]
        }else{} 
        if(k == columnCount){
          reducedRule[k] <- rule[k]#The decision is copied always
        }else if(k == (columnCount + 1)){
          reducedRule[k] <- i
        }else{}
      }
      reducedRulevector <- append(reducedRulevector,reducedRule,after=length(reducedRulevector))
    }
  }
  reducedRuleMatrix <- matrix(reducedRulevector,ncol=(columnCount + 1),byrow = TRUE)
  return(reducedRuleMatrix)
}


### .computeValueReduct
### - theReductWithDecision is a reduct (condition reduct) of a decision table
### - theRuleId is the Id of the rule (row id) in the theReductWithDecision
###
### .computeValueReduct returns a list of numeric matrixes where each column represents a a vector of condition indexs which are attribute reducts of the theRuleId 
###
.computeValueReduct <- function (theRuleId,theReductWithDecision){
  
  rule <- theReductWithDecision[theRuleId,]#Get the rule
  colCount <- length(rule)#Conditions and decision
  ruleCount <- nrow(theReductWithDecision)#Number of rules
  
  #Computes value indescernibility function for al columns (conditions & decision) in theReductWithDecision
  discVecMat <- matrix(NA,nrow = ruleCount,ncol = colCount)
  for(j in 1:colCount){
    cond <- theReductWithDecision[,j]
    discVecMat[,j] <- (cond == rule[j])
  }
  
  #Separates the conditions and decision for further processing
  decisionInd <- discVecMat[,ncol(discVecMat)]
  conditionInd <- discVecMat[,-ncol(discVecMat)]
  
  indDeci <- .boolean2positionVector(decisionInd)#Transform the boolean decision into position indexes
  
  #Test combinations of conditions which may become a value reduct
  lvalueReducts <- list()
  listElementIndex <- 0
  condCount <- ncol(conditionInd)#Condition number
  for(i in 1:condCount){
    combinations <- combn(condCount,i)#All the posible combinations of conditions
    combCount <- ncol(combinations)#Number of possible combinations given a number of conditions
    
    #********************************************************************************************************************************************  
    #PARALELL COMPUTING - SNOWFALL
    #Test every condition combination
    #********************************************************************************************************************************************
    partialResults <- apply(combinations,2,.testValueReduct,theConditionInd=conditionInd,theDecisionInd=indDeci)#Each column is a possible combinations of reduct conditions
    #********************************************************************************************************************************************    
    
    trueCounter <- sum(partialResults)
    #Add the sucessfull (.isValueReduct == TRUE) combinations to a list
    if(trueCounter > 0){
      results <- matrix(combinations[,partialResults],ncol=trueCounter)      
      listElementIndex <- listElementIndex + 1
      lvalueReducts[[listElementIndex]] <- results#Each sucesful combination is a column
      break#Ensures to get only the shortest
    }else{}
  }
  return(lvalueReducts)
}

### .position2booleanVector
### - thePositionVector is a numeric vector indicating the positions on the returned boolean vector which should be TRUE
### - booleanVectorLength is a number indicating the length of the returned boolean vector
###
### .position2booleanVector returns a numeric vector with the indexes of TRUE occurrences
###
.position2booleanVector <- function(thePositionVector,booleanVectorLength){
  booleanVector <- vector(mode = "logical", length = booleanVectorLength)
  booleanVector[thePositionVector] <- TRUE
  return(booleanVector)
}


### .boolean2positionVector
### - theBooleanVector is a boolean vector which TRUE elements are important
###
### .boolean2positionVector returns a numeric vector with the indexes of TRUE occurrences
###
.boolean2positionVector <- function(theBooleanVector){
  len <- length(theBooleanVector)
  positionVector <- vector(mode = "numeric", length = 0)
  for(i in 1:len){
    if(isTRUE(theBooleanVector[i])){
      positionVector[length(positionVector) + 1] <- i
    }else{}
  }
  return(positionVector)
}

### .testValueReduct
### - theCombination is a numeric vector indicating the columns (conditions Id) which are being tested as a value reduct of the Decision Indesnibility function for a given rule
### - theConditionInd is a boolean matrix representing the indiscernibility function of the conditions of a decision table (one per column) of a given rule
### - theDecisionInd is a numeric vector indicating the rule Id of the indescernibility function of the decision of a decision table for a given rule
###
### .testValueReduct returns a boolean indicating if the indiscernibility function of the condition combination tested (theCombination) is a reduct (subset of) the indescernibility function of the decision (theDecisionInd)
###
.testValueReduct <- function(theCombination,theConditionInd, theDecisionInd){
  #formerly known as  .isValueReduct
  testCondition <- theConditionInd[,theCombination]#Builds the condition combination
  if(is.vector(testCondition)){
    testCondition <- matrix(testCondition,ncol=1)#Makes sure it's dealing with matrixes
  }
  indConIntersect <- .computeIndescernibilityConditionIntersection(testCondition)
  answer <- .isSubset(indConIntersect,theDecisionInd)#Identical, isTRUE don't work when the result is a vector of TRUE
  return(answer)
}

### .isSubset
### - theVectorSubset is a numeric vector which is to be tested if it is a subset of theVectorSuperSet
### - theVectorSuperSet is a numeric vector
###
### .isSubset returns a single boolean. If the 2 set are identical, it returns true, if theVectorSubset is a 0 length vector it returns TRUE, even if theVectorSuperSet is 0 length too
###
.isSubset <- function(theVectorSubset,theVectorSuperSet){
  answer <- FALSE
  comparison <- is.element(theVectorSubset, theVectorSuperSet)
  compLen <- length(comparison)
  compSum <- sum(comparison)
  if(compSum == compLen){
    answer <- TRUE
  }else{}
  return(answer)
}

### .computeIndescernibilityConditionIntersection
### - testCondition is a boolean matrix representing on each column the position where the values of a rule(not needed for calculation) condition are duplicated
###
### .computeIndescernibilityConditionIntersection returns a numeric vector indicating the position(rule Id) where the conditions values (for a rule) intersect
###
.computeIndescernibilityConditionIntersection <- function(testCondition){
  conditionCount <- ncol(testCondition)
  indCondIntersect <- vector(mode = "numeric", length = 0)
  for(k in 1:conditionCount){
    indCond <- .boolean2positionVector(testCondition[,k])#Transform the boolean condition into position indexes
    if(k == 1){
      indCondIntersect <- indCond
    }else if(k > 1){
      indCondIntersect <- intersect(indCondIntersect,indCond)
    }else{}
  }
  return(indCondIntersect)
}


### .inflateValueReductMatrix
### - theValueReductMatrix is a numeric matrix representing the value reduct. Its columns match the columns of the ConditionReduct object but not always the columns of the decision table
### - theCondictionReductColumnIds is a numeric vector representing the column ids of the decision table which conform the condition reduct
### - theDecisionTableColumnCounter is the number of columns of the decision table (conditions + decision)
###
### .inflateValueReductMatrix returns a numeric matrix where the conditions of theValueReductMatrix match their positions on the decision table 
###
.inflateValueReductMatrix <- function(theValueReductMatrix, theCondictionReductColumnIds, theDecisionTableColumnCounter){
  
  rowCount <- nrow(theValueReductMatrix)
  colCount <- theDecisionTableColumnCounter
  cids <- theCondictionReductColumnIds
  res <- matrix(NA,nrow = rowCount, ncol = colCount)
  
  for(i in 1:rowCount){
    rule <- theValueReductMatrix[i,]
    counter <- 0
    for(j in 1:colCount){
      if(is.element(j,cids)){
        counter <- (counter + 1)
        res[i,j] <- rule[counter]
      }else{}
    }
  }
  return(res)
}