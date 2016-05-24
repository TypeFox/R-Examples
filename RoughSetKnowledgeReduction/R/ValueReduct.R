setClass (
  Class = "ValueReduct",
  representation = representation(
    conditionReduct = "ConditionReduct",
    valueReduct = "matrix"
  ),
  validity = function(object){
    #cat("~~~ ValueReduct: inspector ~~~ \n")
    
    cr <- object@conditionReduct
    dt <- getConditionReductDecisionTable(cr)
    dtMat <- getDecisionTable(dt)
    
    if(ncol(object@valueReduct) != ncol(dtMat)){
      stop ("[ValueReduct: validation] valueReduct number of columns must match the number of column Ids of the conditionReduct")
    }else{}
    return(TRUE)
  }
)


#*******************************************************
#CONSTRUCTOR

setMethod (
  f="initialize",
  signature="ValueReduct",
  definition=function(.Object,conditionReduct,valueReduct){
    #cat ("~~~~~ ValueReduct: initializator ~~~~~ \n")
    if(!missing(conditionReduct)){
      if(!missing(valueReduct)){
        .Object@conditionReduct <- conditionReduct
        .Object@valueReduct <- valueReduct
        validObject(.Object)# call of the inspector
      }else{
        .Object@valueReduct <- matrix(NA,ncol = 0, nrow = 0)
      }
    }else{
      .Object@conditionReduct <- new(Class="ConditionReduct")
    }
    return(.Object)
  }
)

#CONSTRUCTOR (USER FRIENDLY)
valueReduct <- function(theConditionReduct,theValueReduct){
  #cat ("~~~~~ ValueReduct: constructor ~~~~~ \n")
  new (Class="ValueReduct", conditionReduct = theConditionReduct, valueReduct = theValueReduct)
}


#*******************************************************
#ACCESSORS

#Returns the condition reduct object
setGeneric("getValueReductConditionReduct",function(object){standardGeneric ("getValueReductConditionReduct")})
setMethod("getValueReductConditionReduct","ValueReduct",
          function(object){
            return(object@conditionReduct)
          }
)


#Returns the value reduct as a numeric matrix
setGeneric("getValueReduct",function(object){standardGeneric ("getValueReduct")})
setMethod("getValueReduct","ValueReduct",
          function(object){
            return(object@valueReduct)
          }
)



#*******************************************************
#GENERIC METODS

#Summary

setMethod ("print","ValueReduct",
           function(x,...){
             cat("*** Class ValueReduct, method Print *** \n")
             valueReductMatrix <- x@valueReduct
             
             if(length(valueReductMatrix) != 0){
               print(formatC(valueReductMatrix),quote=FALSE)
             }else{}
             cat("******* End Print (ValueReduct) ******* \n")
           }
)


setMethod("show","ValueReduct",
          function(object){
            cat("*** Class ValueReduct, method Show *** \n")
            cat("* ValueReduct (limited to a matrix 10x10) = \n")
            valueReductMatrix <- object@valueReduct
            if(length(valueReductMatrix) != 0){
              nrowShow <- min(10,nrow(valueReductMatrix))
              ncolShow <- min(10,ncol(valueReductMatrix))
              print(formatC(valueReductMatrix[1:nrowShow,1:ncolShow]),quote=FALSE)
            }else{}
            cat("******* End Show (ValueReduct) ******* \n")
          }
)


#*******************************************************
#METHODS


### removeDuplicatedRulesVR returns a new value reduct object without duplicated rules
setGeneric (name = "removeDuplicatedRulesVR",def = function(object){standardGeneric("removeDuplicatedRulesVR")})
setMethod(
  f = "removeDuplicatedRulesVR",
  signature = "ValueReduct",
  definition = function(object){
    
    vrMat <- object@valueReduct
    vrMat <- unique(vrMat)
    cr <- object@conditionReduct
    newVr <- new(Class="ValueReduct", conditionReduct = cr, valueReduct = vrMat)
    return(newVr)
  }
)


### classifyDecisionTable returns a Decision Table object which rules have the same conditions of input DT object but the rule decisions of the Value Reduct rules where they match
setGeneric (name = "classifyDecisionTable",def = function(object,decisionTable){standardGeneric("classifyDecisionTable")})
setMethod(
  f = "classifyDecisionTable",
  signature = "ValueReduct",
  definition = function(object,decisionTable){
    
    vrMat <- object@valueReduct
    dtMat <- getDecisionTable(decisionTable)
    ruleCountVR <- nrow(vrMat)
    ruleCountDT <- nrow(dtMat)
    columnCountVR <- ncol(vrMat)
    columnCountDT <- ncol(dtMat)
    res <- new(Class="DecisionTable")
    if(columnCountVR != columnCountDT){
      print("ERROR: ValueReduct.classifyDecisionTable: The value reduct and decision table number of conditions do not match")
    }else{
      resMat <- matrix(NA,ncol = columnCountVR, nrow = 0)
      
      for(i in 1:ruleCountVR){
        ruleVR <- vrMat[i,]
        conRuleVR <- ruleVR[-length(ruleVR)]#Gets rule condition
        desRuleVR <- ruleVR[length(ruleVR)]#Gets rule decision
        flag <- FALSE
        
        for(j in 1:ruleCountDT){
          ruleDT <- dtMat[j,]
          conRuleDT <- ruleDT[-length(ruleDT)]#Gets the rule condition
          if(.ruleConditionMatch(ruleVR,ruleDT)){
            classifiedRule <- append(conRuleDT,desRuleVR,after=length(conRuleDT))
            resMat <- rbind(resMat,classifiedRule)
            flag <- TRUE#Indicates there is at least 1 value reduct rule which applies to the decision table rule
          }else{}
        }
        if(flag == FALSE){#The decision table rule cannot be classified according to the reduct value rules
          classifiedRule <- append(conRuleDT,NA,after = length(conRuleDT))
          resMat <- rbind(resMat,classifiedRule)
        }else{}
      }    
      res <- new(Class="DecisionTable",decisionTable = resMat)
    }
    return(res)
  }
)


### classifyDecisionTable returns a numeric matrix which contains the Value Reduct object representation and the support and consistency values of each rule.
setGeneric (name = "computeSupportConsistency",def = function(object,decisionTable){standardGeneric("computeSupportConsistency")})
setMethod(
  f = "computeSupportConsistency",
  signature = "ValueReduct",
  definition = function(object,decisionTable){
    
    vrMat <- object@valueReduct
    dtMat <- getDecisionTable(decisionTable)
    
    ruleCountVR <- nrow(vrMat)
    ruleCountDT <- nrow(dtMat)#Number of Rules in the decision table
    
    consistentCount <- vector(mode="numeric", length = ruleCountVR)
    inconsistentCount <- vector(mode="numeric", length = ruleCountVR)
    
    for(i in 1:ruleCountVR){
      ruleVR <- vrMat[i,]
      decRuleVR <- ruleVR[length(ruleVR)]#Gets rule decision
      ruleVRConsistentCounter <- 0 #Number of times the VR rule conditions appears in the decision table and is decision is consistent
      ruleVRNoConsistentCounter <- 0 #Number of times the VR rule conditions appears in the decision table and is decision is inconsistent
      for(j in 1:ruleCountDT){
        ruleDT <- dtMat[j,]
        decRuleDT <- ruleDT[length(ruleDT)]#Gets rule decision
        if(.ruleConditionMatch(ruleVR,ruleDT)){
          if(decRuleVR == decRuleDT){
            #Rules consistent
            ruleVRConsistentCounter <- (ruleVRConsistentCounter + 1)
          }else{
            #Rules INconsistent
            ruleVRNoConsistentCounter <- (ruleVRNoConsistentCounter + 1)
          }
        }else{}
      }
      consistentCount[i] <- ruleVRConsistentCounter
      inconsistentCount[i] <- ruleVRNoConsistentCounter
    }
    conNoConSumVector <- (consistentCount + inconsistentCount)#Number of times the VR rule appears in the DT
    support <- (conNoConSumVector/ruleCountDT)
    consistency <- (consistentCount/conNoConSumVector)
    res <- matrix(NA,ncol = 0,nrow = ruleCountVR)
    res <- cbind(res,consistentCount)
    res <- cbind(res,inconsistentCount)
    res <- cbind(res,support)
    res <- cbind(res,consistency)
    vrMat <- cbind(vrMat,res)
    return(vrMat)
  }
)





#*******************************************************
#UTIL



### .ruleConditionMatch
### - rule1 is a numeric vector representing simplify/reduced rule. It could include some NA elements
### - rule2 is a numeric vector representing a rule. It must NOT include any NA elements. This is not checked by the function
###
###
### .ruleConditionMatch returns boolean indicating if the condition of rule1 match the conditions of rule2
###
.ruleConditionMatch <- function(rule1,rule2){
  
  res <- TRUE
  condRule1 <- rule1[-length(rule1)]#Gets the rule's condition
  condRule2 <- rule2[-length(rule2)]
  comp <- (condRule1 == condRule2)
  
  for(i in 1:length(comp)){
    el <- comp[i]
    if(!is.na(el)){
      if(el == FALSE){
        res <- FALSE
        break
      }else{}
    }else{}
  }
  return(res)
}
