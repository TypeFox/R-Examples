# MicroStrategyR
# Copyright (c) 2013 MicroStrategy, Inc.
#
#' @title MicroStrategyR Package
#'
#' @description
#' Deploying R Analytics to MicroStrategy made easy.
#'
#' @details
#' The ability to deploy R Analytics to MicroStrategy combines the Best in 
#' Business Intelligence with the world's fastest growing statistical 
#' workbench.  Conceptually, think of the R Analytic as a Black Box with 
#' an R Script inside.  MicroStrategy doesn't need to know much about the 
#' R Script in the black box, it only needs to know how to pass data in 
#' and, after the script has executed in the R environment, how to consume
#' the results out.  
#' 
#' The key for MicroStrategy to execute an R Analytic implemented as an R
#' Script is capturing the R Analytic's \emph{"signature"}, a description of all
#' inputs and outputs including their number, order and nature (data type, 
#' such as numeric or string, and size, scalar or vector).  
#' 
#' The deployR utility analyses an R Script, identifies all potential 
#' variables and allows the user to specify the Analytic's signature, 
#' including the ability to configure additional features such as the 
#' locations of the script and working directory as well as controlling 
#' settings such as how nulls are handled.  
#' 
#' All this information is added to a header block at the top of the R Script.
#' The header block, comprised mostly of non-executing comment lines which are 
#' used by MicroStrategy when executing the script.  The analytic's signature 
#' represented in the header block tells MicroStrategy precisely how to 
#' execute the R Analytic.
#' 
#' Finally, in order to deploy the R analytic to MicroStrategy, the 
#' MicroStrategyR pacakge provides the metric expression of each potential output.  
#' The metric expression can be pasted into any MicroStrategy metric editor 
#' for deploying the R Analytic to MicroStrategy for execution.
#'
#' @references MicroStrategyR. \url{http://RIntegrationPack.codeplex.com/}.
#' @docType package
#' @name MicroStrategyR
#' @aliases MicroStrategyR
NULL    #Required to avoid warnings when publishing documentation

###-------------------------------------------------------------------------+
### Pop-up a message when the package is attached
.onAttach <- function(libname, pkgname) {
  if ("MicroStrategyR" %in% getOption("defaultPackages"))
    deployR()
  else {
    msg <- paste("MicroStrategyR Package\n",
                 "Copyright (c) 2013 MicroStrategy, Inc.\n",
                 "Type 'deployR()' to prepare your R Script and deploy to MicroStrategy.")
    packageStartupMessage(msg)
  }
}  #END-.onAttach

#' @title deployR Function
#'
#' @description
#' \code{deployR} is a utility that prepares an R Script for execution by MicroStrategy.  
#' It captures all the information required to properly execute the R Script into a 
#' comment block added to the beginning of the script.
#'
#' @details
#' This function launches a user interface which allows the user to open 
#' an R Script and capture its \dQuote{signature}, the nature of inputs and outputs 
#' along with other information about how the R analytic should be executed.
#' 
#' The deployR menu bar has the following functions:
#' \itemize{
#'   \item \strong{Open}:  Browse to choose the R Script you want to deploy.
#'   \item \strong{Preview}:  Preview the header block to be added to the top of the script. You'll also be given the option to save as well.
#'   \item \strong{Save}:  Saves the header block to the R Script
#'   \item \strong{Help}:  Opens the on-line help
#'   \item \strong{About}: Provides details about the deployR utility, including the version  
#'   \item \strong{Quit}:  Exits the utility
#' }
#'
#' To get started, click the \strong{Open} button at the upper left to open your R script.
#'
#' After you open your R script, the \code{deployR} utility parses the script to identify all potential variables.  
#' \itemize{
#'   \item If this script contains the MicroStrategy header block at the top then that information will be used to configure the utility and any variables not identified will appear in the \strong{Unused Variables} column.  
#'   \item On the other hand, if there is no MicroStrategy header block, the deployR utility attempts to determine whether a variable is an \strong{Input} or an \strong{Output} based on the first occurrence of that variable in the script. If the variable's first occurrence assigns it a value, it is considered an output; otherwise, it is designated as an input. 
#'   \item For new variables, the default \strong{Data Type} is \emph{Numeric} and the default \strong{Parameter Type} is \emph{Vector}.
#' }
#' The \strong{Use R Script Folder} checkbox controls whether the path used to open the script will be used, or if the R Scripts folder will be used.
#'
#' Next, modify the definition each variable as required to match the function's logic, by draging-and-droping variables to place them in the appropriate columns:
#' \itemize{
#'   \item \strong{Unused Variables} are variables that appear in the R script but are not passed between MicroStraegy and R as either inputs, outputs or parameters.  The order of unused variables does not affect the R script execution.
#'   \item \strong{Input} indicates data imported into R from MicroStrategy.  The order of inputs, from top to bottom, determines the order of arguments passed in from MicroStrategy, from left to right.
#'   \item \strong{Parameter} allows data to be passed as one of the various function parameters available for passing scalar values from MicroStratregy to R.  These function parameters include booleans, numbers, and strings. Parameters are typically used for values that are changed infrequently or values that are not determined from other metrics. Use the Parameter drop-down list to specify which parameter to use. Each parameter can only be used for one variable.  The order of parameters does not affect the R script execution.
#'   \item \strong{Output} indicates data returned from R to MicroStrategy.  If there is more than one output, the first output is considered the default output.  Beyond that, the order of outputs does not affect the signature.
#' }
#'
#' Set the \strong{Data Type} of each variable to one of the following options:
#' \itemize{
#'   \item \emph{Numeric} indicates variables that contain numbers.
#'   \item \emph{String} indicates variables that contain text.
#'   \item \emph{Default} indicates that the data type defined by MicroStrategy should be used. This setting can be used for inputs only. 
#' }
#'
#' Set the \strong{Parameter Type} of each variable to one of the following options:
#' \itemize{
#'   \item \emph{Vector} indicates a variable that contains one or more values.
#'   \item \emph{Scalar} indicates a that variable contains only one value.
#' }
#' 
#' If one or more of the variables form a repeated argument, set the \strong{Repeat Count} to match the number of repeat arguments.  This option identifies an input that can vary in quantity; such variables are known as repeated arguments because they represent a varying number of variables. The Repeat Count value specifies how many of variables are repeated, counting back from the last variable. These variables always occur at the end of the list of arguments.  These variables will appear in the Inputs column with an asterisk (*).  Examples include:
#' \itemize{
#'   \item A predictive analytical function supports one target variable Y (the dependent variable) and an indeterminate number of explanatory variables X (independent variables). Establish this configuration by setting Y as the first variable, setting X as the second variable, and setting Repeat Count to 1. The deployR R Analytic Deploer utility will understand that Y is the first argument passed into the function, followed by one or more X variables.
#'   \item A predictive analytical function supports one target variable Y (the dependent variable) and an indeterminate number of explanatory, independent variable pairs, X1 and X2, with X1 as numeric identifier for an item and X2 as its text description. By configuring Y as the first input, X1 as the second, X2 as the third, and setting the Repeat Count to 2, this defines that Y is the first argument and there can be one or more pairs of X1 and X2 variables passed into the R script.
#'   \item NOTE:  Each individual repeated variable must be of the same data type, either string or numeric.  If default is selected, the data type is determined by the first occurrence.
#' }
#'
#' The \strong{Metric Specification} section contains the following options:
#' \itemize{
#'   \item \strong{Nulls Allowed} controls whether records containing null values are to be passed in as inputs to your analytic. By default this option is selected. If unselected, all records containing null values are eliminated from the analysis
#'   \item \strong{Check Input Count} controls whether MicroStrategy ensures that the number of inputs to the metric matches exactly with the number of inputs specified in the function's signature. By default, the option is selected. If selected and the number of inputs is different, a specific warning message will be returned to the user.  If unchecked and the number of inputs is different, the script execution will attempt to proceed.
#'   \item \strong{Enable Sort By} controls the sorting of records before being passed to R. By default, the option is selected. If this option is selected, the first input must be a vector, since the default behavior sorts records in ascending order by the first input. To specify a particular sorting criterion, you can specify the Sort By value in the field below the check box. If this option is cleared, the order of records passed into R is determined by MicroStrategy automatically.
#'   \item \strong{Specify Working Directory} allows you to specify a working directory for R. By default, this option is cleared. To select a working directory select the check box and specify a working directory in the field below the check box. If this option is selected, MicroStrategy does not alter R's working directory, which is otherwise determined by R.
#'   \item \strong{Output Variable} allows the user to control which output is returned to MicroStrategy. The first output is the default and does not need to be specified in the metric expression.
#' }
#' 
#' After you have configured the variables and specified the metric options, you can save the analytic's signature to the R script by clicking the \strong{Save} button.  Alternatively, you can use the \strong{Preview} button to review the changes before saving.  
#'
#' Finally, after configuring your variables, saving the signature to your script and specifying the metric options, the R Analytic Deployer provides you with the metric expression for your analytic.  The \strong{Metric Expression} section at the bottom-right of the dialog shows the metric expression that defines how MicroStrategy interacts with your function.  Click the \strong{Copy to Clipboard} button and you're ready to paste this metric expression into any MicroStrategy Metric Editor (Desktop, Web and Visual Insight).
#'
#' @references MicroStrategyR. \url{http://RIntegrationPack.codeplex.com/}.
#' @name deployR
#' @aliases deployR
#' @keywords microstrategy
#' @export
#' @examples
#' \dontrun{
#' deployR()}
deployR <- function() {
  VERSION <- "1.0-1"
  COPYRIGHT <- "Copyright (c) 2013 MicroStrategy, Inc."  
  MSTR_WORKINGDIR <- "mstr.WorkingDir"
  #MSTR_EXFLAGIN   <- "mstrExFlag"
  MSTR_EXFLAG     <- "mstr.ExFlag"
  MSTR_ERRMSG     <- "mstr.ErrMsg"
  MSTR_INPUTNAMES <- "mstr.InputNames"
#  EXCLUDED_NAMES <- c(".GlobalEnv", "...", ".", MSTR_EXFLAGIN, MSTR_EXFLAG, MSTR_ERRMSG, MSTR_WORKINGDIR, MSTR_INPUTNAMES)
  EXCLUDED_NAMES <- c(".GlobalEnv", "...", ".", MSTR_EXFLAG, MSTR_ERRMSG, MSTR_WORKINGDIR, MSTR_INPUTNAMES)
  DROP_HANDLERS   <- "dropHandlers"
  WINDOW_SIZE <- c(1000, 600)    
  R_3 <- (as.numeric(R.Version()$major)>=3)   #Flag if executing in R version 3.0.0 or later

  ###-------------------------------------------------------------------------+
  ###Track when a parameter is used or un-used
  TrackParams <- function(pName="", isUsed=TRUE) {
    if(DEBUG<-FALSE) print(paste0("Function enter: TrackParams=", pName, "<-", isUsed))
    if(nchar(pName)>0) {
      pIndex <- (match(tolower(pName), tolower(envRAD$paramList)))    #Get the index of this parameter name, if any
      if(is.na(pIndex)) {                                              #If that's not a valid parameter, return FALSE
        return(-1)
      } else {  
        if(envRAD$paramUsed[pIndex]&&isUsed) {
          return(-2)   #This parameter is already used
        } else {
          envRAD$paramUsed[pIndex] <- isUsed
        }
      }
    }
    return(pIndex)  
  }  #END-TrackParams

  
  ###-------------------------------------------------------------------------+
  ### Show the help menu
  ShowHelp <- function(h, ...) {
    w <- gwindow("deployR Help", visible=FALSE)
    size(w) <- WINDOW_SIZE  
    g <- ggroup(horizontal=FALSE, container=w)
    #g1 <- ggroup(container=g)
    #addSpring(g1)
    #glabel("Help on:", container=g1)
    #e <- gedit("", container=g1)
    helpWidget <- ghelp(container=g, expand=TRUE)
    #addHandlerChanged(e, handler=function(h, ...) {
    #  add(helpWidget, svalue(h$obj))
    #})
    ## add others
    #add(helpWidget,"base:::mean")
    add(helpWidget, list(topic="MicroStrategyR", package="MicroStrategyR"))
    add(helpWidget, list(topic="deployR", package="MicroStrategyR"))
#    add(helpWidget, list(topic="MenuBar", package="MicroStrategyR"))
    #add(helpWidget, "boxplot")
    visible(w) <- TRUE
  }  #END-ShowHelp
  
  ###-------------------------------------------------------------------------+
  ### Get the parameter from a MicroStrategy variable comment line
  GetParam <- function(x) {
    if(DEBUG<-FALSE) print(paste0("Function enter: GetParam=", x))
    pName <- ""                                         #Default result to return
    tokens <- unlist(strsplit(x,"[[:space:]]+"))        #Break up string into tokens:  For a parameter there should only be three, in the form "varName -p/-parameter pName"
    if(length(tokens)>=3) {                             #Must have at least 3 tokens for a parameter
      if(tolower(substr(tokens[2],1,2))=="-p") {        #If the second token starts with "-p"
        pName <- tokens[3]                              #pName is the third token
        i <- TrackParams(pName, TRUE)                   #Track and get index if valid, error code if invalid
        if(i>0) {                                       #If that's a valid parameter that's not already used, return it
          pName <- envRAD$paramList[i]                  #  Get the parameter's name
          return(pName)                                 #  Return it
        } else {                                        #If that's not a valid parameter or its already used, show dialog
          wParam <- gbasicdialog("Invalid Parameter", horizontal=FALSE, handler=function(h, ...){
            pName <- svalue(lParam)                     #  User has selected to replace with an existing parameter
            i <- TrackParams(pName, TRUE)               #  Track that this parameter is used
            envRAD$temp <- envRAD$paramList[i]
          })
          glabel(paste0(ifelse(i==-1, "This invalid parameter was found in the script.", "This parameter is already used.")), container=wParam, anchor=c(-1,0))
          glabel(paste0("\"", pName, "\""), container=wParam, anchor=c(-1,0))
          glabel(paste0(ifelse(i==-1, "Please select the correct parameter for it:", "Please select another parameter for it:")), container=wParam, anchor=c(-1,0))
          lParam <- gdroplist(envRAD$paramList[!envRAD$paramUsed], container=wParam, editable=FALSE, selected=1)
          envRAD$temp=""
          visible(wParam, set=TRUE) # show dialog
          return(envRAD$temp)
        }
      } else {
        return("")
      }
    } else {
      return("")
    }
  } #END-GetParam
  
  SetFunctionID <- function(thisOutput=NA) {    
    if(DEBUG<-FALSE) print("Function enter: SetFunctionID")
    inputCount <- length(envRAD$listIns)                                           #Count the number of inputs
    v <- match(envRAD$listIns, envRAD$dfVar[, 1])                                  #Get the input indexes, in ordef
    inputVectorCount <- sum(envRAD$dfVar[v, 4]=="Vector")                          #Count the number of inputs that are vectors
    bVectorInput <- (inputVectorCount>0)                                           #Flag is any inputs are vectors
    if(is.na(thisOutput)) {                                                        #If the output to use is NA
      if (is.null(svalue(envRAD$lstOut))) {                                        #  Then if the lstOut value is null
        thisOutput <- 1                                                            #    Use the first Output
      } else {                                                                     #  Else if the lstOut value is not null 
        thisOutput <- svalue(envRAD$lstOut, index=TRUE)                            #    Use the selected Output
      }
    }
    v <- match(envRAD$listOuts, envRAD$dfVar[, 1])                                 #Get the output indexes, in order
    outputVectorCount <- sum(envRAD$dfVar[v, 4]=="Vector")                         #Count the number of outputs that are vectors
    outputCount <- length(envRAD$listOuts)                                         #Count the number of outputs
    if ((inputCount==0)||(outputCount==0)) {                                       #Make sure we have at least one input and at least one output
      envRAD$functionID <- 7                                                       #Need to have at least one input or output, ERROR-7
    } else {                                                                       #Ok, We have at least one input and at least one output
      bVectorOutput <- (envRAD$dfVar[v[thisOutput], 4]=="Vector")                  #Flag if the selected output is a vector
      envRAD$functionID <- ifelse(!bVectorInput&&!bVectorOutput, 1,                # VectorInputs VectorOutput  functionID
                                  ifelse(bVectorInput&&bVectorOutput, 2,           #      0           FALSE     Simple-1           
                                      ifelse(bVectorInput&&!bVectorOutput,3, 6)))  #     >=1          TRUE      Relative-2
      #      }                                                                     #     >=1          FALSE     Agg-3
    }                                                                              #      0           TRUE      ERROR-6
    if((envRAD$functionID<2)||(envRAD$functionID>3)) {                             #Handle SortBy Option based on functionID
      enabled(envRAD$chkSort) <- FALSE                                             #Disable SortBy Option if not Relative-2 or Agg-3
    } else {                                                                       #functionID must be Relative-2 or Agg-3
      enabled(envRAD$chkSort) <- TRUE                                              #Enable SortBy 0ption since functionID is Relative-2 or Agg-3 
      if(svalue(envRAD$chkSort)) {                                                 #Determine if SortBy Option is Selected
        inputsParamType <- envRAD$dfVar[(envRAD$dfVar[, 2]=="Input"), 4]           #SortBy is selected, get scalar/vector type of the first input 
        if(envRAD$dfVar[envRAD$listIns[1], 4]=="Vector") {                                         #Determine if the first input is a vector (a requirement for SortBy)
          envRAD$functionID <- envRAD$functionID + 2                               #First input is a vector so use SortBy Function: Relative-2 --> RelativeS-4, Agg-3 --> AggS-5
        } else {                                                                   #Otherwise, first input is not a vector but SortBy is selected
          envRAD$functionID <- 8                                                   #Need to have vector first input to select Sort By --> ERROR-8
        }
      }
    }
  }  #END-SetFunctionID
  
  
  ###-------------------------------------------------------------------------+
  ### Concatenate two strings with a comma, if the first string is not empty
  CommaCat <- function(A, B) {
    return(paste0(A, ifelse((nchar(A)>0), ", ", ""), B))
  }  #END-CommaCat
  
  ### Get string for variables
  GetVarInfo <- function(totalCount, vectorCount) {
    return(paste0(ifelse(totalCount==1,
                         ifelse(vectorCount==0," (scalar)", " (vector)"),
                         ifelse(vectorCount==0, "s (all scalar)", 
                                ifelse(totalCount==vectorCount, "s (all vector)", 
                                       paste0("s (",vectorCount, " vector, ", totalCount-vectorCount, " scalar)"))))))
  }  #END-GetVarInfo
  
  SetOutputs <- function(...) {
    if(!envRAD$allowUpdates) return()                                    #Don't process actions while updates are NOT allowed
    if(DEBUG<-FALSE) print("Function enter: SetOutputs")    
    envRAD$allowUpdates <- FALSE                                         #Temporarily block updates
    selOut <- svalue(envRAD$lstOut)                                      #Remember the currently selected output
    envRAD$lstOut[] <- envRAD$listOuts                                    #Update the gdroplist control with the current outputs
    update(envRAD$lstOut)                                                #Update the control
    if(length(envRAD$lstOut)==0) {                                       #Check if there are no outputs
      enabled(envRAD$lstOut) <- FALSE                                    #When there are no outputs, disable the gdroplist control 
      envRAD$lstOut[] <- c("<-- No Outputs! -->")                        #Add one output so the user will know that there's no outputs
      svalue(envRAD$lstOut, index=TRUE) <- 1                             #Select the first (and only) output in the list
    } else {                                                             #There's at least one output
      enabled(envRAD$lstOut) <- TRUE                                     #Make sure the gdroplist control is enabled when there's at least one output
      envRAD$lstOut[1] <- paste0(envRAD$lstOut[1]," (default)")          #Indicate that the first output is the default
      if(is.null(selOut)) {                                              #If there was no selected output to remember
        svalue(envRAD$lstOut, index=TRUE) <- 1                           #Then select the first output
      } else {                                                           #Otherwise, select the remembered output or, if that one is no longer in the list, select the first output
        svalue(envRAD$lstOut, index=TRUE) <- match(selOut, envRAD$lstOut[], nomatch=1)
      }
    }    
    envRAD$allowUpdates <- TRUE                                          #Re-enable updates
  }  #END-SetOutputs
  
  ###-------------------------------------------------------------------------+
  ### Get parameter value from Metric Expression
  GetParamValue <- function(metExp, pName) {
    if(DEBUG<-FALSE) print(paste0("Function enter: GetParamValue, MetExp=", metExp,", pName=", pName))
    i <- grep(pName, envRAD$paramList)
    if(i<=18) {     #This is a boolean or number
      val <- unlist(strsplit(unlist(strsplit(metExp, paste0(pName,"="), fixed=TRUE))[2], ","))[1]
      val <- unlist(strsplit(val, ">"))[1]
    } else {             #This is a string
      val <- unlist(strsplit(unlist(strsplit(metExp, paste0(pName,"=\""), fixed=TRUE))[2],"\""))[1]
    }
    return(ifelse(is.na(val),"",val))
  }  #END-GetParamValue
  
  
  UpdateExpression <- function(thisOutput=NA) {
    if(!envRAD$allowUpdates) return()    #Don't process actions while updates are NOT allowed
    if(DEBUG<-FALSE) print(paste0("Function enter: UpdateExpression=", thisOutput))    
    SetFunctionID(thisOutput)  #Set the function type
    if ((envRAD$functionID>5)&&(is.na(thisOutput))) {
      svalue(envRAD$txtMetricExp) <- paste0("<-- Function Error: ", envRAD$funcDesc[envRAD$functionID], " -->")
    } else {
      
      #Set the Metric Arguments
      n <- length(envRAD$listIns)
      envRAD$rVarNames <-envRAD$listIns[1]
      if(n>1) {
        for(i in 2:n) envRAD$rVarNames<- paste0(envRAD$rVarNames, ", ", envRAD$listIns[i])
      }
      
      #Set the Function Parameters 
      rScriptName <- ifelse(svalue(envRAD$chkScr), basename(svalue(envRAD$edtScr)), svalue(envRAD$edtScr))
      sParams <- paste0("_RScriptFile=\"", rScriptName, "\", _InputNames=\"", envRAD$rVarNames, "\"")      
      if(!svalue(envRAD$chkCIC))              sParams <- CommaCat(sParams, "_CheckInputCount=False")
      if(!svalue(envRAD$chkNul))              sParams <- CommaCat(sParams, "_NullsAllowed=False")
      if(is.na(thisOutput)) {
        if (!is.null(svalue(envRAD$lstOut))) {
          if(svalue(envRAD$lstOut, index=TRUE)>1) sParams <- CommaCat(sParams, paste0("_OutputVar=\"", svalue(envRAD$lstOut), "\""))
        }
      } else {
        if(thisOutput>1) sParams <- CommaCat(sParams, paste0("[_OutputVar]=\"", envRAD$lstOut[thisOutput], "\""))        
      }
      if(svalue(envRAD$chkDir)) sParams <- CommaCat(sParams, paste0("_WorkingDir=\"", svalue(envRAD$edtDir), "\""))
      if(svalue(envRAD$chkSort)&&(svalue(envRAD$edtSort)!=envRAD$SORT_BY_DEFAULT)&&(nchar(svalue(envRAD$edtSort))>0)) { 
        sParams <- CommaCat(sParams, paste0("SortBy=(", svalue(envRAD$edtSort), ")"))
      }
      vParams <- envRAD$dfVar[envRAD$listPars, ]
      vParams <- as.data.frame(vParams[nchar(vParams[, 6])>0, c(5, 6)])
      n <- nrow(vParams)
      if(n>0) {
        isStr <- (substr(vParams[, 1], 1, 1)=="S")
        p1 <- ifelse((grep(TRUE, isStr)>0), paste0(", ", vParams[isStr, 1],"=\"", vParams[isStr, 2], "\""), "")
        p2 <- ifelse((grep(FALSE, isStr)>0), paste0(", ", vParams[!isStr, 1],"=", vParams[!isStr, 2]))
        p <- append(p1, p2)
        for(i in 1:n) {sParams <- paste0(sParams, p[i])}
        if(substr(sParams, 1, 1) == ",") sParams <- substr(sParams, 3, nchar(sParams)-2)
      }         
      if(nchar(sParams)>0)      sParams <- paste0("<", sParams, ">")
      
      if(is.na(thisOutput)) {
        svalue(envRAD$txtMetricExp) <- paste0(envRAD$funcName[envRAD$functionID], sParams, paste0("(", envRAD$rVarNames, ")"))      
      } else {
        return(paste0(envRAD$funcName[envRAD$functionID], sParams, paste0("(", envRAD$rVarNames, ")"))) 
      }
    }
    
  } #END-UpdateExpression

  GetCommentBlock <- function() {
    if(DEBUG<-FALSE) print("Function enter: GetCommentBlock")
    
    commentBlock <- ""   #Create the comment block
    
    #Inputs first:
    rptCount <- svalue(envRAD$RptCt)
    inputCount <- 0
    numInputs <- length(envRAD$listIns)
    for (v in envRAD$listIns) {
      i = match(v, envRAD$dfVar[, 1])
      inputCount <- inputCount + 1
      commentBlock <- paste0(commentBlock, "#RVAR ", envRAD$dfVar[i, 1], " -", tolower(envRAD$dfVar[i, 2]), 
                             ifelse((envRAD$dfVar[i, 3]=="Default"), "", paste0(" -", tolower(envRAD$dfVar[i, 3]))), 
                             " -", tolower(envRAD$dfVar[i, 4]), 
                             ifelse(inputCount>(numInputs-rptCount)," -repeat", ""), "\n")
    }
    
    numParams <- length(envRAD$listPars)
    if(numParams>0) {
      commentBlock <- paste0(commentBlock, "#\n")
      for (v in envRAD$listPars) {
        i = match(v, envRAD$dfVar[, 1])
        commentBlock <- paste0(commentBlock, "#RVAR ", envRAD$dfVar[i, 1], " -parameter ", envRAD$dfVar[i, 5], "\n")        
      }
    }
    
    outputCount <- 0
    envRAD$functionError <- 0  #Assume that there's no function error
    numOutputs <- length(envRAD$listOuts)
    commentBlock <- paste0(commentBlock, "#\n")
    for (v in envRAD$listOuts) {
      i = match(v, envRAD$dfVar[, 1])
      outputCount <- outputCount + 1
      commentBlock <- paste0(commentBlock, "#RVAR ", envRAD$dfVar[i, 1], " -", tolower(envRAD$dfVar[i, 2]), 
                             ifelse((envRAD$dfVar[i, 3]=="Default"), "Numeric", paste0(" -", tolower(envRAD$dfVar[i, 3]))), 
                             " -", tolower(envRAD$dfVar[i, 4]), 
                             paste0(envRAD$TAG_METRIC_EXP, UpdateExpression(outputCount)), "\n")
      if(envRAD$functionID>5) envRAD$functionError <- envRAD$functionID   #Flag that there's a function error
    }
    
    MSTR_STATEMENT <- paste0("if(exists(\"", MSTR_WORKINGDIR, "\")) setwd(", MSTR_WORKINGDIR, ")")
    MSTR_COMMENT   <- "#Working Directory if executed by MicroStrategy"
    commentBlock <- paste0(commentBlock, MSTR_STATEMENT, "  ", MSTR_COMMENT, "\n")    

    #Return comment block
    return(paste0("#MICROSTRATEGY_BEGIN\n#\n", commentBlock, "#\n#MICROSTRATEGY_END\n"))
    
  }  #END-GetCommentBlock
  
  SaveToFile <- function(...) {
    if(DEBUG<-FALSE) print("Function enter: SaveToFile")

    txtCommentBlock <- GetCommentBlock()  #Get the Comment Block
    if(!OkToSave()) return()
    
    sErrMsg <- tryCatch({  
      saveFile <- file(svalue(envRAD$edtScr), "w")
      cat(txtCommentBlock, file=saveFile, sep="\n")
      cat(envRAD$g.script, file=saveFile, sep="\n")
      close(saveFile)
      envRAD$saved <- TRUE
      gmessage("\nScript saved successfully!\n", icon="info", title="Save")
      sErrMsg <- ""                                           #If we made it here, then there's no errors to report  
    }, error = function(err) {  
      sErrMsg <- geterrmessage()                            #Report error that was caught
    })
    if(nchar(sErrMsg)>0) {
      gmessage(paste0("The script was not able to be saved due to the following error: ","\n\n", sErrMsg), 
               title="Script Not Saved", icon="error")
    }
        
  }  #END-SaveToFile

  OkToSave <- function() {
    if(DEBUG<-FALSE) print("Function enter: OkToSave")
    if(envRAD$functionID < 0) return(FALSE)  #No function ID yet
    if(envRAD$functionError > 5)  {
      gmessage(paste0("Please fix this error before saving:\n\n", envRAD$funcDesc[envRAD$functionError]), title="Cannot Save Due to Error", icon="error")
      return(FALSE)
#    } else {
#      if(envRAD$functionError) {
#        gmessage("One or more outputs (not currently selected) has an error, please fix this problem before saving\n\n", title="Cannot Save Due to Error", icon="error")
#        return(FALSE)
#      }
    }
    return(TRUE)  
  }  #END-OkToSave
  
  ###-------------------------------------------------------------------------+
  ### Function to save the R Script with the Function Signature comment block at the top
  ScriptSave <- function(...) {
    if(DEBUG<-FALSE) print("Function enter: ScriptSave")
    
    txtCommentBlock <- GetCommentBlock()  #Get the Comment Block
    if(!OkToSave()) return()
    
    #Create and show the Save Dialog
    wSave <- gwindow("Save")
    size(wSave) <- WINDOW_SIZE
    gSave <- ggroup(container = wSave, horizontal=FALSE)
    fCB <- gframe("This header block will be added to the top of your script:", container=gSave, anchor=c(-1,0), expand=TRUE)
    tCB <- gtext(txtCommentBlock, container=fCB, expand=TRUE, font.attr=list(family="monospace"))
    fSc <- gexpandgroup("R Script:", container=gSave, expand=FALSE)
    tSc <- gtext(paste0(envRAD$g.script, collapse="\n"), container=fSc, expand=TRUE, font.attr=list(family="monospace"), height=500)
    #addSpring(gSave)
    gButtons <- ggroup(container = gSave, horizontal=TRUE)    ## A group to organize the buttons
    addSpring(gButtons)                                       ## Push buttons to right
    gbutton("Save", container=gButtons, 
            handler = function(h, ...) {
              SaveToFile()
              dispose(wSave)
            })
    gbutton("Cancel", container=gButtons, handler = function(h, ...) dispose(wSave))    
  } #END-ScriptSave
  
  UpdateStatusBar <- function() {
    txt <- svalue(envRAD$statusBar)
    if (nchar(txt)>111) {
      svalue(envRAD$statusBar) <- ">>> Loading"
    } else {
      svalue(envRAD$statusBar) <- paste0(txt, ".")
    }
  }  #END-UpdateStatusBar
  
  ###-------------------------------------------------------------------------+
  ### Set variables in the Environment
  SetEnvVars <- function() {
    envRAD$paramUsed[] <- FALSE
    envRAD$workDr <- NA
    envRAD$sortBy <- NA
    envRAD$args <- NA
    envRAD$functionID <- -1
    envRAD$allowUpdates <- FALSE
    envRAD$listVars <- character(0)
    envRAD$listIns <- character(0) 
    envRAD$listPars <- character(0)
    envRAD$listOuts <- character(0)
  }  #END-SetEnvVars
  
  ###-------------------------------------------------------------------------+
  ### Function to Open a script and parse it to find the potential variables
  ScriptOpen <- function(...) {
    if(DEBUG<-FALSE) print("Function enter: ScriptOpen")
    
    #if (envRAD$loaded) return()
    
    FILE_CHOOSE <- TRUE         #TRUE = Users chooses file; FALSE = Static file name is used
    
    # Get R Script
    if(FILE_CHOOSE)(g.file <- file.choose())                 #For normal use, let's user pick the R Script to open
    else (g.file <- "E:\\Demo\\ARIMA\\ARIMA.R")              #For testing only, opens this specific file without prompting the user

    add(envRAD$winRAD, envRAD$statusBar)
    UpdateStatusBar()
    
    if(!exists("g.file")) return(FALSE)                      #Return if the file doesn't exist

    envRAD$g.file <- g.file                                #Persist the file name      
    envRAD$workDr <- gsub("/", "\\", dirname(g.file), fixed=TRUE)
    tryCatch({  
      scriptErrMsg <- ""                                   #Clear script error message
      envRAD$g.script <- scan(envRAD$g.file, 
                              what=character(), 
                              sep="\n",  
                              quiet=!(DEBUG<-FALSE))       #Persist the script for use when saving
    #Parse the R script
#      if(R_3) {
    scriptParsed <- parse(text=envRAD$g.script)       #Parse the R script using built-in parser
#      } else {
#        scriptParsed <- parser(text=envRAD$g.script)      #Parse the R script using parser package
#      }
    }, error = function(err) {  
      scriptErrMsg <- geterrmessage()                     #Report error that was caught
    })
    
    if(nchar(scriptErrMsg)>0) {                           #There was an error
      gmessage(paste0("The script was not able to be processed due to the following error: ","\n\n", scriptErrMsg), 
               title="Error with Script", icon="error")
      return(FALSE)                                        #Return
    }             
      
    #Process the R Script
        
    envRAD$allowUpdates <- FALSE                         #Don't allow updates while the script is being processed        
    SetEnvVars()
    UpdateStatusBar()

    #Tokenize Script
    #tokens <- attr(scriptParsed, "data")                 #The script tokens are in this data.frame
    if(R_3) {
      tokens <- getParseData(scriptParsed)               #Create data frame with script tokens
      col.token <- "token"  
    } else {
      tokens <- attr(scriptParsed, "data")               #The script tokens are in this data.frame
      col.token <- "token.desc"
    }
    #if(DEBUG<-FALSE) TOKENS <<- tokens                   #For debugging, make a global copy of the TOKENS
    
    #Set columns 
    vars  <- character(0)
    dir   <- character(0)
    dtype <- character(0)
    ptype <- character(0)
    param <- character(0)
    pVals <- character(0)

    #Handle any existing comment blocks

    #Gather all the old style comment lines to them
    mstrVarStart1 <- as.numeric(tokens[grep("[[:space:]]*##*[[:space:]]*MSTR_VAR_START[[:space:]]*", toupper(tokens$text)), 1])        #Detect the start -- old style
    mstrVarEnd1  <- as.numeric(tokens[grep("[[:space:]]*##*[[:space:]]*MSTR_VAR_END[[:space:]]*", toupper(tokens$text)), 1])           #Detect the end -- old style        
    
    mstrVarStart <- as.numeric(tokens[grep("[[:space:]]*##*[[:space:]]*MICROSTRATEGY_BEGIN[[:space:]]*", toupper(tokens$text)), 1])   #Detect the start -- new style
    mstrVarEnd  <- as.numeric(tokens[grep("[[:space:]]*##*[[:space:]]*MICROSTRATEGY_END[[:space:]]*", toupper(tokens$text)), 1])      #Detect the end -- new style
    #mstrVarStart <- c(mstrVarStart1, mstrVarStart2)                                                                        #Combine old and new starts
    #mstrVarEnd <- c(mstrVarEnd1, mstrVarEnd2)                                                                              #Combine old and new end
    
    rptCount <- 0                                   #Set default for Repeat Count
    pValues <- rep("", length(envRAD$paramList))    #Set the parameter values to empty strings
    names(pValues) <- envRAD$paramList              #name the parameter values to match their parameter
    varLines <- numeric(0)                          #Start with no lines to remove
    UpdateStatusBar()                               #Add a tick
    
    if(length(mstrVarStart)!=length(mstrVarEnd)) {  #Handle situation where there are un-matched start and end points
      gmessage(paste0("Number of #MICROSTRATEGY_BEGIN comments (", length(mstrVarStart), 
                      ") does not equal the number of #MICROSTRATEGY_END comments (", length(mstrVarEnd), 
                      "). Please fix and try again."), 
               title="MicroStrategy Variables Comment Block Problems", 
               icon="error")
      return(FALSE)
    }
    
    priorVarsExist <- FALSE      #Boolean flag default -- is that there is no header block
    if(length(mstrVarStart)>0) {  
      #Handle existing comments

      #col.token <- ifelse(R_3, "token", "token.desc")  #MOVED UP UNDER "#Tokenize Script"
      
      mstrVars <- numeric(0)        
      for(i in 1:length(mstrVarStart)) mstrVars <- append(mstrVars, seq(mstrVarStart[i]+1, mstrVarEnd[i]-1))       #Gather all the variables inside the MicroStrategy Comment Block(s)
      mstrVars <- mstrVars[(tokens[mstrVars, col.token]=="COMMENT")]                                            #Get all the Comment lines
      mstrVars <- mstrVars[grep("[[:space:]]*#[[:space:]]*RVAR[[:space:]][[:space:]]*", tokens[mstrVars, "text"], ignore.case=TRUE)]   #Containing #RVAR (whitespace & case tolerant)

      if(length(mstrVars)>0) {                                                                                       #There are variables
        varSpec <- sub("[[:space:]]*#*[[:space:]]*RVAR[[:space:]]*","", tokens[mstrVars, "text"], ignore.case=TRUE)  #Trim off the #RVAR part, leaving the variable spec text
        varText <- sapply(varSpec, function(x) (unlist(strsplit(x, envRAD$TAG_METRIC_EXP)))[1])                      #Split outputs, left is the varText
        metExp <- sapply(varSpec, function(x) (unlist(strsplit(x, envRAD$TAG_METRIC_EXP)))[2])                       #Split ouputs, right is the metric expression
        metExp <- metExp[!is.na(metExp)]                                                                             #Keep just the metric expressions that exist
        if(length(metExp)>0) {                                                                                       #See if we have a metric expression, and if so, get the Working Directory and SortBy values
          envRAD$workDr <- unlist(strsplit(unlist(strsplit(metExp[1], "_WorkingDir=\"", fixed=TRUE))[2],"\""))[1]  #Check the first metric expression for the Working Directory, if any (just need to check the first since they should all be the same, except for the _output parameter)
          envRAD$sortBy <- unlist(strsplit(unlist(strsplit(metExp[1], "SortBy=(", fixed=TRUE))[2],")"))[1]           #Check the first metric expression for the SortBy, if any
          envRAD$args <- unlist(strsplit(unlist(strsplit(metExp[1], ">(", fixed=TRUE))[2],")"))[1]                   #Check the first metric expression for the function arguments, if any    
          pValues <- sapply(envRAD$paramList, function(x) (GetParamValue(metExp[1], x)))                             #Get the values for any parameters
          names(pValues) <- envRAD$paramList
          UpdateStatusBar()
        }

        #Gather the nature of the variables
        bDis <- sapply(varText, function(x) (length(grep("-d", x, ignore.case=TRUE))>0))    #Flag Disableds
        bOut <- sapply(varText, function(x) (length(grep("-o", x, ignore.case=TRUE))>0))    #Flag Outputs
        bPar <- sapply(varText, function(x) (length(grep("-p", x, ignore.case=TRUE))>0))    #Flag Parameters
        bNum <- sapply(varText, function(x) (length(grep("-n", x, ignore.case=TRUE))>0))    #Flag Numerics
        bStr <- sapply(varText, function(x) (length(grep("-st", x, ignore.case=TRUE))>0))   #Flag Strings
        bRep <- sapply(varText, function(x) (length(grep("-r", x, ignore.case=TRUE))>0))    #Flag Repeats
        bSca <- sapply(varText, function(x) (length(grep("-s[^tr]|-s$", x, ignore.case=TRUE))>0)) #Flag Scalars - Improvement from Li Zhang -- thank you Li!
        UpdateStatusBar()
          
        #Handle Repeat Count
        bIn <- !bOut & !bPar                                                                #bIn array wiil have TRUE for any inputs
        firstRpt <- match(TRUE, bRep)                                                       #Find the first repeat, if any
        if(!is.na(firstRpt)) {
          rptCount <- (length(bIn[bIn==TRUE])-firstRpt[1])+1                                #Set repeat count
        } else {
          rptCount <- 0
        }
        UpdateStatusBar()
          
        #Set the spec of the variables
        vars <-sapply(varText, function(x)(unlist(strsplit(sub("[[:space:]]*#[[:space:]]*","",x),"[[:space:]]+")))[1])  #Get variable names -- Improvement from Li Zhang -- thank you Li!                      
        dir <- ifelse(bOut, "Output", ifelse(bPar, "Parameter", ifelse(bDis, "Unused", "Input")))                     #Set direction
        dtype <- ifelse(bPar, "Numeric", ifelse(bStr, "String", ifelse(bNum, "Numeric", "Default")))                           #Set data type
        ptype <- ifelse(bPar, "Vector", ifelse(bSca, "Scalar", "Vector"))                   #Set parameter type
        param <- sapply(varText, function(x) (GetParam(x)))                                 #Get Parameters  
        pVals <- sapply(param, function(x) ifelse(is.na(pValues[x]), "", pValues[x]))       #Set the parameter values, if any
        UpdateStatusBar()
        #Gather all the comment lines to remove and remove them
        mstrVarStart <- c(mstrVarStart1, mstrVarStart)                                                                       #Combine old and new starts
        mstrVarEnd <- c(mstrVarEnd1, mstrVarEnd)                                                                              #Combine old and new end
        for(i in 1:length(mstrVarEnd)) {                                                    #For each header block end line
          while(envRAD$g.script[mstrVarEnd[i]+1]=="") {                                     #Look to see if the next line is an empty line
            mstrVarEnd[i] <- mstrVarEnd[i]+1                                                #If so, include the empty line as part of the header block (so it can be removed)
          }
        }
        for(i in 1:length(mstrVarStart))                                                    #For each header block
          varLines <- append(varLines, seq(mstrVarStart[i], mstrVarEnd[i]))                 #  Collect the lines to remove
        envRAD$g.script <- envRAD$g.script[-varLines]                                       #Remove the comment lines
        
        tokens <- tokens[!(tokens$line1 %in% varLines), ]                                   #Remove all tokens from the MicroStrategy Header Block
        rownames(tokens) <- seq(1:nrow(tokens))                                             #Rename the tokens after we've deleted the header rows, if any
        priorVarsExist <- TRUE                                                              #Boolean flag if there was a header block
      }  #Processed all the existing variables 
    
    } #Done with header block(s)!
    UpdateStatusBar()
    
    #Next, Handle new variables from the R Script
    symbolIds <- tokens[(tokens[, col.token]=="SYMBOL"), "id"]                              #Get the Id for the SYMBOLs in the script, these are potential variables
    
    #Make sure existing symbols are not "member variables" in the form "env$member"
    keptIds <- symbolIds                                                                    #Make a copy of the list of symbols, these we'll keep
    for (i in 1:length(symbolIds)) {                                                        #For every symbol
      memberId <- symbolIds[i]-2                                                            #  Get the member indicator, which would be 2 less then the symbol id
      if(memberId %in% tokens[,"id"]) {                                                     #    If the member Id exists  
        if(tokens[tokens[,"id"]==memberId,"text"]=="$") {                                   #      Then check to see if a member is indicated by "$"
          keptIds <- keptIds[!(keptIds==symbolIds[i])]                                      #        If so, don't keep this symbol
        }
      }
    }
#SID <<- symbolIds
#KID <<- keptIds    
    symbols <- tokens[tokens[, "id"] %in% keptIds, ]                                        #Get a data frame "symbols" with all the non-member variables kept from the step before
#S1 <<- symbols
    symbols <- symbols[!(symbols[, "text"] %in% EXCLUDED_NAMES), ]                          #Exclude certain variables that used only by MicroStrategy, never as inputs, outputs or parameters
#S2 <<- symbols
    symbolIds <- symbols[, "id"]                                                            #Update the list of symbol Ids that we're using
    
#    vars <- append(vars, tokens[symbols, "text"])                                           #Add the variable names
    vars <- append(vars, symbols[, "text"])                                                 #Add the variable names
#V <<- vars
    UpdateStatusBar()
    
    #Determine if the new variables are inputs or outputs, or unused                        
    if(priorVarsExist) {                                                                    #If prior variables exist
      UpdateStatusBar()
      dir <- append(dir, rep("Unused", length(symbolIds)))                                  #  then set any new variables as unused (initially)
    } else {                                                                                #Otherwise, determine the direction of new variables from the script
      bOutputVar <- rep(FALSE, length(symbolIds))                                           #Initialize all output vars to FALSE
#      assignAdderLE <- ifelse(R_3, 1, 2)                                                    #Adder for the id difference between the symbol and any potential assignment operator for LEFT_ASSIGN or EQ_ASSIGN
#      assignAdderR <- ifelse(R_3, -2, -3)                                                   #Adder for the id difference between the symbol and any potential assignment operator for RIGHT_ASSIGN
      assignAdderLE <-  1                                                                   #Adder for the id difference between the symbol and any potential assignment operator for LEFT_ASSIGN or EQ_ASSIGN
      assignAdderR <- -2                                                                    #Adder for the id difference between the symbol and any potential assignment operator for RIGHT_ASSIGN
      for (i in 1:length(symbolIds)) {                                                      #For every symbol
        assignOpId <- symbolIds[i]+assignAdderLE                                            #  Get Id of potential assignment operator
#print(paste0("i=", i, " Sid=", symbolIds[i], " assignOpId_LE=", assignOpId))      
        if(assignOpId %in% tokens[,"id"]) {                                                 #  If that potential assignment operator exists 
          tokenOp <- tokens[tokens[, "id"]==assignOpId, col.token]                          #    Then get it's token
#print(paste0("tokenOp_LE=", tokenOp))
          bOutputVar[i] <- (tokenOp %in% c("LEFT_ASSIGN", "EQ_ASSIGN"))                     #    and Flag if this variable is being LEFT or EQ assigned
        }
        assignOpId <- symbolIds[i]+assignAdderR                                             #  Get Id of potential assignment operator
#print(paste0("i=", i, " Sid=", symbolIds[i], " assignOpId_R=", assignOpId))      
        if(assignOpId %in% tokens[,"id"]) {                                                 #  If that potential assignment operator exists 
          tokenOp <- tokens[tokens[, "id"]==assignOpId, col.token]                          #    Then get it's token
#print(paste0("tokenOp_R=", tokenOp))
          bOutputVar[i] <- bOutputVar[i] || (tokenOp=="RIGHT_ASSIGN")                       #    and Flag if this variable is being RIGHT assigned (being careful not to clear the flag when it was LEFT or EQ assigned)
        }
#print(paste0("bOutputVar[",i,"]=",bOutputVar[i]))      
      }
      dir <- append(dir, ifelse(bOutputVar, "Output", "Input"))                              #Add the direction for each new variable
    }
#    bOutputVar <- ((tokens[symbols+1, col.token]=="LEFT_ASSIGN") 
#                   | (tokens[symbols+1, col.token]=="EQ_ASSIGN")  
#                   | (tokens[abs(symbols-2), col.token]=="RIGHT_ASSIGN"))                   #Flag varaiables as outputs
    dtype <- append(dtype, rep("Numeric", length(symbolIds)))                               #Add the data type for each new variable
    ptype <- append(ptype, rep("Vector", length(symbolIds)))                                #Add the parameter type for each new variable
    param <- append(param, rep("", length(symbolIds)))                                      #Add the parameter for each new variable
    pVals <- append(pVals, rep("", length(symbolIds)))                                      #Add the parameter value for each new variable
    
    UpdateStatusBar()
    badNames <- grep("$", vars, fixed=TRUE)
    if(length(badNames)>0) {
      sNames <- paste(unlist(vars[badNames]), collapse=", ")
      gmessage(paste0("\"$\" detected in ", length(badNames), ifelse(length(badNames)>1," variables:\n", " variable:\n"), 
                      sNames, "\n\nVariables that use \"$\" are often not in the global environment and only variables in the global environment can be used as inputs and outputs.  Note that variables can be promoted to global by using the \"<<-\" operator or the \"assign\" function."),
               title="Non-Global Variable Warning", icon="warning")
    }
    
    #Now that we've collected all the variables (old and new), put them into a data frame
    envRAD$dfVar <- data.frame(vars, dir, dtype, ptype, param, pVals, 
                               stringsAsFactors=FALSE)                                      #Create variable data frame
    envRAD$dfVar <- envRAD$dfVar[!duplicated(envRAD$dfVar[ , 1]), ]                         #Delete any duplicate variables, keeping the first one
    rownames(envRAD$dfVar) <- envRAD$dfVar[, 1]                                             #Set the names of the rows in the variables table
    colnames(envRAD$dfVar) <- c("varName", "dir", "dType", "pType", "param", "pValue")      #Set the column headers for the variables table
    envRAD$deletedVars <- envRAD$dfVar[1, ]                                                 #Create a data frame for holding deleted variables
    envRAD$delVarMin <- nrow(envRAD$deletedVars)                                            #Remember row count at this time: this value means the deleted var list can be considered empty (no variables to un-delete)
    
    #Add the script and variables data frame to the GUI
    CreateDialog2()                                                                         #Create the dialog
    if(is.na(envRAD$workDr)) {                                                              #Did the script contain an existing Working Directory?
      svalue(envRAD$chkDir) <- FALSE                                                        #  No, so reset the Working Directory to it's default (FALSE)
      enabled(envRAD$edtDir) <- FALSE                                                       #      and Disable the Working Direcoty edit box
      svalue(envRAD$edtDir) <- gsub("/", "\\", dirname(envRAD$g.file), fixed=TRUE)          #      and set the directory to the default of the R Script's directory
    } else {
      svalue(envRAD$chkDir) <- TRUE                                                         #  Yes, so set the Working Directory to TRUE
      enabled(envRAD$edtDir) <- TRUE                                                        #      and Enable the Working Direcoty edit box
      svalue(envRAD$edtDir) <- envRAD$workDr                                                #      and set the directory 
    }
    if(is.na(envRAD$sortBy)) {                                                              #Did the script contain an existing SortBy value?
      svalue(envRAD$edtSort) <- envRAD$SORT_BY_DEFAULT                                      #  No, so use the default
    } else {
      svalue(envRAD$edtSort) <- envRAD$sortBy                                               #  Yes, so set it
    }
    svalue(envRAD$RptCt) <- rptCount                                                        #Set the repeat count
    svalue(envRAD$edtScr) <- envRAD$g.file                                                  #Display the script file location
    envRAD$allowUpdates <- TRUE                                                             #Allow updates now that the script has been processed
    SetOutputs()                                                                            #Set the outputs for the first time
    UpdateExpression(NA)                                                                    #Update the metric expression
    envRAD$saved <- TRUE                                                                    #Set the saved flag, we haven't modified any results (yet)
    envRAD$loaded <- TRUE                                                                   #Set the loaded flag, we have loaded an R Script 
    delete(envRAD$winRAD, envRAD$statusBar)                                                 #Remove the status bar, we don't need it after the dialog is loaded
  } #END-ScriptOpen
  
  ###-------------------------------------------------------------------------+
  ### Create main dialog
  CreateDialog2 <- function() {
    if(DEBUG<-FALSE) print(paste0("Function enter: CreateDialog2"))
    
    PARAM_DEFAULT         <- "{Default}"
    PARAM_DEFAULT_BOOLEAN <- "{FALSE}"
    PARAM_DEFAULT_NUMERIC <- "{Default=0}"
    PARAM_DEFAULT_STRING  <- "{Default=\"\"}"
    
    ###-------------------------------------------------------------------------+
    ### Persist the variables information after user actions
    UpdateVars <- function(varCtrl) {
      if(DEBUG<-FALSE) print(paste0("Function enter: UpdateVars, varName=", id(varCtrl)))
      envRAD$dfVar[(envRAD$dfVar[, 1]==id(varCtrl)), ] <- c(id(varCtrl), 
                               tag(varCtrl, "dir"), 
                               tag(varCtrl, "dType"), 
                               tag(varCtrl, "pType"), 
                               tag(varCtrl, "param"), 
                               tag(varCtrl, "pvalue"))
      UpdateExpression(NA)
      envRAD$saved <- FALSE 
    }

    ###-------------------------------------------------------------------------+
    ### Label the order of inputs and if they're repeated
    LabelInputs <- function(rptCt) {
      if(DEBUG<-FALSE) print("Function enter: LabelInputs")      
      inputCt <- length(envRAD$listIns)
      if(inputCt==0) return()
      nonRptCt <- inputCt-rptCt
      for(i in 1:inputCt) {
        z <- match(envRAD$listIns[i], rownames(envRAD$dfVar))
        svalue(envRAD$listButtons[[z]]) <- paste0(i, ") ", envRAD$listIns[i], 
                                                  ifelse((i>nonRptCt), "*", " "))
        font(envRAD$listButtons[[z]]) <- c(color="darkgreen", weight="bold")
      }
    }  #End-LabelInputs

    ###-------------------------------------------------------------------------+
    ### Label the first output as the Default
    LabelOutputs <- function() {
      if(DEBUG<-FALSE) print("Function enter: LabelOutputs")      
      outputCt <- length(envRAD$listOuts)
      for(i in 1:outputCt) {
        z <- match(envRAD$listOuts[i], rownames(envRAD$dfVar))
        if(i==1) {
          svalue(envRAD$listButtons[[z]]) <- paste0(" ", envRAD$listOuts[i], " (default)")
        } else {
          svalue(envRAD$listButtons[[z]]) <- paste0(" ", envRAD$listOuts[i], " ")
        }
        font(envRAD$listButtons[[z]]) <- c(color="tomato", weight="bold")
      }
    }  #END-LabelOutputs
      
    ###-------------------------------------------------------------------------+
    ### Set the font for the variable name
    ColorButton <- function(btn, dir) {
      if(DEBUG<-FALSE) print(paste0("Function enter: ColorButton, var=", svalue(btn), ", dir=", dir))
      spec <- switch(dir, 
                     Unused=c(color="black", weight="normal"), 
                     Input=c(color="darkgreen", weight="bold"),
                     Parameter=c(color="royalblue", weight="bold"), 
                     Output=c(color="tomato", weight="bold"))
      font(btn) <- spec
    }  #END-ColorButton
    
    ###-------------------------------------------------------------------------+
    ### Set the value for a parameter, making sure it's appropriate to the data type (Boolean, Numeric, String)
    SetParamValue <- function(edtPvalue, param) {
      if(DEBUG<-FALSE) print(paste0("Function enter: SetParamValue, param=", param, ", pvalue=", svalue(edtPvalue)))
      pvalue <- svalue(edtPvalue)
      if(pvalue %in% c(PARAM_DEFAULT, PARAM_DEFAULT_BOOLEAN, PARAM_DEFAULT_NUMERIC, PARAM_DEFAULT_STRING)) {
        pvalue=""
      }
      paramType <- substr(param, 1, 1)
      if(paramType=="B") {
        pText <- ifelse(tolower(substr(pvalue, 1, 4))=="true", "TRUE", PARAM_DEFAULT_BOOLEAN) 
      } else {
        if(paramType=="N") {
          pText <- ifelse(is.na(x <- as.numeric(pvalue)), PARAM_DEFAULT_NUMERIC, x)
        } else {
          pText <- gsub("[[:space:]]*", "", pvalue)
          if(length(pText)==0) {
            pText <- PARAM_DEFAULT_STRING
          } else {
            if(nchar(pText)==0) pText <- PARAM_DEFAULT_STRING
          }  
        }
      }
      svalue(edtPvalue) <- pText
      if(pText %in% c(PARAM_DEFAULT, PARAM_DEFAULT_BOOLEAN, PARAM_DEFAULT_NUMERIC, PARAM_DEFAULT_STRING)) {
        return("")
      } else {
        return(pText)
      }
    }  #END-SetParamValue
    
    ###-------------------------------------------------------------------------+
    ### Convert a parameter name to it's index in it's list, or return 1 if it's not there
    SetParamIndex <- function(param, pList) {
      if(DEBUG<-FALSE) print(paste0("Function enter: SetParamIndex, param=", param))
      return(ifelse(param %in% pList, grep(param, pList), 1))
    }  #END-SetParamIndex
    
    ###-------------------------------------------------------------------------+
    ### Set the list of avalable parameters
    SetParamCtrl <- function(dlParam, newParam, oldParam=NULL) {
      if(DEBUG<-FALSE) print(paste0("Function enter: SetParamCtrl, newParam=", newParam, ", oldParam=", oldParam))
      tempVal <- envRAD$allowUpdates
      envRAD$allowUpdates <- FALSE
      if(!is.null(oldParam)) {
        envRAD$paramUsed[grep(oldParam, envRAD$paramList)] <- FALSE         #Mark the old one as unused 
      }
      tempParamUsed <- envRAD$paramUsed
      if(!(is.null(newParam)||nchar(newParam)==0)) {
        if(newParam %in% envRAD$paramList) {
          envRAD$paramUsed[grep(newParam, envRAD$paramList)] <- TRUE        #Mark the new one as used 
          tempParamUsed[grep(newParam, envRAD$paramList)]<-FALSE            #Temporarily mark it as unused so it displays in the dropdown list
        }
      }
      dlParam[] <- envRAD$paramList[!tempParamUsed]                         #Update the parameter list
      svalue(dlParam) <- dlParam[SetParamIndex(newParam, dlParam[])]        #Set the selected parameter
      envRAD$allowUpdates <- tempVal
      return(svalue(dlParam))
    }  #END-SetParamCtrl
    
    ###-------------------------------------------------------------------------+
    ### Update the list of parameters for each variable's Parameter Dropdown List
    UpdateParamLists <- function() {
      if(DEBUG<-FALSE) print("Function enter: UpdateParamLists")      
      paramCt <- length(envRAD$listPars)
      if(paramCt==0) return()
      for(i in 1:paramCt) {
        z <- match(envRAD$listPars[i], rownames(envRAD$dfVar))
        dlParam <- envRAD$listParams[[z]]       
        SetParamCtrl(dlParam, tag(tag(dlParam, "varCtrl"), "param"))
      }
    }  #End-UpdateParamLists
    
    ###-------------------------------------------------------------------------+
    ### Set the maximum repeat count to match the number of inputs
    SetMaxRptCt <- function() {
      if(DEBUG<-FALSE) print("Function enter: SetMaxRptCt")
      n <- length(envRAD$listIns)                              #Get the number of inputs
      if(svalue(envRAD$RptCt)>n) {
        svalue(envRAD$RptCt) <- n     #Possibly lower the max Rpt Count next, make sure we can
      }
      envRAD$RptCt[] <- seq(0, n, by=1 )                       #Change the max Rpt Count
      LabelInputs(svalue(envRAD$RptCt))
    }  #END-SetMaxRptCt
    
    ###-------------------------------------------------------------------------+
    ### Update the appearance of a variable control when it's moved
    UpdateControl <- function(varCtrl) {
      if(DEBUG<-FALSE) print(paste0("Function enter: UpdateControl, varCtrl=", id(varCtrl)))
      envRAD$allowUpdates <- FALSE
      dir <- tag(varCtrl, "dir")
      ColorButton(tag(varCtrl, "btnVar"), dir)                  #Color the variable name appropriately
      if(dir!=envRAD$dirList[1]) {
        if(dir==envRAD$dirList[3]) {
          #tag(varCtrl, "param") <- SetParamCtrl(tag(varCtrl, "dlParam"), tag(varCtrl, "param"))
          #tag(varCtrl, "pvalue") <- SetParamValue(tag(tag(varCtrl, "dlParam"), "edtPvalue"), tag(varCtrl, "param"))
        } else {
          dType <- tag(varCtrl, "dType")
          if(dir==envRAD$dirList[4]) {
            dtList <- envRAD$dTypeList[-3]
          } else {
            dtList <- envRAD$dTypeList
          }
          i <- ifelse((dType %in% dtList), grep(dType, dtList), 1)
          dlDataType <- tag(varCtrl, "dlDataType")
          dlDataType[] <- dtList
          svalue(dlDataType, index=TRUE) <- i
        }  
      }
      envRAD$allowUpdates <- TRUE
    }  #END-UpdateControl
    
    ###-------------------------------------------------------------------------+
    ### Move a variable in the GUI
    MoveVar<- function(varCtrl, from, to) {
      if(DEBUG<-FALSE) print(paste0("Function enter: MoveVar, varCtrl=", id(varCtrl), " from=", from, " to=", to))
      visible(varCtrl, set=FALSE)                       #Hide this variable control
      envRAD$saved <- FALSE                             #Clear the saved flag
      outMoved <- FALSE                                 #Assume we're not moving an output, at least to start
      
      if(is.null(from)) from <- tag(varCtrl, "dir")     #From=NULL means get the "from" from the varCtrl
      
      #First, remove the variable from where it was from
      if(!is.na(from)) {                                #If there is a "from" (either passed in or converted from NULL)?
        if(from==envRAD$dirList[1]) {                   #  Is this variable currently unused
          i <- match(id(varCtrl), envRAD$listVars)      #    Find the index of this variable in list
          envRAD$listVars <- envRAD$listVars[-i]        #    Remove variable from list           
          delete(frame.vars, varCtrl)                   #    Remove variable from its frame
        } else {                                        #If there is a "from"
          if(from==envRAD$dirList[2]) {                 #  Is this variable currently an input
            i <- match(id(varCtrl), envRAD$listIns)     #    Find the index of this variable in list
            svalue(envRAD$listButtons[[match(id(varCtrl), 
                                             rownames(envRAD$dfVar))]]) <- paste0(" ", id(varCtrl), " ") #Clear any asterisk         
            envRAD$listIns <- envRAD$listIns[-i]        #    Remove variable from list           
            SetMaxRptCt()                               #    Set max repeat count
            delete(varCtrl, tag(varCtrl, "gVar"))       #    Remove Variable Controls
            delete(frame.ins, varCtrl)                  #    Remove variable from its frame
            } else {                                    #  Is this variable NOT currently an input
            if(from==envRAD$dirList[3]) {               #    Is this variable currently a parameter
              i <- match(id(varCtrl), envRAD$listPars)  #      Find the index of this variable in list
              envRAD$listPars <- envRAD$listPars[-i]    #      Remove variable from list
              delete(varCtrl, tag(varCtrl, "gPar"))     #      Remove Parameter Controls
              delete(frame.pars, varCtrl)               #      Remove variable from its frame
            } else {                                    #    Is this variable an output
              outMoved <- TRUE                          #      Flag that we're moving an ouput
              i <- match(id(varCtrl), envRAD$listOuts)  #      Find the index of this variable in list
              svalue(envRAD$listButtons[[match(id(varCtrl), 
                                               rownames(envRAD$dfVar))]]) <- paste0(" ", id(varCtrl), " ") #Clear any asterisk         
              envRAD$listOuts <- envRAD$listOuts[-i]    #      Remove variable from list
              delete(varCtrl, tag(varCtrl, "gVar"))     #      Remove Variable Controls
              delete(frame.outs, varCtrl)               #      Remove variable from its frame
            }
          }
        }
      }
      if(to==envRAD$dirList[1]) {                                  #If moving to Unused
        envRAD$listVars <- append(envRAD$listVars, id(varCtrl))    #  Add variable to end of Unused list
        add(frame.vars, varCtrl)                                   #  Add variable to Unused frame
      } else {                                                     #If not moving to Unused
        if(to==envRAD$dirList[2]) {                                #  If moving to Inputs
          envRAD$listIns <- append(envRAD$listIns, id(varCtrl))    #    Add variable to end of list
          SetMaxRptCt()                                            #    Update the repeat count as needed
          n <- length(envRAD$listIns)                              #    Get the number of inputs
          if(svalue(envRAD$RptCt)>n) svalue(envRAD$RptCt) <- n     #    Possibly lower the max Rpt Count next, make sure we can
          attr(envRAD$RptCt, "to") <- n                            #    Change the max Rpt Count
          add(frame.ins, varCtrl)                                  #    Add variable to it's frame
          add(varCtrl, tag(varCtrl, "gVar"))                       #    Configure control
        } else {                                                   #  If not moving to Inputs
          if(to==envRAD$dirList[3]) {                              #    If moving to Parameters
            envRAD$listPars <- append(envRAD$listPars, id(varCtrl))#      Add variable to end of list
            add(frame.pars, varCtrl)                               #      Add variable to it's frame
            add(varCtrl, tag(varCtrl, "gPar"))                     #      Configure control
            AddParamChangedHandler(tag(varCtrl, "dlParam"))        #      Add a handler for when the parameter is changed
            AddParamValueHandler(tag(varCtrl, "edtPvalue"))        #      Add a handler for when the parameter value is changed
          } else {                                                 #    If moving to Outputs
            envRAD$listOuts <- append(envRAD$listOuts, id(varCtrl))#      Add variable to end of list
            add(frame.outs, varCtrl)                               #      Add variable to it's frame
            add(varCtrl, tag(varCtrl, "gVar"))                     #      Configure control
            LabelOutputs()                                         #      Label the output
            outMoved <- TRUE                                       #      Flag that an output was moved
          }
        }
      }
      if(!is.na(from)) {                                           #Handle parameters that get moved
        if((from==envRAD$dirList[3])||(to==envRAD$dirList[3])) {   #  Did we move a parameter
          if(from!=to) {                                           #    Did that parameter change direction
            if(from==envRAD$dirList[3]) {                          #      Was it a parameter but now it's not
              envRAD$paramUsed[tag(varCtrl, "param")] <- FALSE     #        Mark this parameter as unused
              if(!(tag(varCtrl, "dType") %in% envRAD$dTypeList)) { #        If there's no valid data type
                tag(varCtrl, "dType") <- envRAD$dTypeList[1]       #          Set the data type to Numeric
              }
              if(!(tag(varCtrl, "pType") %in% envRAD$pTypeList)) { #        If there's no valid data type
                tag(varCtrl, "pType") <- envRAD$pTypeList[1]       #          Set the data type to Numeric
              }
            } else {                                               #      It's now parameter but before it wasn't
              param <- tag(varCtrl, "param")                       #        Get any previous parameter
              if (!(param %in% envRAD$paramList)) {                #          If this is not a valid parameter
                param <- names(envRAD$paramUsed[!envRAD$paramUsed][1])   #      Select the first unused parameter to use
                tag(varCtrl, "param") <- param                     #            Save this as the parameter
              } else {                                             #          This is a valid parameter
                if(envRAD$paramUsed[param]) {                      #            If this parameter is currently used
                  param <- names(envRAD$paramUsed[!envRAD$paramUsed][1]) #        Select the first unused parameter to use
                  tag(varCtrl, "param") <- param                   #              Save this as the parameter
                  tag(varCtrl, "pvalue") <- SetParamValue(tag(varCtrl, "edtPvalue"), param)  #              Set the value appropriately
                 }
              }
              envRAD$paramUsed[param] <- TRUE                      #        Mark this parameter as used                } else
            }
            UpdateParamLists()                                     #      Update the Parameter Lists that are in use           
          }
        }
      } else {
        if(to==envRAD$dirList[3]) UpdateParamLists()               #Update the Parameter Lists that are in use           
      }
      tag(varCtrl, "dir") <- to                                    #Persist new direction
      tag(tag(varCtrl, "btnVar"), "dir") <- to                     #Set the button's direction
      UpdateControl(varCtrl)                                       #Update the control so it matches it's new destination
      if(outMoved) SetOutputs()                                    #Update the Outputs
      if(!is.na(from)) UpdateVars(varCtrl)                         #Update variable, only if it was moved
      visible(varCtrl, set=TRUE) # show dialog                     #Make the control visible
    } #END-MoveVar
    
    ###-------------------------------------------------------------------------+
    ### Add a handler to the list of handlers to be dropped when dialog is destroyed
    AddDropHandler <- function(ctrl, obj, id) {
      dropHandlers <- tag(ctrl, DROP_HANDLERS)
      dropHandlers[[length(dropHandlers)+1]] <- list(obj=obj, id=id)
      tag(ctrl, DROP_HANDLERS) <- dropHandlers
    }  #END-AddDropHandler
    
    ParamValueHandler <- function(h, ...) {
      if(envRAD$allowUpdates) {         #Don't process actions while updates are NOT allowed
        if(DEBUG<-FALSE) print(paste0("Handler enter: edtPvalue, param=", svalue(h$obj)))
        if(envRAD$allowUpdates) {
          envRAD$allowUpdates <- FALSE
          tag(tag(h$obj, "varCtrl"), "pvalue") <- SetParamValue(h$obj, svalue(tag(h$obj, "dlParam")))
          envRAD$allowUpdates <- TRUE
          UpdateVars(tag(h$obj, "varCtrl"))
        }
      }
      return(FALSE)
    }  #END-ParamValueHandler

    AddParamValueHandler <- function(edtPvalue) {
      if(DEBUG<-FALSE) print(paste0("Function enter: AddParamValueHandler, varName=", id(tag(edtPvalue, "varCtrl"))))
      h.edtPvalue.Blur <- addHandlerBlur(edtPvalue, handler=ParamValueHandler)    #Add handler after selection is made
      AddDropHandler(tag(edtPvalue, "varCtrl"), edtPvalue, h.edtPvalue.Blur)
      tag(edtPvalue, "h.dlParam.Changed") <- h.edtPvalue.Blur
      if(DEBUG<-FALSE) print(paste0("Function exit: AddParamValueHandler, varName=", id(tag(edtPvalue, "varCtrl")), 
                                   " h.ID=", tag(edtPvalue, "h.edtPvalue.Blur")))
    }  #END-AddParamValuedHandler
    
    ###-------------------------------------------------------------------------+
    #Add handler so that if parameter changes, make sure the pValue is valid
    ParamChangedHandler <- function(h, ...) {  #Add handler after selection is made
      if(DEBUG<-FALSE) print(paste0("HANDLER_NEW: dlParam-Changed=", id(tag(h$obj,"varCtrl"))))
      if(envRAD$allowUpdates) {   #Don't process actions while updates are NOT allowed
        if(DEBUG<-FALSE) print(paste0("--> HANDLER_ENTER: dlParam-Changed=", id(tag(h$obj,"varCtrl"))))
        envRAD$allowUpdates <- FALSE
        tag(tag(h$obj, "varCtrl"), "param") <- SetParamCtrl(h$obj, svalue(h$obj), tag(tag(h$obj, "varCtrl"), "param"))
        tag(tag(h$obj, "varCtrl"), "pvalue") <- SetParamValue(tag(h$obj, "edtPvalue"), svalue(h$obj))
        envRAD$allowUpdates <- TRUE
        UpdateVars(tag(h$obj, "varCtrl"))
        UpdateParamLists()                                     #Update the Parameter Lists that are in use
        if(DEBUG<-FALSE) print(paste0("--> HANDLER_EXIT: dlParam-Changed=", id(tag(h$obj,"varCtrl"))))
      }
    }  ##END-ParamChangedHandler
    
    AddParamChangedHandler <- function(dlParam) {
      if(DEBUG<-FALSE) print(paste0("Function enter: AddParamChangedHandler, varName=", id(tag(dlParam, "varCtrl"))))
      h.dlParam.Changed <- addHandlerChanged(dlParam, handler=ParamChangedHandler)  #Add handler after selection is made
      AddDropHandler(tag(dlParam, "varCtrl"), dlParam, h.dlParam.Changed)
      tag(dlParam, "h.dlParam.Changed") <- h.dlParam.Changed
      if(DEBUG<-FALSE) print(paste0("Function exit: AddParamChangedHandler, varName=", id(tag(dlParam, "varCtrl")), 
                                   " h.ID=", tag(dlParam, "h.dlParam.Changed")))
    }  #END-AddParamChangedHandler
    
    ###-------------------------------------------------------------------------+
    ### Create a variable control
    CreateVarCtrl <- function(varName, dir, dType, pType, param, pvalue) {
      if(DEBUG<-FALSE) print(paste0("Function enter: CreateVarCtrl, varName=", varName))
      UpdateStatusBar()
      
      varCtrl <- ggroup(varName, horizontal=FALSE, spacing=0)
      id(varCtrl) <- varName
      tag(varCtrl, "dir") <- dir
      tag(varCtrl, "dType") <- dType
      tag(varCtrl, "pType") <- pType
      tag(varCtrl, "param") <- param
      tag(varCtrl, "pvalue") <- pvalue
      addHandlerUnrealize(varCtrl, handler = function(h,...) {      #Remove drop handlers when varCtrl is unrealized
        dropHandlers <- tag(varCtrl, DROP_HANDLERS)
        if(length(dropHandlers) > 0) {
          for(i in 1:length(dropHandlers)) {
            removehandler(dropHandlers[[i]]$obj, dropHandlers[[i]]$id)
          }
        }
      })
      
      #Add Upper Button
      btnVar <- gbutton(paste0(" ", varName, " "), container=varCtrl, expand=TRUE, icon="") 
      adddropsource(btnVar, targetType="object", 
                    handler=function(h, ...) getToolkitWidget(btnVar)$getParent(), varCtrl)
      tag(btnVar, "dir") <- dir
      tag(varCtrl, "btnVar") <- btnVar
      x <- match(varName, rownames(envRAD$dfVar))     
      envRAD$listButtons[x] <- btnVar
      
      #Create Lower Group for Variables
      gVar <- ggroup(container=NULL, spacing=0)
      tag(varCtrl, "gVar") <- gVar
      
      #Create Data Type dropdown
      dlDataType <- gdroplist(envRAD$dTypeList, container=gVar, selected=ifelse((dType %in% envRAD$dTypeList), match(dType, envRAD$dTypeList), 1), 
                              expand=TRUE, anchor=c(-1,0),
                              handler=function(h, ...) {
                                if(envRAD$allowUpdates) {    #Don't process actions while updates are NOT allowed
                                  tag(varCtrl, "dType") <- svalue(h$obj)
                                  UpdateVars(varCtrl)
                                }
                              })    
      tag(varCtrl, "dlDataType") <- dlDataType
      
      #Create Parameter Type dropdown
      pType <- tag(varCtrl, "pType")
      dlSize <- gdroplist(envRAD$pTypeList, selected=ifelse((pType %in% envRAD$pTypeList), match(pType, envRAD$pTypeList), 1), 
                          container=gVar, anchor=c(1,0), expand=TRUE,
                          handler=function(h, ...) {
                            if(envRAD$allowUpdates) {    #Don't process actions while updates are NOT allowed
                              tag(varCtrl, "pType") <- svalue(h$obj)
                              UpdateVars(varCtrl)    
                            }
                          })
      tag(varCtrl, "dlSize") <- dlSize
      
      #Create Lower Group for Parameters
      gPar <- ggroup(container=NULL, spacing=0)
      tag(varCtrl, "gPar") <- gPar
      
      dlParam <- gdroplist(envRAD$paramList[!envRAD$paramUsed], expand=TRUE, container=gPar)
      envRAD$listParams[x] <- dlParam
      tag(varCtrl, "dlParam") <- dlParam
      tag(dlParam, "varCtrl") <- varCtrl
      
      edtPvalue <- gedit(text=tag(varCtrl, "pvalue"), initial.msg=PARAM_DEFAULT, container=gPar, editable=TRUE, expand=TRUE, width=10, font.attr=c(size=8),
                         handler=function(h,...) {
                           if(envRAD$allowUpdates) {         #Don't process actions while updates are NOT allowed
                             if(DEBUG<-FALSE) print(paste0("edtPvalue Handler: ", svalue(dlParam), "=", svalue(h$obj)))
                             envRAD$allowUpdates <- FALSE
                             tag(tag(h$obj, "varCtrl"), "pvalue") <- SetParamValue(h$obj, svalue(dlParam))
                             envRAD$allowUpdates <- TRUE
                             UpdateVars(tag(h$obj, "varCtrl"))
                           }
                         }, action=dlParam)      
      tag(edtPvalue, "dlParam") <- dlParam
      tag(edtPvalue, "varCtrl") <- varCtrl
      tag(varCtrl, "edtPvalue") <- edtPvalue
      tag(dlParam, "edtPvalue") <- edtPvalue
      
      #Add handler so that if parameter changes, make sure the pValue is valid
      AddParamChangedHandler(dlParam)
      
      #Add handler so that if parameter value changes, it gets update when the text box loses focus
      AddParamValueHandler(edtPvalue) 
      
      return(varCtrl)
    } #END-CreateVarCtrl

    ###-------------------------------------------------------------------------+
    ### Create Primary User Interface 
    
    if (envRAD$loaded) {              #If the dialog was previously loaded,
      delete(envRAD$winRAD, envRAD$gRAD)     #  Then delete the gRAD group
      dispose(envRAD$gRAD)                   #  And dispose of it
    }  

    envRAD$allowUpdates <- FALSE
    
    envRAD$gRAD <- ggroup(horizontal=FALSE, container=NULL, use.scrollwindow=TRUE)     #Group container for main dialog
    add(envRAD$winRAD, envRAD$gRAD)
    visible(envRAD$gRAD) <- FALSE                                                           #Hide this group while it's being created
    addHandlerUnrealize(envRAD$gRAD, handler = function(h,...) {      #Remove drop handlers when gRAD is unrealized
      dropHandlers <- tag(envRAD$gRAD, DROP_HANDLERS)
      if(length(dropHandlers) > 0) {
        for(i in 1:length(dropHandlers)) {
          removehandler(dropHandlers[[i]]$obj, dropHandlers[[i]]$id)
        }
      }
    })

    fScript <- gframe("R Script File", container=envRAD$gRAD) #Frame to hold the Script Name
    envRAD$edtScr <- gedit("", editable=TRUE, container=fScript, expand=TRUE, fill="both")         #edtScr is the edit box with the Script name    
    id <- addHandlerKeystroke(envRAD$edtScr, handler=function(h, ...) {UpdateExpression(NA)})    #Add handler for editing the script file location
    AddDropHandler(envRAD$gRAD, envRAD$edtScr, id)
    envRAD$chkScr <- gcheckbox("Use R Scripts Folder", checked=FALSE, container=fScript,  
                               handler=function(h, ...) {                            #Handler for this checkbox
                                 envRAD$saved <- FALSE                               #  Remember this change
                                 UpdateExpression(NA)                                #  Update the expression
                               })                                                    #Checkbox, when checked, _RScriptFile is saved without a path

    top <- ggroup(container=envRAD$gRAD, expand=TRUE, fill="both")
    pg4 <- gpanedgroup(container=top, expand=TRUE)     #Group container for main dialog
    frame.vars <- gframe("Unused Variables", horizontal=FALSE, expand=FALSE, container=pg4, raise.on.dragmotion=TRUE)
    #gseparator(horizontal=FALSE, container=top, spacing=0, fill="y")
    pg3 <- gpanedgroup(container=pg4, expand=TRUE)     #Group container for main dialog
    frame.ins <- gframe("Inputs", horizontal=FALSE, expand=FALSE, container=pg3, fill="y", raise.on.dragmotion=TRUE)
    pg2 <- gpanedgroup(container=pg3, expand=TRUE)     #Group container for main dialog
    frame.pars <- gframe("Parameters", horizontal=FALSE, expand=FALSE, container=pg2, fill="y", raise.on.dragmotion=TRUE)
    pg1 <- gpanedgroup(container=pg2, expand=TRUE)     #Group container for main dialog
    frame.outs <- gframe("Outputs", horizontal=FALSE, expand=FALSE, container=pg1, fill="y", raise.on.dragmotion=TRUE)
    frame.metExp <- gframe("Metric Expression", horizontal=FALSE, expand=TRUE, container=pg1, fill="both")
    
    adddroptarget(frame.vars, targetType="object", handler=function(h,...) {
      #print(paste0("DROPPED on Vars"))
      MoveVar(h$dropdata, NULL, envRAD$dirList[1])
    })
    adddroptarget(frame.ins, targetType="object", handler=function(h,...) {
      #print(paste0("DROPPED on Inputs"))
      MoveVar(h$dropdata, NULL, envRAD$dirList[2])
    })
    adddroptarget(frame.pars, targetType="object", handler=function(h,...) {
      #print(paste0("DROPPED on Pars"))
      MoveVar(h$dropdata, NULL, envRAD$dirList[3])
    })
    adddroptarget(frame.outs, targetType="object", handler=function(h,...) {
      #print(paste0("DROPPED on Outs"))
      MoveVar(h$dropdata, NULL, envRAD$dirList[4])
    })

    ctl <- ggroup(container=frame.metExp, expand=TRUE, horizontal=FALSE) 
    gTop <- ggroup(container=ctl)
    gTopLeft <- ggroup(container=gTop, horizontal=FALSE)
    gseparator(container=gTop, horizontal=FALSE)
    gTopRight <- ggroup(container=gTop, horizontal=FALSE)
    envRAD$chkNul <- gcheckbox("Nulls Allowed", checked=TRUE, anchor=c(-1,0), container=gTopLeft, 
                               handler=function(h, ...) {
                                 envRAD$saved <- FALSE                                #Remember this change
                                 UpdateExpression(NA)})                               #Checkbox for Allowing Nulls
    gseparator(container=gTopLeft)
    envRAD$chkCIC <- gcheckbox("Check Input Count", checked=TRUE, anchor=c(-1,0), container=gTopLeft, 
                               handler=function(h, ...) {
                                 envRAD$saved <- FALSE                                #Remember this change
                                 UpdateExpression(NA)})                               #Checkbox for Check Input Count
    glabel("Repeat Count", container=gTopRight, editable=FALSE, anchor=c(-1,1))       #Two word label for Repeat Count
#    glabel("Repeat", container=gTopRight, editable=FALSE, anchor=c(-1,1))             #1st work label for Repeat Count
#    glabel("Count", container=gTopRight, editable=FALSE, anchor=c(-1,1))              #2nd word label for Repeat Count
    envRAD$RptCt <- gspinbutton(from=0, to=100, by=1, value=0, anchor=c(-1,-1), container=gTopRight, expand=TRUE, 
                                handler=function(h, ...) {
                                  LabelInputs(svalue(h$obj))
                                  svalue(lblLower) <- ifelse((svalue(h$obj)>0), "* = Repeated", "")
                                })                                                    #Spin button control for the Repeat Count
    lblLower <- glabel("", container=gTopRight, font=c(color="darkgreen", weight="bold"))
    gseparator(container=ctl)
    envRAD$chkSort <- gcheckbox("Enable Sort By", checked=TRUE, container=ctl,  
                                handler=function(h, ...) {                            #Handler for SortBy Checkbox
                                  enabled(envRAD$edtSort) <- svalue(envRAD$chkSort)   #  Enable SortBy edit box if checked, disable otherwise
                                  envRAD$saved <- FALSE                               #  Remember this change
                                  UpdateExpression(NA)                                #  Update the expression
                                })                                                    #Checkbox for SortBy
    envRAD$edtSort <- gedit(envRAD$SORT_BY_DEFAULT, editable=TRUE, expand=FALSE, container=ctl,
                            handler=function(h, ...) {
                              envRAD$saved <- FALSE                                    #  Remember this change
                              UpdateExpression(NA)})                                   #Edit box for SortBy value
    id <- addHandlerKeystroke(envRAD$edtSort, handler=function(h, ...) {UpdateExpression(NA)})    #Add handler for editing the working directory
    AddDropHandler(envRAD$gRAD, envRAD$edtSort, id)
    
    gseparator(container=ctl)
    envRAD$chkDir <- gcheckbox("Specify Working Directory", expand=FALSE, container=ctl, 
                               handler=function(h, ...) {                             #Handler for Working Directory Checkbox
                                 enabled(envRAD$edtDir) <- svalue(envRAD$chkDir)      #  Enable directory edit box if checked, disable otherwise
                                 envRAD$saved <- FALSE                                #  Remember this change
                                 UpdateExpression(NA)                                 #  Update the expression
                               })                                                     #Checkbox for setting the Working Directory and it's handler
    envRAD$edtDir <- gedit("", editable=TRUE, expand=FALSE, container=ctl,  
                           handler=function(h, ...) {UpdateExpression(NA)})           #Edit box for Working Directory
    addHandlerKeystroke(envRAD$edtDir, handler=function(h, ...) {UpdateExpression(NA)})    #Add handler for editing the working directory
    AddDropHandler(envRAD$gRAD, envRAD$edtDir, id)
    
    gseparator(container=ctl)
    glabel("Output Variable", container=ctl, anchor=c(-1,0))       #Label for Output Variable
    envRAD$lstOut <- gdroplist(envRAD$listOuts, editable=FALSE, container=ctl,
                        handler=function(h, ...) {UpdateExpression(NA)})       #Dropdown container for Outputs (list to be udpated as Variables are processed and changed)  
    fME <- gframe("Metric Expression", expand=TRUE, horizontal=FALSE, container=ctl)       #Container Frame for Metric Expression
    gbutton("Copy to Clipboard", container=fME, expand=FALSE, handler = function(h, ...) {          #Copy ME to Clipboard Handler
#      writeClipboard(svalue(envRAD$txtMetricExp))                                     #  If changes have been made but the script has not been saved, warn the user
      cat(svalue(envRAD$txtMetricExp), file = (con <- file( description="clipboard", open="w", encoding = "UTF-8")))
      close(con)
      if(!envRAD$saved) galert("Be sure to save your script!", title="Unsaved Changes", delay=3, widget=fME, icon="info")
    })                                                                              #Button to copy Metric Expression to Clipboard, and it's handler
    envRAD$txtMetricExp <- gtext("", container=fME, editable=FALSE, expand=TRUE, fill="both", do.autoscroll=TRUE)                      #Text box control for Metric Expression
    
    n <- nrow(envRAD$dfVar)
    envRAD$listButtons <- vector(length=n, mode="list")
    envRAD$listParams <- vector(length=n, mode="list")
    apply(envRAD$dfVar, 1, function(var) {
      varCtrl <- CreateVarCtrl(var[1], var[2], var[3], var[4], var[5], var[6])
      MoveVar(varCtrl, NA, var[2])
    })
    SetMaxRptCt()                          #Set the initial repeat count
    visible(top) <- TRUE
    visible(envRAD$gRAD) <- TRUE                  #UnHide after it's been created
    envRAD$allowUpdates <- TRUE
    return(window)
  }  #END-CreateDialog2  
  
  ###-------------------------------------------------------------------------+
  ### Check for required packages and, if not present, install them
  CheckInstallPackage <- function(packageName) {
    if(is.na(match(packageName, installed.packages()))) {   #is the package installed?
      install.packages(packageName)                         #install the package
    }  
  }  #END-CheckInstallPackage
  
  ###-------------------------------------------------------------------------+
  ### Main function
  
  #Configure required packages
#  if(!R_3){                     #If not R version 3.0.0 or later, then the parser package is required
#    CheckInstallPackage("parser")
#    require(parser)              #Require the Parser Package
#  }
  CheckInstallPackage("gWidgetsRGtk2")
  require(gWidgetsRGtk2)       #Require the gWidgetstRGtk2 Package
  options(guiToolkit = "RGtk2")
      
  ### Create R Analytic Deployer Dialog 
  envRAD <- new.env() 
  evalq({

	  #Global Variables
    TAG_METRIC_EXP        <- "  #Metric Expression: "
    SORT_BY_DEFAULT       <- "{Default=First Input}"

    #Member variables  
    dirListOLD <- c("Input", "Output", "Parameter", "Unused")
    dirList <- c("Unused", "Input", "Parameter", "Output")
    dTypeList <- c("Numeric", "String", "Default")
    pTypeList <- c("Vector", "Scalar")
    paramList <- c(paste("BooleanParam", seq(1, 9), sep = ""), 
                   paste("NumericParam", seq(1, 9), sep = ""), paste("StringParam", seq(1, 9), sep = ""))
    paramUsed <- rep(FALSE, length(paramList))
    names(paramUsed) <- paramList
    funcNameOLD <- c("RScriptSimple", "RScriptRelative", "RScriptAgg", 
                  "RScriptRelativeS", "RScriptAggS", "<--INVALID-->", 
                  "<--INVALID-->", "<--INVALID-->", "<--INVALID-->")    
    funcName <- c("RScriptSimple", "RScriptU", "RScriptAggU", 
                     "RScript", "RScriptAgg", "<--INVALID-->", 
                     "<--INVALID-->", "<--INVALID-->", "<--INVALID-->")    
    funcList <- c("Simple", "Relative (Unsorted)", "Aggregation (Unsorted)", 
                  "Relative", "Aggregation", "<--INVALID-->", 
                  "<--INCOMPLETE-->", "<--SORT BY ERROR-->", "<--OUTPUTS ERROR-->")
    funcDesc <- c("Scalar inputs & output (row at a time)", 
                  "Vector inputs & output (table at a time)", 
                  "Vector inputs & scalar output (grouping)",
                  "Vector inputs & output (table at a time)", 
                  "Vector inputs & scalar output (grouping)", 
                  "Scalar inputs & Vector output are not allowed", 
                  "Need at least one input & one output",
                  "Sort By requires that the first input is a vector",
                  "Outputs must be either all scalar or all vector")
    saved <- TRUE
    loaded <- FALSE
    
    if(DEBUG<-FALSE) print("RAD: Member variables Created")
      
    winRAD <- gwindow(paste0("deployR -- ver. ", VERSION), spacing=2, parent=c(0,0))       #Container for R Analytic Deployer Window
    size(winRAD) <- WINDOW_SIZE                                                            #Inital size 
    statusBar <- gstatusbar(">>> Loading")
    gtoolbar(list(                                                           #Define Toolbar
      open=gaction("Open", icon="open", , handler=function(...) {
        newOpen <- FALSE
        if(envRAD$saved) {
          newOpen <- TRUE
        } else {
          confirmDialog <- gconfirm("Changes have been made that have not been saved to your R Script.\n\nIf you really want to open a new script and lose these changes, click 'Ok' -- otherwise, click 'Cancel'.", 
                     title="Open without Saving?", icon="warning", 
                                    handler=function(h,...) {
                                      envRAD$saved <- TRUE
                                    })
        }
        if(envRAD$saved) ScriptOpen()
      }), 
      preview=gaction("Preview", icon="gtk-page-setup", handler=ScriptSave), 
      save=gaction("Save", icon="save", handler=SaveToFile), 
      help=gaction("Help", icon="help", handler=ShowHelp), 
      about=gaction("About", icon="about", handler=function(...) {
                      gmessage(title="deployR", 
                               paste("This utility captures the signature of your R analytic into a header block that's used by MicroStrategy when executing the R Script.", 
                                     "", "After capturing the signature and saving, simply copy the metric expression into the MicroStrategy metric editor to deploy.",
                                     "", COPYRIGHT, paste0("Version ", VERSION), sep = "\n"))
      }), 
      quit=gaction("Quit", icon="quit", handler=function(...) {
        if(!envRAD$saved) {
          confirmDialog <- gconfirm("Changes have been made that have not been saved to your R Script.\n\nIf you really want to quit and lose these changes, click 'Ok' -- otherwise, click 'Cancel'.", 
                                    title="Exit without Saving?", icon="warning", handler=dispose(winRAD))
        } else dispose(winRAD)
      })
    ), container=winRAD)
  }, envRAD)

  return(envRAD)  
  
}  #END-deployR

if(DEBUG<-FALSE) deployR()   #Auto-launch when in DEBUG mode
