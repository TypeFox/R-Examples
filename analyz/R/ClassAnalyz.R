# -- CONSTRUCTOR
analyz <- function(){
  return(new("Analyz")) 
}
#' Class Analyz
#' 
#' @title Class to manage analysis steps described in a CSV file.
#'
#' @slot steps A data frame attribute with the steps loaded from a CSV file.
#' @slot nrRows A numeric attribute with the quantity of steps.
#' @slot nrColumns A numeric attribute with the quantity of items of the largest step definition.
#' @slot stepItems A list attribute with the items of a step.
#' @slot results A list attribute with the execution result of the steps.
#' 
#' @import methods
#' 
setClass( Class="Analyz",
          representation( steps      = "data.frame",
                          nrRows     = "numeric",
                          nrColumns  = "numeric",
                          stepItems  = "list",
                          results    = "list"),
          validity=function(object){ 
            #--- INSPECTOR
           if(TRUE){
             return(TRUE)
           }else{
             return(FALSE)
           }
          }
)
setMethod( f="initialize",
           signature="Analyz",
           definition=function(.Object, 
                               steps, 
                               nrRows, 
                               nrColumns,
                               stepItems,
                               results){ 
            #--- INITIALIZER
            # -- Set the attibutes with the defaults
            .Object@steps     <- data.frame()
            .Object@nrRows    <- 0
            .Object@nrColumns <- 0
            .Object@stepItems <- list()
            .Object@results   <- list()
            # -- Class inspection
            validObject(.Object)
            return(.Object) }
)
# STEP 1
#' Method Analyz.loadSteps
#' 
#' Method for reading the analysis file and fill the "steps" class spot.
#'
#' @param object Object instance.
#' @param path Path and analysis file name.
#' @return object Object instance.
#' 
#' @export
#' @docType methods
#' @rdname Analyz.loadSteps-methods
#' 
#' @examples
#' obj <- new("Analyz")
#' v_path <- vector()
#' Analyz.loadSteps(obj, v_path)
#' @export
#'
setGeneric("Analyz.loadSteps",
           function(object, path){standardGeneric("Analyz.loadSteps")})
#' @rdname Analyz.loadSteps-methods
#' @aliases Analyz.loadSteps,Analyz,Analyz-method
setMethod("Analyz.loadSteps",
          "Analyz",
          function(object, path)
  {
    # -- Constants   
    HEADER  <- FALSE
    FACTORS <- FALSE
    # -- BODY
    vSteps   <- NULL

vSteps  <- tryCatch(read.csv (
                              file      = path, 
                              header    = HEADER,
                              stringsAsFactors = FACTORS
                              ),
              error   = function(e) message(path, '\n', e),
              warning = function(w) message(path, '\n', w)
    )
    if(!is.null(vSteps)){
      object@nrRows    <- nrow(vSteps)
      object@nrColumns <- ncol(vSteps)
      
      object@steps<-vSteps                      
    }
    return(object)
  }
)
# STEP 2
#' Method Analyz.getNrRows 
#' 
#' Method for returning the number of rows read from the analysis file.
#' 
#' @param object Object instance.
#' @return nrRows Number of read rows from the analysis file.
#' 
#' @export
#' @docType methods
#' @rdname Analyz.getNrRows-methods
#' 
#' @examples
#' obj <- new("Analyz")
#' v_NrRow <- Analyz.getNrRows(obj)
#' @export
#' 
setGeneric("Analyz.getNrRows",
           function(object){ standardGeneric("Analyz.getNrRows") })
#' @rdname Analyz.getNrRows-methods
#' @aliases Analyz.getNrRows,Analyz,Analyz-method
setMethod("Analyz.getNrRows",
          "Analyz",
          function(object){ return(object@nrRows) }
)
# STEP 3
#' Method Analyz.getNrColumns
#' 
#' Method for returning the number of columns read from the analysis file.
#' 
#' @param object Object instance.
#' @return nrColumns Number of read columns from the analysis file.
#' 
#' @export
#' @docType methods
#' @rdname Analyz.getNrColumns-methods
#' 
#' @examples
#' obj <- new("Analyz")
#' v_nrCol <- Analyz.getNrColumns(obj)
#' @export
#' 
setGeneric("Analyz.getNrColumns",
           function(object){standardGeneric("Analyz.getNrColumns")})
#' @rdname Analyz.getNrColumns-methods
#' @aliases Analyz.getNrColumns,Analyz,Analyz-method
setMethod("Analyz.getNrColumns",
          "Analyz",
          function(object){ return(object@nrColumns) }
)
# STEP 4
#' Method Analyz.setStepItems
#' 
#' Method for interpret an analysis line and fill a class spot.
#' 
#' @param object Object instance.
#' @param index Index number of the step line to be read.
#' @return object Object instance.
#' 
#' @export
#' @docType methods
#' @rdname Analyz.setStepItems-methods
#' 
#' @examples
#' obj <- new("Analyz")
#' v_index <- numeric()
#' Analyz.setStepItems(obj, v_index)
#' @export
#' 
setGeneric("Analyz.setStepItems",
           function(object, index){standardGeneric("Analyz.setStepItems")})
#' @rdname Analyz.setStepItems-methods
#' @aliases Analyz.setStepItems,Analyz,Analyz-method
setMethod("Analyz.setStepItems",
          "Analyz", 
          function(object, index)
      {
      
      vTitle      <- c()
      vCommand    <- c()
      vStepSize   <- numeric()
      vParameters <- list()
      
      vStepLine   <- object@steps[index,]
      if(length(vStepLine) > 0 ){
        
        vTitle      <- vStepLine[[1]]
        vCommand    <- vStepLine[[2]]
        vStepSize   <- length(vStepLine)
        vParameters <- vStepLine[3:vStepSize]
        
        # Remove NAs
        vNAtab <- is.na(vParameters)
        vParameters <- vParameters[!vNAtab]  
        
      }      
      object@stepItems <- list("TITLE"=vTitle,
                               "COMMAND"=vCommand,
                               "PARAMETERS"=vParameters)
      
      return(object)
      }
)
# STEP 5
#' Method Analyz.getStepTitle
#' 
#' Method for returning the current step title.
#' 
#' @param object Object instance.
#' @return stepTitle Description of the current step title.
#' 
#' @export
#' @docType methods
#' @rdname Analyz.getStepTitle-methods
#' 
#' @examples
#' obj <- new("Analyz")
#' v_Title <- Analyz.getStepTitle(obj)
#' @export
#' 
setGeneric("Analyz.getStepTitle",
           function(object){standardGeneric("Analyz.getStepTitle")})
#' @rdname Analyz.getStepTitle-methods
#' @aliases Analyz.getStepTitle,Analyz,Analyz-method
setMethod("Analyz.getStepTitle",
          "Analyz",
          function(object) { return(unlist(object@stepItems["TITLE"])) })

# STEP 6
#' Method Analyz.getStepCommand
#' 
#' Method for returning the current step command.
#' 
#' @param object Object instance.
#' @return stepCommand  Description of the current step command.
#' 
#' @export
#' @docType methods
#' @rdname Analyz.getStepCommand-methods
#' 
#' @examples
#' obj <- new("Analyz")
#' v_Command <- Analyz.getStepCommand(obj)
#' @export
#' 
setGeneric("Analyz.getStepCommand",
           function(object){standardGeneric("Analyz.getStepCommand")})
#' @rdname Analyz.getStepCommand-methods
#' @aliases Analyz.getStepCommand,Analyz,Analyz-method
setMethod("Analyz.getStepCommand",
          "Analyz",
          function(object) { return(unlist(object@stepItems["COMMAND"])) })

# STEP 7
#' Method Analyz.getStepParameters
#' 
#' Method for returning the current step command parameters.
#' 
#' @param object Object instance.
#' @return stepParameters Description of the current step command parameters.
#' 
#' @export
#' @docType methods
#' @rdname Analyz.getStepParameters-methods
#' 
#' @examples
#' obj <- new("Analyz")
#' v_Parameters <- Analyz.getStepParameters(obj)
#' @export
#' 
setGeneric("Analyz.getStepParameters",
           function(object){standardGeneric("Analyz.getStepParameters")})
#' @rdname Analyz.getStepParameters-methods
#' @aliases Analyz.getStepParameters,Analyz,Analyz-method
setMethod("Analyz.getStepParameters",
          "Analyz",
          function(object) 
  { 
    vItems      <- unlist(object@stepItems["PARAMETERS"])
    vSize       <- length(vItems)
    vIsType     <- 1
    vType       <- c()
    vParameters <- list()
    vCount      <- 1
    vEnvironment <- ".GlobalEnv"
    
    if(vSize > 0){
      
      for(x in 1:vSize){
        if (vItems[x] != ""){
          
          if(vIsType ==1){ # Item is the parameter type
            vType   <- vItems[x]  
            vIsType <- 0
            
          }else{ # Item is the parameter itself
            switch(vType,
                   # At (get the value from the results)
                   "@"={
                     vParameters[vCount]<-list(Analyz.getResult( object, as.numeric(vItems[x]) ))
                   },
                   # Dolar (get value from the Globals)
                   "$"={
                     vParameters[vCount]<-list(get(vItems[x], envir=as.environment(vEnvironment)))
                   },
                    # Default (coerce the value)
                   {
                     vParameters[vCount]<-list(Analyz.coerceType( object, vItems[[x]], vType ))
                  }
              )
            vCount  <- vCount+1
            vIsType <- 1
          }
        }else{ NULL }
      }  
    }
    return(vParameters)
  })

#' Method Analyz.runAnalysis 
#' 
#' Method for step execution.
#' 
#' @param object Object instance.
#' @param command Description of the step command to be executed.
#' @param parameters Paramater(s) for the informed step command.
#' @return result Result of the step execution.
#' 
#' @export
#' @docType methods
#' @rdname Analyz.runAnalysis-methods
#' 
#' @examples
#' obj <- new("Analyz")
#' Analyz.runAnalysis(obj, 'mean', list(c(1,2,3)))
#' 
setGeneric("Analyz.runAnalysis",
           function(object, command, parameters){standardGeneric("Analyz.runAnalysis")})
#' @rdname Analyz.runAnalysis-methods
#' @aliases Analyz.runAnalysis,Analyz,Analyz-method
setMethod("Analyz.runAnalysis",
          "Analyz",
          function(object, command, parameters)
  { 
    vResult <- vInfo <- c()
    
    vInfo <- tryCatch(
        vResult <- do.call(command, parameters),
        error = function(e) return(e),
        warning = function(w) return(w) 
      )
    if(length(vResult) > 0){
      return( vResult )      
    }else{
      return(print(vInfo))
    }
  }
)
#' Method Analyz.setResult<-
#' 
#' Method for storing the informed result to the "results" class spot.
#' 
#' @param object Object instance.
#' @param value Result value of an execution.
#' @return object Object instance.
#' 
#' @export
#' @docType methods
#' @rdname Analyz.setResult-methods
#' 
#' @examples
#' obj <- new("Analyz")
#' v_result <- vector()
#' Analyz.setResult(obj) <- v_result
#' @export
setGeneric("Analyz.setResult<-",
           function(object, value){standardGeneric("Analyz.setResult<-")})
#' @rdname Analyz.setResult-methods
#' @aliases Analyz.setResult,Analyz,Analyz-method
setReplaceMethod( f="Analyz.setResult",
                  signature="Analyz",
                  definition=function(object, value){
                    # -- BODY
                    vIndex <- length(object@results)+1
                    object@results[vIndex] <- list(value)
                    return(object)
                  }
)
#' Method Analyz.getResult 
#' 
#' Method for getting the informed result from "results" class spot.
#' 
#' @param object Object instance.
#' @param index Index number of the results to be read.
#' @return result Result value of the informed index.
#' 
#' @export
#' @docType methods
#' @rdname Analyz.getResult-methods
#' 
#' @examples
#' obj <- new("Analyz")
#' vIdx <- numeric()
#' v_result <- Analyz.getResult(obj, vIdx)
#' @export
setGeneric("Analyz.getResult",
           function(object, index){standardGeneric("Analyz.getResult")})
#' @rdname Analyz.getResult-methods
#' @aliases Analyz.getResult,Analyz,Analyz-method
setMethod("Analyz.getResult",
          "Analyz",
          function(object, index)
  { 
    if( length(index) > 0 ){
      if(index <= object@nrRows){
        return( (object@results[[index]]) )
      }else{
        return( NULL )
      } 
    }else{
      return( NULL )
    }
  }
)
#' Method Analyz.coerceType
#' 
#' Method for variable type coercion.
#' 
#' @param object Object instance.
#' @param variable Variable to be coerced.
#' @param type Data type to coerce the informed variable.
#' @return result Coerced variable.
#' 
#' @export
#' @docType methods
#' @rdname Analyz.coerceType-methods
#' 
#' @examples
#' obj <- new("Analyz")
#' v_numeric <- Analyz.coerceType(obj, '2', 'numeric')
#' 
setGeneric("Analyz.coerceType",
           function(object, variable, type){standardGeneric("Analyz.coerceType")})
#' @rdname Analyz.coerceType-methods
#' @aliases Analyz.coerceType,Analyz,Analyz-method
#' 
setMethod("Analyz.coerceType",
          "Analyz",
          function(object, variable, type){
            result   <- FALSE
            # -- Coerce a variable to a defined type
            switch(type,
                   null={
                     tryCatch(result <- as.null(variable),
                              error = function(e) FALSE,
                              warning = function(w) FALSE)
                   },
                     numeric={
                       tryCatch(result <- as.numeric(variable),
                                error = function(e) FALSE,
                                warning = function(w) FALSE)
                     },
                     double={
                       tryCatch(result <- as.double(variable),
                                error = function(e) FALSE,
                                warning = function(w) FALSE)
                     },
                     character={
                       tryCatch(result <- as.character(variable),
                                error = function(e) FALSE,
                                warning = function(w) FALSE)
                     },
                     logical={
                       tryCatch(result <- as.logical(variable),
                           error = function(e) FALSE,
                           warning = function(w) FALSE)
                     },
                     vector={
                       tryCatch(result <- c(variable),
                                error = function(e) FALSE,
                                warning = function(w) FALSE)
                     },
                     list={
                       tryCatch(result <- list(variable),
                                error = function(e) FALSE,
                                warning = function(w) FALSE)
                     },
                     result <- variable
            )
            return( result )
          }
)