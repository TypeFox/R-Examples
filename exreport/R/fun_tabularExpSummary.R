#' Summarize the experiment with a table for given outputs
#'
#' This function generates a table for the given outputs of the experiment, 
#' comparing all methods for each one of the problems. In addition the function
#' can highlight the best results for each problem as well as display a range of
#' parameters for the posterior renderization.
#'
#' @export
#' @param exp The ource experiment to generate the table from
#' @param outputs A given variable or list of them to be the target of the table
#' @param boldfaceColumns Indicate ("none","max" or "min") to highlight the
#' method optimizing the variables for each problem.
#' @param format Indicates the format of the numeric output using C formatting 
#' styles. Defaults to 'f'
#' @param digits The number of decimal digits to include for the numeric output.
#' @param tableSplit Indicates the number of parititions of the table that
#' will be rendered. Usefull when the the table is excessivelly wide.
#' @param rowsAsMethod Display the methods as the rows of the table, indicate 
#' FALSE for a transpose table.
#' @return An extabular object
#' 
#' @examples
#' 
#' # This example plots the distribution of the trainingTime variable in the 
#' # wekaExperiment problem.
#' 
#' # First we create the experiment from the problem.
#' experiment <- expCreate(wekaExperiment, name="test", parameter="fold")
#' 
#' # Next we must process it to have an unique parameter configuration:
#' # We select a value for the parameter featureSelection:
#' experiment <- expSubset(experiment, list(featureSelection = "yes"))
#' # Then we reduce the fold parameter:
#' experiment <- expReduce(experiment, "fold", mean)
#' # Finally we remove unary parameters by instantiation:
#' experiment <- expInstantiate(experiment, removeUnary=TRUE)
#' 
#' # Generate the default table:
#' tabularExpSummary(experiment, "accuracy")
#' 
tabularExpSummary <- function(exp, outputs, boldfaceColumns="none", format="f", digits=4, tableSplit=1, rowsAsMethod=TRUE) {
  
  if( !is.experiment(exp) )
    stop(.callErrorMessage("wrongParameterError", "exp", "experiment"))
  
  if( !is.character(outputs) )
    stop(.callErrorMessage("wrongParameterError", "outputs", "character"))
  
  if (!all(outputs %in% exp$outputs))
    stop(.callErrorMessage("variableNotPresentError", outputs))
  
  if( !is.character(boldfaceColumns) )
    stop(.callErrorMessage("wrongParameterError", "boldfaceColumns", "character"))
  
  if( !all(boldfaceColumns %in% c("none", "max", "min")) )
    stop(.callErrorMessage("parameterInRangeError","boldfaceColumns", "\"max\",\"min\" or \"none\""))
  
  if( !is.character(format) )
    stop(.callErrorMessage("wrongParameterError", "format", "character"))
  
  if( !all(format %in% c("a","A","d","i","f","e","E","g","G","o","s","x","X")) )
    stop(.callErrorMessage("parameterInRangeError","format", "[aAdifeEgGosxX]"))
  
  if( !is.numeric(digits) )
    stop(.callErrorMessage("wrongParameterError","digits", "numeric"))
  
  if( !is.numeric(tableSplit) )
    stop(.callErrorMessage("wrongParameterError", "tableSplit", "numeric"))
  
  if( tableSplit<1 || tableSplit==Inf)
    stop(.callErrorMessage("parameterInRangeError","tableSplit", "[1...Inf)"))
  
  if( !is.logical(rowsAsMethod) )
    stop(.callErrorMessage("wrongParameterError", "rowsAsMethod", "logical"))
  
  
  # Check if instantiation of parameters is needed
  if( length(exp$parameters) != 0 )
    stop(.callErrorMessage("requireInstantiationError"))
  
  tables <- list()
  formats <- list()
  
  for(i in 1:length(outputs)){
    table <- reshape2::dcast(exp$data, paste(exp$method, "~", exp$problem), value.var = outputs[i])
    # We set the index value for boldfaceColumns, format and digits as they may be single numbers or arrays.
    # In case of index out of range, we use the last value for the rest of elements.
    ibf <- i
    if(i>length(boldfaceColumns))
      ibf <- length(boldfaceColumns)
    ifmt <- i
    if(i>length(format))
      ifmt <- length(format)
    idig <- i
    if(i>length(digits))
      idig <- length(digits)
      
    frank <- function(...) { max(...,na.rm = TRUE) }
    mode  <- "eq"
    if (boldfaceColumns[ibf] == "min")
      frank <- function(...) { min(...,na.rm = TRUE) }
    else if (boldfaceColumns[ibf] == "none")
      mode <- "none"
    
    tableFormat <- apply(table[,-1,drop=FALSE],2,FUN = function(v){
      if(mode!="none"){
        #In that case we specify the number of digits
        if(format[i] %in% c("f","e","E","g","G","a","A"))
          ifelse(v==frank(v),paste0("b%.",digits[idig],format[ifmt]),paste0("%.",digits[idig],format[ifmt]))
        else
          ifelse(v==frank(v),paste0("b%",format[ifmt]),paste0("%.",digits[idig],format[ifmt]))
      }
      else{
        #In that case we specify the number of digits
        if(format[ifmt] %in% c("f","e","E","g","G","a","A"))
          rep(paste0("%.",digits[idig],format[ifmt]),length(v))   
        else
          rep(paste0("%",format[ifmt]),length(v))   
      }
    })
    tableFormat <- data.frame(matrix(tableFormat,nrow = nrow(table)))
    #We introduce the first column for method names.
    tableFormat <- cbind(rep("%s",nrow(tableFormat)),tableFormat)
    colnames(tableFormat) <- colnames(table)
    
    tables[[outputs[i]]] <- table
    formats[[outputs[i]]] <- tableFormat
  }
  
  title <- sprintf("Detailed results for output(s) %s", paste0("\"",outputs,"\"", collapse = ", "))  
  
  tab <- .exTabular(tables = tables, 
                    target = toString(outputs),
                    formats=formats, 
                    tableSplit=tableSplit, 
                    tableType="plain", 
                    title = title,
                    alias = "ExpSummaryTable",
                    tags = exp$tags)
  
  if(!rowsAsMethod)
    tab <- .mlTableTranspose(tab, exp$problem)
  
  tab
}