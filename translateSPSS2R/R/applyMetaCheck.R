#' Apply Meta Checks
#'
#' Checks the attributes for eventual merging the data
#'
#' Helper Function for Merging Attributedatasets
#' 
#' @param x a (non-empty) data.frame, data.table object or input data of class \code{"xpssFrame"}. 
#' @return Output is the input data merged on the basis of the condition in the attributes of the actual dataset.
#' @author Bastian Wiessner
#' @importFrom data.table is.data.table
#' @keywords internal
#' @export

applyMetaCheck <- function(x, pos = 1, envir = as.environment(pos)) {


  functiontype <- get("functiontype", parent.frame())
    
  if(functiontype == "ME") {
    if(!is.element("xpssFrame",class(x))) {
      message("actual data got coerced to xpssFrame, due missing of necessary attributes")
      x <- as.xpssFrame(x)
    } 
  }
  
  if(functiontype == "DM" || functiontype == "SB") {
    if(is.element("factor",lapply(x,class))){
      var <- names(which(lapply(x,class) == "factor"))
      for(i in 1:length(var)) {
        eval(parse(text =paste("x$",var[[i]]," <- as.character(x$",var[[i]],")",sep="")))  
      }
    }
  }
  
  
  if((length(attributes(x)$FILTER)>0) && (length(attributes(x)$TEMPORARY)>0) && (length(attributes(x)$SPLIT_FILE)>0) && (length(attributes(x)$DO_IF)>0) && (length(attributes(x)$SELECT_IF)>0) && (length(attributes(x)$WEIGHTS)>0)) {
    
    if(!is.element("xpssFrame",class(x))) {
      message("The actual data isn't a  xpssFrame object. The functionality of some functions is severely limited")
    }
    
    
    if(!is.null(attributes(x)$FILTER) && functiontype == "DM"){
      if(attributes(x)$FILTER != FALSE && "DM" %in% class(x)) {
        x <- rbind(x,attributes(x)$FILTERED_DATA)
        attBack <- attributesBackup(x)
        x <- x[order(as.numeric(rownames(x))),]      
        x <- applyAttributes(x, attBack)
        
        attributes(x)$FILTERED_DATA <- NULL
      }
    }
    
    if((attributes(x)$TEMPORARY == TRUE) && (functiontype == "AN")) { 
       
        attribut_backup <- attributesBackup(x)
        dataname <- get("dataname", parent.frame())
        globalVariables(c(dataname),package="translateSPSS2R",add=F)        
        if(getRversion() >= "3.1.0") utils::suppressForeignCheck(c(dataname))
        assign(x=dataname,value=attributes(x)$ORIGIN,envir= envir)
        eval(parse(text = paste0("attributes(", dataname,")$FILTER <- FALSE")), envir = .GlobalEnv)
        eval(parse(text = paste0("attributes(", dataname,")$TEMPORARY <- FALSE")), envir = .GlobalEnv)
        eval(parse(text = paste0("attributes(", dataname,")$SPLIT_FILE <- FALSE")), envir = .GlobalEnv)
        eval(parse(text = paste0("attributes(", dataname,")$DO_IF <- FALSE")), envir = .GlobalEnv)
        eval(parse(text = paste0("attributes(", dataname,")$SELECT_IF <- FALSE")), envir = .GlobalEnv)
        eval(parse(text = paste0("attributes(", dataname,")$WEIGHTS <- FALSE")), envir = .GlobalEnv)
        eval(parse(text = paste0("class(", dataname,") <- c('data.frame','xpssFrame')")), envir = .GlobalEnv)
        eval(parse(text = paste0("attributes(", dataname,")$ORIGIN <- NULL")), envir = .GlobalEnv)
    }
  } else {
    if(functiontype == "AN") {
      warning("actual data got coerced to xpssFrame, due missing of necessary attributes", call. =F)
      x <- as.xpssFrame(x)
    }
    warning("Essential attributes are missing in the actual object. This problem caused, whether you haven't created a xpssFrame object, due this issue the functionality of some functions is severely limited", call. = F)
  }
    return(x)
}
