#' Defines missing values for variables.
#'
#' R implementation of the SPSS \code{MISSING VALUES} function. xpssMissingValues defines values as missing and replaces them with \code{NA}. Position and Value are stored in the attributes of the specific variables.
#' 
#' xpssMissingValues specifies values for missing data for the selected variables. Those variables which match the terms of beeing a missing data get treated as \code{NA}. In most cases, variables which contain \code{NA} receive a special treatment in data management, case selection, and descriptive, respectively inductive statistics.  
#' User-missing values and system-missing values get treated as exactly one kind of missing data. The only difference in those missing values are that system missings get automatically assigned by the program when no legal value can be produced (e.g. character input at a numeric varibale, failed datatransformation) and user-defined missings, which are missing user data (e.g. the respondent forgot to answer, or skipped the question). \cr \cr Common is that this empty spaces are filled with \emph{-9 till -999} (for e.g. refusal to respond, inability to respond, Non-contact).
#'
#' The \code{as.missing} statement indicates the handling of values which are matched by the as.missing statement. Input format is a list with the arguments \code{range} to determine a range of values with the arguments \code{from} and \code{to} as NA or \code{singlevalues} to specify one more singlevalues as missing.
#'
#' \strong{NOTE:} The special arguments \code{lo} and  \code{hi} can be used to determine the lowest and highest value of a numeric value, wheter a missing \code{range} gets indexed.
#'
#' @usage xpssMissingValues(x, variables = NULL, 
#' as.missing = list(range = c(from=NULL,to=NULL),singlevalues = NULL), 
#' append = FALSE)
#' 
#' @param x a (non-empty) data.frame, data.table object or input data of class \code{"xpssFrame"}. 
#' @param variables atomic character or character vector with the names of the variables.
#' @param as.missing numeric list containing range and singlevalues.
#' @param range numeric vector containing a missing range from i to n.
#' @param singlevalues atomic numeric or numeric vectors containing singlevalues which determine missing values.
#' @param append logical. Indicating, if the existing missings should get overwritten or not.
#' @return a \code{xpssFrame} object with \code{NAs} located at the position where the specified values in as.missing used to be. In the attributes of the object the position and the value itself is stored. 
#' @author Andreas Wygrabek
#' @examples 
#' data(fromXPSS)
#' 
#' fromXPSS <- xpssMissingValues(fromXPSS, 
#' variables = "V6", 
#' as.missing = list(range=c(from="lo",
#' to=45)))
#' 
#' fromXPSS <- xpssMissingValues(fromXPSS, 
#' variables = "V3", 
#' as.missing = list(singlevalues=c(1,
#' 2)))
#' 
#' fromXPSS <- xpssMissingValues(fromXPSS, 
#' variables = "V6", 
#' as.missing = list(singlevalues="lo",
#' range=c(from="50",
#' to="hi")))
#' 
#' @export
xpssMissingValues <- function(x, variables = NULL, as.missing = list(range = c(from=NULL,to=NULL), singlevalues = NULL), append = FALSE){
  
  missing_range <- list()
  ####################################################################
  ####################################################################
  
  functiontype <- "DM"
  x <- applyMetaCheck(x)
  
  ####################################################################
  ####################################################################
  ####################################################################
  
  
  
  
  for(i in 1:length(variables)){
    x[,variables[i]] <- computeValue(x,variables[i]) 
  } 
      LEN <- length(variables)
        for(i in 1:length(variables)){
        VAR <- paste("x$",variables[i], sep = "")
        evalVAR <- eval(parse(text = paste(VAR)))
        
        
        #abfangen der special characters
        if(!is.null(as.missing$range)) {
          if("lo" %in% as.missing$range) {
            pos <- which(as.missing$range %in%  "lo" | as.missing$range %in% "lowest")  
            as.missing$range[[pos]] <- min(evalVAR,na.rm=T)
          }
          if("hi" %in% as.missing$range) {
            pos <- which(as.missing$range %in% "hi" |as.missing$range %in% "highest")  
            as.missing$range[[pos]] <- max(evalVAR,na.rm=T)
          }
        }
        
        #abfangen der special characters
        if(!is.null(as.missing$singlevalues)) {
          if("lo" %in% as.missing$singlevalues) {
            pos <- which(as.missing$singlevalues %in%  "lo" | as.missing$singlevalues %in% "lowest")  
            as.missing$singlevalues[[pos]] <- min(evalVAR,na.rm=T)
          }
          if("hi" %in% as.missing$singlevalues) {
            pos <- which(as.missing$singlevalues %in% "hi" |as.missing$singlevalues %in% "highest")  
            as.missing$singlevalues[[pos]] <- max(evalVAR,na.rm=T)
          }
        }
        
        
        if(!is.null(as.missing$range)) {
          missing_min <- as.numeric(as.missing$range[[1]])
          missing_max <- as.numeric(as.missing$range[[2]])
          for(j in 1:length(evalVAR)){
            if((is.na(evalVAR[j]) == FALSE) && (round(missing_min,digits=5) <= round(evalVAR[j],digits=5)) && (round(evalVAR[j],digits=5) <= round(missing_max,digits=5)))
            {
              missing_range[[j]] <- evalVAR[j]
            }
          }
        }
        missings <- c(as.missing$singlevalues,unlist(missing_range))
              

        if(eval(parse(text = paste("is.null(attributes(x$",variables[i],")$defined.MIS) | (!is.null(attributes(x$",variables[i],")$defined.MIS) & !append)", sep = "")))){
          
                
          if(length(as.missing$singlevalues)> 0 && length(as.missing$range)> 0) {
            eval(parse(text=paste("attr(x$",variables[i],",'defined.MIS') <- list(values = as.missing$singlevalues,range= as.missing$range)", sep ="")))
            
          } else if((length(as.missing$singlevalues)> 0)  && (is.null(as.missing$range))) {
           
            eval(parse(text=paste("attr(x$",variables[i],",'defined.MIS') <- list(values= as.missing$singlevalues)", sep ="")))
           
          } else {
            eval(parse(text=paste("attr(x$",variables[i],",'defined.MIS') <- list(range= as.missing$range)", sep ="")))
          }
          
        } else if(eval(parse(text = paste("!is.null(attributes(x$",variables[i],")$defined.MIS) & append", sep = "")))){ 
          
          
          if(length(as.missing$singlevalues)> 0 && length(as.missing$range)> 0) {
            eval(parse(text=paste("attr(x$",variables[i],",'defined.MIS') <- c(attributes(x$",variables[i],")$defined.MIS,list(values=unique(as.missing$singlevalues), range= as.missing$range))", sep ="")))
            
          } else if((length(as.missing$singlevalues)> 0)  && (is.null(as.missing$range))) {
                    
                    eval(parse(text=paste("attr(x$",variables[i],",'defined.MIS') <- c(attributes(x$",variables[i],")$defined.MIS,list(values=unique(as.missing$singlevalues)))", sep ="")))
                    
          } else {
            eval(parse(text=paste("attr(x$",variables[i],",'defined.MIS') <- c(attributes(x$",variables[i],")$defined.MIS,list(range=unique(as.missing$range)))", sep ="")))
          }
          
        } 
        
        
        
        if(is.null(eval(parse(text = paste("attributes(x$",variables[i],")$MIS", sep = ""))))){ 

            POS <- which(is.element(evalVAR,missings))
            VAL <- evalVAR[which(is.element(evalVAR,missings))]
            newMIS <- cbind(POS,VAL)
            
            eval(parse(text = paste("attr(x","$",variables[i],",'MIS') <- newMIS", sep = "")))
  
         
            } else if(append){
                
                eval(parse(text = paste("OBJ <- attributes(x$",variables[i],")$MIS", sep = "")))
                
                POS <- which(is.element(evalVAR,missings))
                VAL <- evalVAR[which(is.element(evalVAR,missings))]
                newMIS <- cbind(POS,VAL)
                MIS <- eval(parse(text=paste("rbind(attributes(x$",variables[i],")$MIS, newMIS)", sep = "")))
                MIS <- unique(MIS)
                MIS <- MIS[order(MIS[,1]),]
           
                eval(parse(text = paste("attr(x$",variables[i],",'MIS') <- MIS", sep = "")))
                
            } else {
                
                eval(parse(text = paste("x$",variables[i],"[attributes(x$",variables[i],")$MIS[,1]] <- attributes(x$",variables[i],")$MIS[,2]", sep ="")))
                evalVAR <- eval(parse(text = paste(VAR)))
                
  
                POS <- which(is.element(evalVAR,missings))
                VAL <- evalVAR[which(is.element(evalVAR,missings))]
                newMIS <- cbind(POS,VAL)
                  
                eval(parse(text = paste("attr(x","$",variables[i],",'MIS') <- newMIS", sep = "")))
            }
}
  
   for(i in 1:length(variables)){
       
       eval(parse(text = paste("x$",variables[i],"[attributes(x$",variables[i],")$MIS[,1]] <- NA", sep = "")))
       
   }


x <- applyAttributeDemerge(x)
  
  return(x)
}


