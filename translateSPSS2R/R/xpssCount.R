#' Counts frequencies of specific observations
#'
#' R implementation of the SPSS \code{$COUNT} system variable.
#'
#' Count displays the frequencies of observations matching the count statement.
#'
#' @usage xpssCount(x, variables, count)
#'
#' @param x a (non-empty) data.frame or input data of class \code{xpssFrame}. 
#' @param variables atomic character or character vector with the name of the variables.
#' @param count atomic character or atomic numeric pattern.
#' @return A vector of the same length as x.
#' @author Bastian Wiessner
#' @examples
#' 
#' data(fromXPSS)
#' 
#' xpssCount(x=fromXPSS, 
#'    variables = "V1", count=list(exact="Nissan"))
#' 
#' xpssCount(x=fromXPSS, 
#'    variables = "V5", count=list(exact=2))
#'
#' xpssCount(fromXPSS, 
#'    variables = span(fromXPSS, from = "V5",
#'                              to = "V7_2"), 
#'                    count = list(from = "lo",
#'                                   to = 100))
#' @seealso Related Functions \code{\link{xpssAny}} \code{\link{\%in\%}} \code{\link{is.element}}
#' @export

xpssCount <- function(x,
                      variables,
                      count)
{
  
  for(i in 1:length(variables)){
    if(!(is.element(variables[i],names(x)))) {
      stop("The selected variables are not in the dataset")
    }
  }
  if(length(unlist(unique(lapply(x[,variables],class))))>1){
    stop("Occurence of strings or numeric values only")
  }
    
  functiontype <- "SB"
  x <- applyMetaCheck(x)
  
  x <- computeValue(x)
  
  #abfangen der special characters
  mintemp  <- maxtemp<- vector()
  if(!is.null(count)) {
    if("lo" %in% count) {
      for(i in 1:length(variables)){
        evalVAR <- eval(parse(text = paste("x$",variables[[i]],sep="")))
        mintemp <- c(min(evalVAR,na.rm=T),mintemp)
      }
      pos <- which(count %in% "lo" |count %in% "lowest")  
      count[[pos]] <- min(mintemp)
      
    }
    if("hi" %in% count) {
      for(i in 1:length(variables)){
        evalVAR <- eval(parse(text = paste("x$",variables[[i]],sep="")))
        maxtemp <- c(max(evalVAR,na.rm=T),mintemp)
      }
      pos <- which(count %in% "hi" |count %in% "highest")  
      count[[pos]] <- max(maxtemp,na.rm=T)
    }
  }
  x$count <- 0
  temp <- x
  Varvec <- vector()
  
  counter <- vector(length=length(x[,variables[1]]))
  #schleife mit den variablen einleiten
  if(length(count$exact)>0){
    for(l in 1:length(variables)){
      if(is.character(count)){
        if(count == "sysmis"|| count =="missing" || count == "nmiss"){
          counter[which(x[,variables[l]] %in% count)] <- counter[which(is.na(x[,variables[l]]))] +1  
        } 
        if(count == "nvalid"){
          counter[which(x[,variables[l]] %in% count)] <- counter[which(!(is.na(x[,variables[l]])))] +1  
        }        
      } else{
        counter[which(x[,variables[l]] %in% count)] <- counter[which(x[,variables[l]] %in% count)] +1  
      }
    }
    count <- counter
  }
  else {    
    count_min <- as.numeric(count$from)
    count_max <- as.numeric(count$to)
    for(l in 1:length(variables)){      
      evalVAR <- eval(parse(text = paste("x$",variables[[l]],sep="")))
      for(j in 1:length(evalVAR)){
        if((is.na(evalVAR[j]) == FALSE) && (round(count_min,digits=5) <= round(as.numeric(evalVAR[j]),digits=5)) && (round(as.numeric(evalVAR[j]),digits=5) <= round(count_max,digits=5))){
          x$count[[j]] <- x$count[[j]]+1
        }
      }
    }
    count <- x$count
  }
  return(count)
}
