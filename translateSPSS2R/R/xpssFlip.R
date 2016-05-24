#' Flips variables
#'
#' R Implementation of the SPSS \code{FLIP} Function. 
#' 
#' @param x a (non-empty) data.frame or input data of class "xpssFrame". 
#' @param variables atomic character or character vector with the names of the variables to flip
#' @param names atomic character with the name of the variable for coloumn names.
#' @return A flipped, respectively transposed xpssFrame object.
#' @author Bastian Wiessner
#' @seealso \code{\link{t}}
#' @export 
#' @examples 
#' data(fromXPSS)
#' xpssFlip(x=fromXPSS,variables=c("V4","V5","V6"),names="V1")

xpssFlip <- function(x, variables = "all", names = NULL){
  #variables <- c("V3","V4","V5","V6")
  #names <- "V6"
  functiontype <- "DM"
  x <- applyMetaCheck(x)
  options(warn = -1)
  if(variables == "all")    {
    ta <- t(x) 
    for(i in 1:length(x[[1]])){
   colnames(ta)[i] <- paste("var",i,sep="")
    }
  } else {
    for(i in 1:length(variables)) {
      if(!(is.element(variables[[i]],names(x)))) {
        stop("The selected variables has to be in the dataset")
      }  
    }
    ta <- t(x[variables])    
  }    
  ta <- as.xpssFrame(ta)
  
  #---------------------
  
  if(is.null(names)) {
    for(i in 1:length(x[[1]])){
      names(ta)[i] <- paste("var",i,sep="")
    }
  } else {
    evalVAR <- eval(parse(text=paste("x$",names,sep="")))
    for(i in 1:length(x[[1]])) {
      if(is.numeric(evalVAR)) {
        names(ta)[i] <- paste("K_", evalVAR[i],sep="")  
      } else {
        names(ta)[i] <- paste(evalVAR[i],sep="")  
      }      
    }
    hauef <- table(names(ta))
    for(i in 1:length(hauef)){
      id <- which(names(ta) %in% names(hauef[i]))
      for(j in 1:length(id)) {
        names(ta)[id][j] <- paste(names(ta)[id][j],"_",j,sep="")
      }
    }
    if(is.element(names,variables)){
      id <- which(!(is.element(variables,names)))
      ta <- ta[id,]      
    }
  }
  options(warn = 0)
  return(ta)  
}
