#' @title Possible configurations 
#' @export pos_cfg_cfa
#' @description Calculates all possible configuartions for some variables with different numbers of categories.
#' @details No details
#' 
#' @param kat a numerical vector containing kardinal numbers, giving the number of categories for each variable. 
#' So the length of this numerical vector represents the number of variables.
#' @param fact logical, default is \code{(fact=FALSE)}. If this argument is set to \code{(fact=TRUE)} the result is coerced to a data.frame with factor variables. 

#' @return An object of class "matrix" or "data.frame" (depending on the argument \code{fact}) containing all possible configurations for \code{lenght(kat)} variables with the respective number of categories given as kardinal numbers in the vector \code{kat}.
#' @references No references in the moment 
#' @examples #######################################
#' # possible configurations for ...
#' # three variables with two categories each (Linert LSD example).
#' pos_cfg_cfa(kat=c(2,2,2))
#' #######################################

############### start of function definition ##################
pos_cfg_cfa<-function(kat,fact=FALSE){
# func. by joerg-henrik heine jhheine(at)googlemail.com  
###############################################################
laenge=prod(kat)
  b <- matrix(0, laenge, length(kat))
if(length(names(kat))!=0){colnames(b)<-names(kat)}
 
for (i in 1:length(kat)){b[,i]<-rep(1:kat[i], each=(laenge/prod(kat[1:i])) ,length.out=laenge)}
if(fact==T){b<-as.data.frame(apply(b,2,factor))
           if(length(names(kat))==0){names(b)<-paste("V",1:length(kat),sep="")} }

return(b)
}
