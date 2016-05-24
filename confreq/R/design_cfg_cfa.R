#' @title Designmatrix for log linear CFA models 
#' @export design_cfg_cfa
#' @exportClass design_cfg_cfa
#' @description Calculates the designmatrix corresponding to a dataset with \code{length(kat)} columns (variables). 
#' @details This function internaly calls the function \code{pos_cfg_cfa}. 
#' 
#' For further information on designmatrices see decription on function \code{model.matrix} in the package \code{stats}. 
#'    
#' @param kat a numerical vector containing kardinal numbers, giving the number of categories for each variable of a dataset (in the respective order of the variables in such a dataset) which corresponds to the requested designmatrix. So the length of this numerical vector represents the number of variables.
#' @param form a character string which can be coerced into a model formulae with the function \code{as.formula} in the package \code{stats}. If this argument is left empty the function \code{design_cfg_cfa()} will return a designmatrix coding only main effects and no interactions -- for a designmatrix refering to  three variables for example, leaving the argument \code{form} empty will be equivalent to assigning the character \code{"~ V1 + V2 + V3"} to the argument (\code{form="~ V1 + V2 + V3"}).
#' 
#' A special Case is to define a null-model or rather a cfa model of order zero. In such a model no (main) effects are considered. This can be achieved bei passing the character expression \code{"null"} to the argument \code{form} -- so: \code{form = "null"} 
#' @param ... additional parameters passed through to function \code{model.matrix} in package \code{stats}.
#' @return A designmatrix - an object of class \code{c("matrix","design_cfg_cfa")} - for the formula therm given in argument\code{form}.
#' @references No references in the moment 
#' @examples #######################################
#' # designmatrix with three main effects.
#' # three variables with two categories each.
#' design_cfg_cfa(kat=c(2,2,2))
#' # two variables with two categories each and one variable
#' # with 7 categories (Linert LSD example).
#' design_cfg_cfa(kat=c(2,2,7))
#' ###########
#' # designmatrix with three main effects an three interactions.
#' # three variables with two categories each.
#' design_cfg_cfa(kat=c(2,2,2),form="~ V1 + V2 + V3 + V1:V2 + V1:V3 + V2:V3")
#' # two variables with two categories each and one variable
#' # with 7 categories (Linert LSD example).
#' design_cfg_cfa(kat=c(2,2,7),form="~ V1 + V2 + V3 + V1:V2 + V1:V3 + V2:V3")
#' #######################################

############### start of function definition ##################
design_cfg_cfa<-function(kat,form=paste("~", paste(paste("V",1:length(kat),sep=""),collapse=" + ")), ...){ 
# depends: function: pos_cfg_cfa  by joerg-henrik heine  
# func. by joerg-henrik heine jhheine(at)googlemail.com  
###############################################################
#  'desigmatrix' for cfa model of order zero
if (form== "null"){
  des<-matrix(rep(1,dim(pos_cfg_cfa(kat,fact=F))[1]),ncol=1)
}

#  'designmatrix' for other cfa models
if (form!= "null"){ 
old.o <- options(contrasts = c("contr.sum","contr.sum")) # set  and save some options first
d<-(pos_cfg_cfa(kat,fact=T))
form<-as.formula(form)
des<-model.matrix(form, d, ...)
options(old.o)# restore some options previously saved
}
class(des) <- c("matrix", "design_cfg_cfa")
return(des)
}
