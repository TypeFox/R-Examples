#' @title Degrees of freedom 
#' @export df_des_cfa
#' @description Calculates the degrees of freedom based on an designmatrix for a (log liniear) CFA model. 
#' #' @details No details
#' 
#' @param des a designmatrix (object of class "matrix") as returned by function \code{design_cfg_cfa}.
#' 
#' @return An object of class "integer" giving the degrees of freedom for the designmatrix defined in argument \code{des}.
#' @references No references in the moment
#' @examples #######################################
#' # degrees of freedom for designmatrix with three main effects.
#' # three variables with two categories each.
#' df_des_cfa(design_cfg_cfa(kat=c(2,2,2)))
#' # two variables with two categories each and one variable
#' # with 7 categories (Linert LSD example).
#' df_des_cfa(design_cfg_cfa(kat=c(2,2,7)))
#' ###########
#' # degrees of freedom for designmatrix with three main effects
#' # and three 'two by two'interactions.
#' # and tripple interaction --> saturated model --> df=0
#' # three variables with two categories each.
#' df_des_cfa(design_cfg_cfa(kat=c(2,2,2),form="~ V1 + V2 + V3 + V1:V2 + V1:V3 + V2:V3 + V1:V2:V3"))
#' ####################################### 

############### start of function definition ##################
df_des_cfa<-function(des){
# func. by joerg-henrik heine jhheine(at)googlemail.com  
###############################################################
cfg_row<-dim(des)[1] #how many different configurations in'des'
cfg_col<-dim(des)[2] #how many maineffects (variables) + interactions in'des'  
df<-  (cfg_row) - (cfg_col)
#names(df)<-"df"
return(df)
}
