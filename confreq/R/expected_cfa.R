#' @title Expected frequencies with glm 
#' @export expected_cfa
#' @description Calculates the expected frequencies of counts using log liniear model. 
#' @details No details
#' 
#' @param des a designmatrix (object of class \code{"matrix"}) as returned by function \code{design_cfg_cfa}.
#' @param observed a integer vector with \code{lenght(observed) == dim(des)[1]}. WARNING: The observed frequencies counts must be in an order corresponding to the coding sheme in designmatix (see argument \code{des}).
#' @param family argument passed to \code{\link{glm.fit}} with default set to \code{poisson()}
#' @param intercept argument passed to \code{glm.fit} with default set to \code{FALSE}
#' @param ... aditional arguments optional passed to \code{\link{glm.fit}} 
#' 
#' @return An vector object giving the expected counts.
#' @references No references in the moment
#' @examples #######################################
#' # expected counts for LienertLSD data example.
#' designmatrix<-design_cfg_cfa(kat=c(2,2,2)) # generate an designmatrix (only main effects)
#' data(LienertLSD) # load example data
#' observed<-LienertLSD[,4] # extract observed counts
#' expected_cfa(des=designmatrix, observed=observed) # calculation of expected counts
#' ####################################### 

############### start of function definition ##################
expected_cfa<-function(des,observed,family=poisson(), intercept=FALSE,...){
# func. by joerg-henrik heine jhheine(at)googlemail.com  
###############################################################
FIT<-glm.fit(x=des, y=observed ,family=family, intercept = intercept, ... ) # verglichen mit cfa.exe manual von Eye --> OK
exp.freq<-FIT$fitted.value

#cat("expected frequencies:", "\n")
return(exp.freq)
}
