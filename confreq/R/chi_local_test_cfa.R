#' @title Local Chi-Square Test 
#' @export chi_local_test_cfa
#' @description Calculates the local chi-square test based on obseved and expected frequencies. 
#' @details No details in the moment.
#' 
#' @param expected a vector giving the expected frequencies.
#' @param observed a vector giving the observed frequencies.
#' @return a list with chi-square statistic and corresponding degrees of freedom an p-value.
#' @references No references in the moment
#' @examples #######################################
#' # first calculate expected counts for LienertLSD data example.
#' designmatrix<-design_cfg_cfa(kat=c(2,2,2)) # generate an designmatrix (only main effects)
#' data(LienertLSD) # load example data
#' observed<-LienertLSD[,4] # extract observed counts
#' expected<-expected_cfa(des=designmatrix, observed=observed) # calculation of expected counts
#'  chi_local_test_cfa(observed,expected)
#' ####################################### 

############### start of function definition ##################
chi_local_test_cfa<-function(observed,expected){
# func. by joerg-henrik heine jhheine(at)googlemail.com  
###############################################################
loc.chi.square<-((observed-expected)^2)/expected
loc.df<-rep(1,length(loc.chi.square))
loc.chi.square.p<-1-pchisq(loc.chi.square,loc.df) 

erg<-list(Chi=loc.chi.square,df=loc.df,pChi=loc.chi.square.p)
####
#cat("local chi-square:")
return(erg)
}
