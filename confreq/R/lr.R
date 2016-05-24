#' @title Likelihood Ratio Chi-square (LR) 
#' @export lr
#' @description Calculates the likelihod ratio chi-square statistic based on observed and expected counts. 
#' @details No details in the moment.
#' 
#' @param expected a vector giving the expected frequencies.
#' @param observed a vector giving the observed frequencies.
#' @return numeric giving the likelihood ratio chi-square statistic.
#' @references Stemmler, M. (2014). \emph{Person-Centered Methods -- Configural Frequency Analysis (CFA) and Other Methods for the Analysis of Contingency Tables}. Cham Heidelberg New York Dordrecht London: Springer.
#' 
#' @examples #######################################
#' ######### some examples ########
#' data(newborns)
#' newborns
#' designmatrix <- design_cfg_cfa(kat=c(2,2)) # generate an designmatrix (only main effects)
#' observed <- newborns[,3] # extract observed counts
#' expected <- expected_cfa(des=designmatrix, observed=observed) # calculation of expected counts
#' lr(observed,expected) # calculation of the likelihood ratio chi-square statistic
#' 

lr <- function(observed,expected){
  erg <- 2*(sum(observed*(log(observed/expected)),na.rm=T)) # check this for observed=0
  return(erg)
}