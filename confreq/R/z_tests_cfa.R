#' @title Two z-Approximation Tests
#' @export z_tests_cfa
#' @description Calculates the Chi-square approximation to the z-test and the binomial approximation to the z-test. 
#' @details An continuity correction can be applied to the binomial approximation -- see argument \code{ccor}.
#' @param expected a vector giving the expected frequencies.
#' @param observed a vector giving the observed frequencies.
#' @param ccor either a logical (TRUE / FALSE) determining wether to apply a continuity correction or not. When set to \code{ccor=TRUE} continuity correction is applied for expected values 5 =< expected =< 10. For \code{ccor=FALSE} no continuity correction is applied. Another option is to set \code{ccor=c(x,y)} where x is the lower and y the upper bound for expected values where continuity correction is applied. So \code{ccor=c(5,10)} is equivalent to \code{ccor=TRUE}. 
#' @param ntotal optional a numeric giving the total number of observations. By default ntotal is calculated as \code{ntotal=sum(observed)}. 
#' @return a list with z an p-values.
#' @references No references in the moment
#' @examples #######################################
#' # expected counts for LienertLSD data example.
#' designmatrix<-design_cfg_cfa(kat=c(2,2,2)) # generate an designmatrix (only main effects)
#' data(LienertLSD) # load example data
#' observed<-LienertLSD[,4] # extract observed counts
#' expected<-expected_cfa(des=designmatrix, observed=observed) # calculation of expected counts
#' z_tests_cfa(observed,expected)
#' ####################################### 

############### start of function definition ##################
z_tests_cfa<-function(observed,expected,ccor=FALSE,ntotal=sum(observed)){
# func. by joerg-henrik heine jhheine(at)googlemail.com  
###############################################################

# 3.4 Chi-Square Approximation to the z-Test
z.Chi<-(observed - expected) /sqrt(expected)
z.pChi<- 1 - pnorm(abs(z.Chi))

# 3.5 Binomial Approximation to the z-Test
if(class(ccor)=="numeric"){
  stopifnot(length(ccor)==2)
  continuity<-(expected<= ccor[2] & expected >= ccor[1])/2  
}
if(class(ccor)=="logical"){
  if(ccor==TRUE){ continuity <- (expected<= 10 & expected >= 5)/2 }
  if(ccor==FALSE){ continuity <- rep(0,length(expected))}
}
ccor.Binomial<-continuity!=0
z.Bin<-(observed - expected - continuity )/ sqrt(expected * (1 - expected/sum(observed)))
z.pBin <- 1 - pnorm(abs(z.Bin))

erg<-list(z.Chi=z.Chi,z.pChi=z.pChi,z.Bin=z.Bin,z.pBin=z.pBin,cor.=ccor.Binomial)

# cat("z and p-values:", "\n")
return(erg)
}
