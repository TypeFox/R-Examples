# 2010-05-25 CJS fixed error in se of Petersen estimator

SimplePetersen <- function( n1, m2, u2) {
#
#   Estimate abundance of unmarked fish at second sample occasion using the SimplePetersen estimator
#   n1 - number of animals marked
#   m2 - number of recaptures of those in n1
#   u2 - number of animals unmarked
#   All three elements can be vectors in which case a vector of estimates (and se) are returned.
#
#   NOTE that the estimate is for UNMARKED at the second sample occasion. Many people want
#   the total abundance at the second occasion in which case you need to add the n1 to the total.
# 
#   The Chapman estimator is returned with its se
# 
#   Output: a list with elements
#        $est and $se
#

   est <- (n1+1)*(u2+1)/(m2+1) - 1
   se  <- sqrt((n1+1)*(m2+u2+1)*(n1-m2)*(u2)/(m2+1)^2/(m2+2))
   list(est=est, se=se)
} #end of function


TestIfPool <- function(n1, m2){
#
#   Test if can use a pooled Petersen experiment by seeing if the recapture rate is
#   equal across all strata. This is a simple chi-square test to see if the proportion
#   recovered are equal.
#
#   Input:
#       n1  - vector of marked fish released
#       m2  - vector of marks recovered
# 
#   Output
#       $chi - an object of type "htest" with components $chi$statistic with the X2 statistic
#              and $chi$o,value containing the p.value. [More components are returned.]
#       $fisher - an object of type "htest" with components $fisher$p.value with the p.value (not currently done)
#
#
#   browser()
    options(warn=-1)  # turn off warnng display about small Chi-equare values  
    temp <- cbind(n1-m2, m2)
    chi <- chisq.test(temp)
    fisher <- NULL   # fisher.test(temp, simulate.p.value=TRUE)
    list(chi=chi, fisher=fisher)
} # end of function

