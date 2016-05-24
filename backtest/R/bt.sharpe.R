################################################################################
##
## $Id: bt.sharpe.R 1300 2008-08-27 21:01:11Z zhao $
##
## Returns a one-dimensional array of sharpe ratios
##
################################################################################

.bt.sharpe <- function(object){

 stopifnot(object@natural)

 ## creates the array to store the sharpe ratios

 sr <- array(dim = c(1, length(object@in.var)),
             dimnames = list("sharpe", object@in.var))

 ## handles single or multiple "in.var"

 for(i in object@in.var){
   x <- object@results[object@ret.var, i, , ,"means"]
   spreads <- x[,"high"] - x[,"low"]
   sr[1,i] <- mean(spreads, na.rm = TRUE)/sd(spreads, na.rm = TRUE)
 }
 sr
}
