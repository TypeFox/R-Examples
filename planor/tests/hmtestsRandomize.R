library("planor")
#---------------------------------------------------------------------------
# EXAMPLES TO CHECK THE RANDOMIZATION FUNCTION
#---------------------------------------------------------------------------
# Example 1: classical block design
#---------------------------------------------------------------------------
cat("\n")
cat("***************** RANDOMIZATION EXAMPLE 1 *****************\n")
cat("\n")
cat("A simple block design: 4 blocks and 7 treatments\n")
cat("Model: bloc+treatment\n")
cat("N=28\n")
#

randK <- planor.designkey(factors=c("bloc","treatment"),
                          nlevels=c(4,7),
                          model=~bloc+treatment, nunits=28,
                       base=~bloc, verbose=T)
randP <- planor.design(key=randK)@design
 planor.randomize(~bloc, randP) 
 planor.randomize(~bloc/UNITS, randP) 
#---------------------------------------------------------------------------
# Example 2: Latin square
#---------------------------------------------------------------------------
cat("\n")
cat("***************** RANDOMIZATION EXAMPLE 2 *****************\n")
cat("\n")
cat("A Latin square: 3 rows, 3 columns, 3 treatments\n")
cat("Model: row + column + treatment\n")
cat("N=9\n")
#
lsK <- planor.designkey(factors=c("row","col","treatment"),
                        nlevels=rep(3,3) ,
                        model=~row+col+treatment,
                        nunits=9,
                       base=~row+col, verbose=T)
lsP <- planor.design(key=lsK)@design
 planor.randomize(~row*col, lsP) 

