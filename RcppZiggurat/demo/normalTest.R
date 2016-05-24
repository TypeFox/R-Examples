
## This modifies the approach in (R)DieHarder which does
##   take N draws from a U(0,1)
##   repeat M times
##   and for large enough N, then the sum of all N draws goes to
##       mean   --> N/2
##       stddev --> sqrt(N/12)
##   which is known as the Irwin-Hall distribution
##   then for each of these M values use the inverse of normal to obtain a p-value
##   that p value should be uniformly distributed across these M draws
##   so use Kuiper's K/S test variant to test for uniform U(0,1)
##
## Here we don't need Irwin-Hall: the sum of N vars drawn as N(0,1) will be N(0,sqrt(N))
## So we compute a p value from that and assemple M such p values

library(RcppZiggurat)

norres <- RcppZiggurat:::normalTest(N=1e5,      		# individual draws
                                    M=1e2,  		    # repeats pre draw
                                    seed=123456789,
                                    generators=c("Ziggurat", "MT", "LZLLV", "GSL", "V1", "QL"),
                                    showplot=interactive())



