#################
## multi-split ##
#################

library(hdi)

set.seed(123)

x <- matrix(rnorm(100*100), nrow = 100, ncol = 100)
y <- x[,1] + x[,2] + rnorm(100)

set.seed(3) ; fit.mult <- multi.split(x, y)
set.seed(3) ; fit.tmp <- multi.split(x, y, verbose = TRUE)

## dput(fit.mult$pval.corr)
stopifnot(all.equal(fit.mult$pval.corr,c(2.19211217621905e-10, 2.63511914584751e-08, 1, 1, 1, 1, 1, 
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
1, 1, 1, 1, 1, 1, 1, 1, 1)))

## dput(fit.mult$lci)
stopifnot(all.equal(fit.mult$lci, c(0.845556485400509, 0.592654394654161, -Inf, -Inf, -0.559330600307058, 
-Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, 
-Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -0.494058185775476, 
-Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -0.935987254184296, 
-0.686212365897482, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, 
-Inf, -0.477928536514776, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, 
-Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, 
-Inf, -Inf, -Inf, -0.740160526334972, -Inf, -Inf, -Inf, -0.558056531182565, 
-Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, 
-Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -0.476688695088987, 
-Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -0.353795495366226, 
-Inf)))

## dput(fit.mult$uci)
stopifnot(all.equal(fit.mult$uci, c(1.48387111041928, 1.26688746702801, Inf, Inf, 0.919205974134296, 
Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, 
Inf, Inf, Inf, Inf, Inf, Inf, Inf, 0.64068728225218, Inf, Inf, 
Inf, Inf, Inf, Inf, Inf, 0.782508357141595, 0.549789868734234, 
Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, 0.624958772229364, Inf, 
Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, 
Inf, Inf, Inf, Inf, Inf, Inf, 0.221415168229707, Inf, Inf, Inf, 
0.545315747796322, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, 
Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, 0.701205415738981, 
Inf, Inf, Inf, Inf, Inf, Inf, Inf, 0.637286715741768, Inf)))
 

                                         
