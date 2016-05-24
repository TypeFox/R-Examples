library(expm)

## Missing REPROTECT(), till 2014-09-03 [because 'A' is *integer*]:
set.seed(17)
n <- 300
A <- matrix(rbinom(n^2, size=1, prob=0.1), n,n)
A2 <- A %^% 2
for(i in 1:100) {
    A. <- A %^% 2
    if(!isTRUE(all.equal(A2, A.)))
        cat("not equal; i=",i,"\n")
}
## MM: On nb-mm3, I get a different error which shows memory corruption:
##     REAL() can only be applied to a 'numeric', not a 'character'
## or  REAL() can only be applied to a 'numeric', not a 'NULL'
