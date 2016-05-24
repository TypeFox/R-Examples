
#### Testing  medcouple  and related functions

### Strict (and timing) tests are in ./mc-strict.R
###                                   ~~~~~~~~~~~~
### Here, we produce output which is *compared* with ./mc-etc.Rout.save

library(robustbase)
source(system.file("xtraR/mcnaive.R", package = "robustbase"))# mcNaive()

## This is somewhat interesting {diff the output !}
## particularly since *most* give  the  'not found' diagnostic
set.seed(17)
for(n in 1:100) {
    cat(sprintf("n =%3d:\n------\n", n))
    mcval <- mc(rlnorm(n), trace=TRUE, doRefl=FALSE)
    cat(sprintf(" --> mc(rlnorm(%d)) = %.6f\n", n, mcval))
}

allEQ <- function(x,y) all.equal(x,y, tolerance = 1e-12)

x3 <- c(-2, rep(-1,4), rep(0,6), 2, 2, 2:4)
mcNaive(x3,"h.use") # 1/3
mcNaive(x3,"simple")#  0

try( mc(x3, doRefl = FALSE, maxit = 15, trace = 3)) ## "non-convergence" (32-bit)
str(robustbase:::mcComp(-x3, doRefl = FALSE, maxit = 15, trace = 4))

### And here is the "real" problem of the whole 'eps' idea:

x4 <- c(1:5,7,10,15,25, 1e15) ## this is also in mc-strict.R (but differently)
mcNaive(x4,"h.use") # 0.5833333
mcNaive(x4,"simple")# == " == 7/12
try( mc(x4) )# not converged  !!
str(robustbase:::mcComp( x4, doRefl= FALSE, maxit = 15, trace= 3))## = 0: conv.quickly
str(robustbase:::mcComp(-x4, doRefl= FALSE, maxit = 15, trace= 3)) # *not* conv!

## ## a much more extreme eps seems the cure:
## str(robustbase:::mcComp( x4, doRefl= FALSE, eps1=.Machine$double.xmin))
## str(robustbase:::mcComp(-x4, doRefl= FALSE, eps1=.Machine$double.xmin))

### Examples "like x3" (non-convergence on 32-bit)
xClist <- list(## length 5 :
               c(0,0, 1, 3,3),
               c(0,0, 1, 3:4),
               ##
               ## length 6 :
               c(0,0, 2, 4:6),    c(0,0, 2, 3, 4, 6), c(0,0, 4, 5, 7, 8),
               c(0, 1,1, 2, 6,6), c(0, 3,3, 4, 6,6),
               c(0,0, 1, 3, 5,5),
               c(0,0, 1, 4,4, 6), c(0,0, 1, 4,4, 7),  c(0,0, 1, 5,5, 6),

               ## n = 9 :
               c(-2,-2,-2, -1,-1, 1,1,1, 3),
               c(-3,-1,-1,  0, 1, 2,2,2,2)
               )

rlis <- lapply(xClist, function(x)
               try(mc(x, maxit=9), silent=TRUE))
table(sapply(rlis, class))
## if(R.version$arch == "x86_64") {
    print(unlist(rlis))
    rl2 <- lapply(xClist, mc, maxit=9) ##, eps1= 1e-10)
    stopifnot(allEQ(rlis, rl2),
              allEQ(unlist(rlis), sapply(xClist, mcNaive)))
##}


set.seed(47)
for(n in 3:60) {
    cat(" ")
    x <- round(2 * rnorm(n)) # many ties, often at median -- not converging
    ## if(R.version$arch == "x86_64") {
        ## non-convergence BUG  rarely and only on 32-bit (solved, MK)
        mc1 <- mc(x)
        mc2 <- mcNaive(x, method = "simple")
        mc3 <- mcNaive(x, method = "h.use")
        stopifnot(allEQ(mc1, mc3))
        if(mc2 != mc3) {
            cat("d"); if(!isTRUE(allEQ(mc2, mc3))) cat("!!")
        }
    ## }
    cat(".")
};  cat("\n")


cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''

quit('no')
##  ------

## Find short example of non-convergence (32-bit) --> for above xClist
n <- 9
for(ii in 1:100) {
    x <- round(2 * rnorm(n)) # many ties, often at median -- not converging
    mc1 <- mc(x)
}
##
x5 <- c(-3, -3, -2, -1, -1, 0, 0, 1, 2, 2, 3, 4)
x6 <- c(-5, -2, -1, -1, -1, 0, 0, 0, 2, 2, 2, 4)

