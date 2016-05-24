library(fracdiff)
if(FALSE) # manual testing
library(fracdiff, lib="/u/maechler/R/Pkgs/fracdiff.Rcheck")

.ptime <- proc.time()
## Test if the default  'n.start' is ok, i.e., if the
## "burn in" period is long enough :

n <- 512

set.seed(101) ; ok <- TRUE
for(i in 1:2000) {
    r <- fracdiff.sim(n, ar = -0.9, ma = NULL, d = 0.3)$series
    if(max(abs(r)) > 10) {
        cat("OOps : indices ",str(ibig <- which(big <- abs(r) > 10)),
            "\n are > 10\n")
        if(any(ibig < 200) && (length(ibig) > 5 || abs(r)[big] > 20)) {
            cat("Some have index < 200 --> BREAK\n")
            ok <- FALSE
            break
        }
    }
    if(i %% 100 == 0) {
        cat(i,": ACF = \n")
        print(acf(r, plot=FALSE))
    }
}
if(!ok) {
    cat("i=",i," gave series \n")
    print(head(r)) ; cat(".......\n")
    plot(as.ts(r)) ## clearly did show problem {when we had bug}
}

## Try to find an example more quickly with setting `one seed':
.AR <- c(-.75, -.9)
.MA <- c(0.2, 0.1)
ok <- TRUE
set.seed(1)
r0 <- fracdiff.sim(100, d = 0.3)
r1 <- fracdiff.sim(100, ar = .AR, d = 0.25)
r2 <- fracdiff.sim(100, ar = .AR, ma = .MA, d = 0.2)
for(i in 1:1000) {
    set.seed(1)# yes; identical ones
    r0i <- fracdiff.sim(100, d = 0.3)
    r1i <- fracdiff.sim(100, ar = .AR, d = 0.25)
    r2i <- fracdiff.sim(100, ar = .AR, ma = .MA, d = 0.2)
    stopifnot(identical(r0, r0i),
              identical(r1, r1i),
              identical(r2, r2i))
}

## Last Line:
cat('Time elapsed: ', proc.time() - .ptime,'\n')
