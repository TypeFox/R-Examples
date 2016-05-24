## VT::15.09.2013 - this will render the output independent
##  from the version of the package
suppressPackageStartupMessages(library(rrcov))

library(MASS)
dodata <- function(nrep = 1, time = FALSE, full = TRUE) {
    domest <- function(x, xname, nrep = 1) {
        n <- dim(x)[1]
        p <- dim(x)[2]
        mm <- CovMest(x)
        crit <- log(mm@crit)
##        c1 <- mm@psi@c1
##        M <- mm$psi@M

        xres <- sprintf("%3d %3d %12.6f\n", dim(x)[1], dim(x)[2], crit)
        lpad <- lname-nchar(xname)
        cat(pad.right(xname,lpad), xres)

        dist <- getDistance(mm)
        quantiel <- qchisq(0.975, p)
        ibad <- which(dist >= quantiel)
        names(ibad) <- NULL
        nbad <- length(ibad)
        cat("Outliers: ",nbad,"\n")
        if(nbad > 0)
            print(ibad)
        cat("-------------\n")
        show(mm)
        cat("--------------------------------------------------------\n")
    }

    options(digits = 5)
    set.seed(101) # <<-- sub-sampling algorithm now based on R's RNG and seed

    lname <- 20

    data(heart)
    data(starsCYG)
    data(phosphor)
    data(stackloss)
    data(coleman)
    data(salinity)
    data(wood)
    data(hbk)

    data(Animals, package = "MASS")
    brain <- Animals[c(1:24, 26:25, 27:28),]
    data(milk)
    data(bushfire)

    tmp <- sys.call()
    cat("\nCall: ", deparse(substitute(tmp)),"\n")

    cat("Data Set               n   p       c1        M     LOG(det)       Time\n")
    cat("======================================================================\n")
    domest(heart[, 1:2], data(heart), nrep)
    domest(starsCYG, data(starsCYG), nrep)
    domest(data.matrix(subset(phosphor, select = -plant)), data(phosphor), nrep)
    domest(stack.x, data(stackloss), nrep)
    domest(data.matrix(subset(coleman, select = -Y)), data(coleman), nrep)
    domest(data.matrix(subset(salinity, select = -Y)), data(salinity), nrep)
    domest(data.matrix(subset(wood, select = -y)), data(wood), nrep)
    domest(data.matrix(subset(hbk,  select = -Y)), data(hbk), nrep)


    domest(brain, "Animals", nrep)
    domest(milk, data(milk), nrep)
    domest(bushfire, data(bushfire), nrep)
    cat("======================================================================\n")
}

#   generate contaminated data using the function gendata with different
#   number of outliers and check if the M-estimate breaks - i.e. the
#   largest eigenvalue is larger than e.g. 5.
#   For n=50 and p=10 and d=5 the M-estimate can break for number of
#   outliers grater than 20.
dogen <- function(){
    eig <- vector("numeric",26)
    for(i in 0:25) {
        gg <- gendata(eps=i)
        mm <- CovMest(gg$x, t0=gg$tgood, S0=gg$sgood, arp=0.001)
        eig[i+1] <- ev <- getEvals(mm)[1]
        # cat(i, ev, "\n")

        stopifnot(ev < 5 || i > 20)
    }
    # plot(0:25, eig, type="l", xlab="Number of outliers", ylab="Largest Eigenvalue")
}

#
# generate data 50x10 as multivariate normal N(0,I) and add
# eps % outliers by adding d=5.0 to each component.
#   - if eps <0 and eps <=0.5, the number of outliers is eps*n
#   - if eps >= 1, it is the number of outliers
# - use the center and cov of the good data as good start
# - use the center and the cov of all data as a bad start
#   If using a good  start, the M-estimate must iterate to
#   the good solution: the largest eigenvalue is less then e.g. 5
#
gendata <- function(n=50, p=10, eps=0, d=5.0){

    if(eps < 0 || eps > 0.5 && eps < 1.0 || eps > 0.5*n)
        stop("eps is out of range")

    library(MASS)

    x <- mvrnorm(n, rep(0,p), diag(p))
    bad <- vector("numeric")
    nbad = if(eps < 1) eps*n else eps
    if(nbad > 0){
        bad <- sample(n, nbad)
        x[bad,] <- x[bad,] + d
    }
    cov1 <- cov.wt(x)
    cov2 <- if(nbad <= 0) cov1 else cov.wt(x[-bad,])

    list(x=x, bad=sort(bad), tgood=cov2$center, sgood=cov2$cov, tbad=cov1$center, sbad=cov1$cov)
}

pad.right <- function(z, pads)
{
    ## Pads spaces to right of text
    padding <- paste(rep(" ", pads), collapse = "")
    paste(z, padding, sep = "")
}


## -- now do it:
dodata()
dogen()
#cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
