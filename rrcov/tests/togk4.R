## VT::15.09.2013 - this will render the output independent
##  from the version of the package
suppressPackageStartupMessages(library(rrcov))

dodata <- function(nrep=1, time=FALSE, short=FALSE, full=TRUE, method = c("FASTMCD","MASS")){
    domcd <- function(x, xname, nrep=1){
        n <- dim(x)[1]
        p <- dim(x)[2]

        mcd<-CovOgk(x)

        xres <- sprintf("%3d %3d\n", dim(x)[1], dim(x)[2])

        lpad<-lname-nchar(xname)
        cat(pad.right(xname,lpad), xres)

        dist <- getDistance(mcd)
        quantiel <- qchisq(0.975, p)
        ibad <- which(dist >= quantiel)
        names(ibad) <- NULL
        nbad <- length(ibad)
        cat("Outliers: ",nbad,"\n")
        if(nbad > 0)
            print(ibad)
        cat("-------------\n")
        show(mcd)
        cat("--------------------------------------------------------\n")
    }

    lname <- 20

    ## VT::15.09.2013 - this will render the output independent
    ##  from the version of the package
    suppressPackageStartupMessages(library(rrcov))

    method <- match.arg(method)

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

    cat("Data Set               n   p  Half LOG(obj)        Time\n")
    cat("========================================================\n")
    domcd(heart[, 1:2], data(heart), nrep)
    domcd(starsCYG,data(starsCYG), nrep)
    domcd(data.matrix(subset(phosphor, select = -plant)), data(phosphor), nrep)
    domcd(stack.x,data(stackloss), nrep)
    domcd(data.matrix(subset(coleman, select = -Y)), data(coleman), nrep)
    domcd(data.matrix(subset(salinity, select = -Y)), data(salinity), nrep)
    domcd(data.matrix(subset(wood, select = -y)), data(wood), nrep)
    domcd(data.matrix(subset(hbk,  select = -Y)), data(hbk), nrep)

    domcd(brain, "Animals", nrep)
    domcd(milk, data(milk), nrep)
    domcd(bushfire, data(bushfire), nrep)
    cat("========================================================\n")
}

pad.right <- function(z, pads)
{
### Pads spaces to right of text
    padding <- paste(rep(" ", pads), collapse = "")
    paste(z, padding, sep = "")
}

dodata()
