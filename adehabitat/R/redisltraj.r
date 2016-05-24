"redisltraj" <- function(l, u, burst = NULL, samplex0 = FALSE,
                         addbit = FALSE, nnew=5)
{
    ## Verifications
    if (!inherits(l, "ltraj"))
        stop("l should be of class 'ltraj'")
    if (is.null(burst)) {
        burst <- unlist(lapply(l, function(x) attr(x, "burst")))
    }
    ml <- l

    ## The function to be applied to each burst
    foo <- function(bu) {

        ## remove the missing values
        l <- ml[burst=bu]
        x <- na.omit(l[[1]]$x)
        y <- na.omit(l[[1]]$y)
        dat <- as.numeric(l[[1]]$dat[!is.na(l[[1]]$x)])

        ## Should the first relocation be sampled along the first step?
        if (samplex0) {

            ## Computation of the slope and intercept for the first step
            pente <- (y[2] - y[1]) / (x[2] - x[1])
            ori <- y[2] - pente*x[2]

            ## sample x and y at t=0
            if (x[1] <= x[2])
                x1 <- runif(1,x[1],x[2])
            if (x[1] > x[2])
                x1 <- runif(1,x[2],x[1])
            y1 <- pente*x1 + ori
            x0 <- c(x1,y1)

            ## linear interpolation for the date
            di1 = sqrt(  (x0[1] - x[1])^2 + (x0[2] - y[1])^2 )
            R = sqrt((x[2] - x[1])^2 + (y[2] - y[1])^2)
            di2 = dat[2] - dat[1]
            dat0 = dat[1] + (di1 * di2 / R)
        } else {

            ## if not sampled, take the first one
            x0 <- c(x[1],y[1])
            dat0 <- dat[1]
        }
        n <- length(x)
        nn <- nnew*n

        ## External call to the C function "discretrajr"
        toto <- .C("discretrajr", as.double(x), as.double(y), as.double(dat),
                   double(nn), double(nn), as.integer(n),
                   as.integer(nn), double(nn), as.double(x0[1]),
                   as.double(x0[2]), as.double(u), as.double(dat0), integer(1),
                   PACKAGE = "adehabitat")

        ## Number of steps
        neff <- toto[[13]] - 1
        if (neff >= (nn-1))
            stop("too small rediscretization step length. Try to increase \"nnew\"")

        ## The coordinates and the dates
        x <- toto[[4]][1:neff]
        y <- toto[[5]][1:neff]
        dat <- toto[[8]][1:neff]
        class(dat) <- c("POSIXt", "POSIXct")

        ## Should the final fragment of step be added
        if (addbit) {
            x <- c(x, l[[1]]$x[length(l[[1]]$x)])
            y <- c(y, l[[1]]$y[length(l[[1]]$y)])
            dat <- c(dat, l[[1]]$date[length(l[[1]]$date)])
        }

        ## Converts to ltraj
        opt <- options(warn=-1)
        if (!attr(ml, "typeII")) {
            nl <- as.ltraj(data.frame(x,y), id=attr(l[[1]], "id"),
                           burst = paste(attr(l[[1]], "burst"),
                           ".R",u,sep=""), typeII=FALSE)
        } else {
            nl <- as.ltraj(data.frame(x,y), id=attr(l[[1]], "id"),
                           date = dat, burst = paste(attr(l[[1]], "burst"),
                           ".R",u,sep=""))
        }
        nl[[1]]$rel.ang[is.na(nl[[1]]$rel.ang)] <- 0

        ## Output
        class(nl) <- "ltraj"
        options(opt)
        return(nl)
    }

    ## applies the function to all bursts and pool the results
    nl <- do.call("c.ltraj",lapply(burst,foo))
    return(nl)
}

