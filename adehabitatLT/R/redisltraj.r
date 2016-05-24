redisltraj <- function (l, u, burst = NULL, samplex0 = FALSE, addbit = FALSE,
                        nnew = 5, type = c("space","time"))
{
    if (!inherits(l, "ltraj"))
        stop("l should be of class 'ltraj'")
    type <- match.arg(type)

    if (type=="space") {
        if (is.null(burst)) {
            burst <- unlist(lapply(l, function(x) attr(x, "burst")))
        }
        ml <- l
        foo <- function(bu) {
            l <- ml[burst = bu]
            x <- na.omit(l[[1]]$x)
            y <- na.omit(l[[1]]$y)
            dat <- as.numeric(l[[1]]$dat[!is.na(l[[1]]$x)])
            if (samplex0) {
                pente <- (y[2] - y[1])/(x[2] - x[1])
                ori <- y[2] - pente * x[2]
                if (x[1] <= x[2])
                    x1 <- runif(1, x[1], x[2])
                if (x[1] > x[2])
                    x1 <- runif(1, x[2], x[1])
                y1 <- pente * x1 + ori
                x0 <- c(x1, y1)
                di1 = sqrt((x0[1] - x[1])^2 + (x0[2] - y[1])^2)
                R = sqrt((x[2] - x[1])^2 + (y[2] - y[1])^2)
                di2 = dat[2] - dat[1]
                dat0 = dat[1] + (di1 * di2/R)
            }
            else {
                x0 <- c(x[1], y[1])
                dat0 <- dat[1]
            }
            n <- length(x)
            nn <- nnew * n
            toto <- .C("discretrajr", as.double(x), as.double(y),
                       as.double(dat), double(nn), double(nn), as.integer(n),
                       as.integer(nn), double(nn), as.double(x0[1]), as.double(x0[2]),
                       as.double(u), as.double(dat0), integer(1), PACKAGE = "adehabitatLT")
            neff <- toto[[13]] - 1
            if (neff >= (nn - 1))
                stop("too small rediscretization step length. Try to increase \"nnew\"")
            x <- toto[[4]][1:neff]
            y <- toto[[5]][1:neff]
            dat <- toto[[8]][1:neff]
            class(dat) <- c("POSIXt", "POSIXct")
            attr(dat, "tzone") <- attr(l[[1]]$date, "tzone")

            if (addbit) {
                x <- c(x, l[[1]]$x[length(l[[1]]$x)])
                y <- c(y, l[[1]]$y[length(l[[1]]$y)])
                dat <- c(dat, l[[1]]$date[length(l[[1]]$date)])
            }
            opt <- options(warn = -1)
            if (!attr(ml, "typeII")) {
                nl <- as.ltraj(data.frame(x, y), id = attr(l[[1]],
                                                 "id"), burst = paste(attr(l[[1]], "burst"), ".R",
                                                        u, sep = ""), typeII = FALSE)
            }
            else {
                nl <- as.ltraj(data.frame(x, y), id = attr(l[[1]],
                                                 "id"), date = dat, burst = paste(attr(l[[1]],
                                                                    "burst"), ".R", u, sep = ""))
            }
            nl[[1]]$rel.ang[is.na(nl[[1]]$rel.ang)] <- 0
            class(nl) <- "ltraj"
            options(opt)
            return(nl)
        }
        nl <- do.call("c.ltraj", lapply(burst, foo))
        return(nl)
    } else {
        nl <- do.call("c.ltraj", lapply(1:length(l), function(i) {
            oo <- ld(l[i])[,1:3]
            ii <- as.data.frame(.Call("redistime", oo, as.double(u), as.double(samplex0)), PACKAGE="adehabitatLT")
            df <- ii[ii[,3]>0,]
            da <- df[,3]
            class(da) <- c("POSIXct","POSIXt")
            attr(da, "tzone") <- attr(oo[,3], "tzone")
            as.ltraj(df[,1:2], da, id=id(l)[i], burst=burst(l)[i], typeII=attr(l, "typeII"))
        }))
        return(nl)
    }
}


