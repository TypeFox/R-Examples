.mat <- function(ref, obs)
{
    obscons <- obs
    ## beginning after
    rea <- sapply(1:length(ref), function(i) {
        if (i>1) {
            rr <- ref[-c(1:(i-1))]
        } else {
            rr <- ref
        }
        sum(rr[1:min(c(length(rr), length(obs)))]==obs[1:min(c(length(rr), length(obs)))])
    })

    ## beginning before
    reb <- sapply(1:length(obs), function(i) {
        if (i>1) {
            rr <- obs[-c(1:(i-1))]
        } else {
            rr <- obs
        }
        sum(rr[1:min(c(length(rr), length(ref)))]==ref[1:min(c(length(rr), length(ref)))])
    })

    b <- which.max(c(max(rea), max(reb)))
    if (b==1) {
        rajnadeb <- which.max(rea)-1
        indicedebobs <- 1
        if (rajnadeb>0)
            obs <- c(ref[1:rajnadeb], obs)
    } else {
        rajnadeb <- 0
        indicedebobs <- which.max(reb)
        obs <- obs[indicedebobs:length(obs)]
    }

    ## Et pour la fin:
    rajnafin <- max(length(ref)-length(obs), 0)
    indicefinobs <- length(obscons)-max(length(obs)-length(ref),0)


    return(list(rajna=c(rajnadeb, rajnafin),
                indice=c(indicedebobs,indicefinobs)))
}




set.limits <- function(ltraj, begin, dur, pattern,
                       units=c("sec", "min", "hour", "day"),
                       tz = "", ...)
{
    if (!inherits(ltraj,"ltraj"))
        stop("ltraj should be of class \"ltraj\"")
    if (!is.regular(ltraj))
        stop("ltraj should be regular")
    units <- match.arg(units)
    dur <- .convtime(dur, units)

    res <- strsplit(pattern,"")[[1]]
    sep <- res[c(1,diff(as.numeric(res=="%")))==0][1]
    ordre <- unlist(lapply(split(res,cumsum(as.numeric(res=="%"))),
                           function(x) x[2]))
    dt <- ltraj[[1]]$dt[1]
    beginP <- strptime(begin,pattern, tz=tz)
    beg <- as.numeric(as.POSIXct(beginP, tz))
    en <- beg+dur
    if ((en-beg)%%dt!=0)
        stop("uncorrect duration")
    dat <- seq(beg, en, by=dt)
    class(dat) <- c("POSIXct", "POSIXt")
    datrefc <- dat
    dat <- as.POSIXlt(dat, tz)
    class(dat) <- "list"

    napo <- c("sec","min","hour","mday","mon","year","wday","yday")
    nbpo <- c("S","M","H","d","m","Y","w","j")

    amodif <- napo[nbpo%in%ordre]
    reftraj <- do.call("paste", dat[amodif])

    res <- lapply(ltraj, function(x) {
        da <- as.POSIXlt(x$date, tz)
        infol <- attr(x, "infolocs")
        dact <- as.numeric(as.POSIXct(da, tz))
        class(da) <- "list"
        foctraj <- do.call("paste", da[amodif])
        licom <- .mat(reftraj, foctraj)

        ## update the traj
        arajdeb <- as.data.frame(matrix(NA, ncol=ncol(x), nrow=licom$rajna[1]))
        arajfin <- as.data.frame(matrix(NA, ncol=ncol(x), nrow=licom$rajna[2]))
        names(arajdeb) <- names(arajfin) <- names(x)
        re <- rbind(arajdeb, x[licom$indice[1]:licom$indice[2],], arajfin)
        re$date <- datrefc
        ## update the infol
        if (!is.null(infol)) {
            arajdeb <- as.data.frame(matrix(NA, ncol=ncol(infol), nrow=licom$rajna[1]))
            arajfin <- as.data.frame(matrix(NA, ncol=ncol(infol), nrow=licom$rajna[2]))
            names(arajdeb) <- names(arajfin) <- names(infol)
            infol2 <- rbind(arajdeb, infol[licom$indice[1]:licom$indice[2],, drop=FALSE], arajfin)
        }
        attr(re,"id") <- attr(x,"id")
        attr(re,"burst") <- attr(x,"burst")
        if (!is.null(infol)) {
            attr(re,"infolocs") <- infol2
        }
        return(re)
    })
    class(res) <- c("ltraj","list")
    attr(res, "typeII") <- TRUE
    attr(res, "regular") <- TRUE
    res <- rec(res,...)
    return(res)
}
