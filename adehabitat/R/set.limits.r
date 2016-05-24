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
    dat <- as.POSIXlt(dat, tz)
    class(dat) <- "list"

    napo <- c("sec","min","hour","mday","mon","year","wday","yday")
    nbpo <- c("S","M","H","d","m","Y","w","j")

    amodif <- napo[nbpo%in%ordre]
    reftraj <- do.call("paste", dat[amodif])

    res <- lapply(ltraj, function(x) {
        da <- as.POSIXlt(x$date, tz)
        dact <- as.numeric(as.POSIXct(da, tz))
        class(da) <- "list"

        foctraj <- do.call("paste", da[amodif])
        deds <- (reftraj%in%foctraj)
        re <- x

        if (any(!deds)) {
            avt <- cumsum(deds)
            avt <- length(avt[avt==0])
            if (avt > 0) {
                av <- as.data.frame(matrix(NA,ncol=ncol(x), nrow=avt))
                names(av) <- names(x)
                re <- rbind(av,re)
                foctraj <- c(reftraj[1:avt],foctraj)
                dde <- dact[1]
                aj <- seq(dde,dde-avt*dt, by=-dt)
                aj <- aj[(length(aj)):2]
                dact <- c(aj, dact)
            }

            deds <- deds[(length(deds)):1]
            apr <- cumsum(deds)
            apr <- length(apr[apr==0])

            if (apr > 0) {
                ap <- as.data.frame(matrix(NA,ncol=ncol(x), nrow=apr))
                names(ap) <- names(x)
                re <- rbind(re,ap)
                foctraj <- c(foctraj,
                             reftraj[((length(reftraj)-
                                       apr+1):length(reftraj))])
                fi <- dact[length(dact)]
                aj <- seq(fi, fi+dt*apr, by=dt)[-1]
                dact <- c(dact,aj)
            }
        }

        ## le contraire
        deds <- (foctraj%in%reftraj)
        if (any(!deds)) {
            avt <- cumsum(deds)
            avt <- length(avt[avt==0])
            if (avt > 0) {
                re <- re[-c(1:avt),]
                dact <- dact[-c(1:avt)]
            }
            deds <- deds[(length(deds)):1]
            apr <- cumsum(deds)
            apr <- length(apr[apr==0])

            if (apr > 0) {
                re <- re[-c((nrow(re)-apr+1):nrow(re)),]
                dact <- dact[-c((length(dact)-apr+1):length(dact))]
            }
        }
        class(dact) <- c("POSIXct", "POSIXt")
        re$date <- dact
        attr(re,"id") <- attr(x,"id")
        attr(re,"burst") <- attr(x,"burst")
        return(re)
    })
    class(res) <- c("ltraj","list")
    attr(res, "typeII") <- TRUE
    attr(res, "regular") <- is.regular(res)
    res <- rec(res,...)
    return(res)
}
