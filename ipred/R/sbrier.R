# $Id: sbrier.R,v 1.5 2009/03/27 16:18:38 hothorn Exp $

getsurv <- function(obj, times)
{
    # get the survival probability for times from KM curve `obj'

    if (!inherits(obj, "survfit")) stop("obj is not of class survfit")
    # <FIXME: methods may have problems with that>
    class(obj) <- NULL
    # </FIXME>
    lt <- length(times)
    nsurv <- times

    # if the times are the same, return the km-curve

    if(length(times) == length(obj$time)) {
        if (all(times == obj$time)) return(obj$surv)
    }

    # otherwise get the km-value for every element of times separatly
 
    inside <- times %in% obj$time
    for (i in (1:lt)) {
        if (inside[i])
            nsurv[i] <- obj$surv[obj$time == times[i]]
        else  {
            less <- obj$time[obj$time < times[i]]
            if (length(less) == 0) 
                nsurv[i] <- 1
            else 
                nsurv[i] <- obj$surv[obj$time == max(less)]
        }
    }
    nsurv
}

sbrier <- function(obj, pred, btime = range(obj[,1]))
{
    if(!inherits(obj, "Surv"))
        stop("obj is not of class Surv")

    # check for right censoring

    # <FIXME>
    class(obj) <- NULL
    # </FIXME>
    if (attr(obj, "type") != "right")
        stop("only right-censoring allowed")
    N <- nrow(obj)	

    # get the times and censoring of the data, order them with resp. to time

    time <- obj[,1]
    ot <- order(time)
    cens <- obj[ot,2]
    time <- time[ot]

    # get the times to compute the (integrated) Brier score over

    if (is.null(btime)) stop("btime not given")
    if (length(btime) < 1) stop("btime not given")

    if (length(btime) == 2) {
        if (btime[1] < min(time)) warning("btime[1] is smaller than min(time)")
        if (btime[2] > max(time)) warning("btime[2] is larger than max(time)")
        btime <- time[time >= btime[1] & time <=
                                          btime[2]]
    }

    ptype <- class(pred)
    # <begin> S3 workaround
    if (is.null(ptype)) {
      if (is.vector(pred)) ptype <- "vector"
      if (is.list(pred)) ptype <- "list"
    }
    # <end>
    if (ptype == "numeric" && is.vector(pred)) ptype <- "vector"

    survs <- NULL
    switch(ptype, survfit = {
        survs <- getsurv(pred, btime)
        survs <- matrix(rep(survs, N), nrow=length(btime))
    }, list = {
        if (!inherits(pred[[1]], "survfit")) stop("pred is not a list of survfit objects") 
        if (length(pred) != N) stop("pred must be of length(time)")
        pred <- pred[ot]
        survs <-  matrix(unlist(lapply(pred, getsurv, times = btime)),
                                nrow=length(btime), ncol=N)
    }, vector = {
        if (length(pred) != N) stop("pred must be of length(time)")
        if (length(btime) != 1) stop("cannot compute integrated Brier score with pred")
        survs <- pred[ot]
    }, matrix = {
        # <FIXME>
        if (all(dim(pred) == c(length(btime), N)))
            survs <- pred[,ot]
        else
            stop("wrong dimensions of pred")
        # </FIXME>
    })
    if (is.null(survs)) stop("unknown type of pred")

    # reverse Kaplan-Meier: estimate censoring distribution

    ### deal with ties
    hatcdist <- prodlim(Surv(time, cens) ~ 1,reverse = TRUE)
    csurv <- predict(hatcdist, times = time, type = "surv")
    csurv[csurv == 0] <- Inf
    # hatcdist <- survfit(Surv(time, 1 - cens) ~ 1)
    # csurv <- getsurv(hatcdist, time)
    # csurv[csurv == 0] <- Inf

    bsc <- rep(0, length(btime))
    
    # compute Lebesque-integrated Brier score

    if (length(btime) > 1) {
        for (j in 1:length(btime)) {
            help1 <- as.integer(time <= btime[j] & cens == 1)
            help2 <- as.integer(time > btime[j])
            bsc[j] <-  mean((0 - survs[j,])^2*help1*(1/csurv) +
                            (1-survs[j,])^2*help2*(1/csurv[j]))
        }

        ### apply trapezoid rule
        idx <- 2:length(btime)
        RET <- diff(btime) %*% ((bsc[idx - 1] + bsc[idx]) / 2)
        RET <- RET / diff(range(btime))

        ### previously was
        #diffs <- c(btime[1], btime[2:length(btime)] -
        #                     btime[1:(length(btime)-1)])
        #RET <- sum(diffs*bsc)/max(btime)
        names(RET) <- "integrated Brier score"
        attr(RET, "time") <- range(btime)

    # compute Brier score at one single time `btime'
 
    } else {
        help1 <- as.integer(time <= btime & cens == 1)
        help2 <- as.integer(time > btime)
        cs <- predict(hatcdist, times=btime, type = "surv")
        ### cs <- getsurv(hatcdist, btime)
        if (cs == 0) cs <- Inf
        RET <-  mean((0 - survs)^2*help1*(1/csurv) +
                     (1-survs)^2*help2*(1/cs))
        names(RET) <- "Brier score"
        attr(RET, "time") <- btime
    }
    RET
}
