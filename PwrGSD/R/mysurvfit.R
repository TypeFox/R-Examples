Surv <- survival::Surv

"mysurvfit" <- 
function (formula = formula(data), data = parent.frame(), subset, 
    na.action = na.fail) 
{
    # call
    #
    # mysurvfit(cbind(ti, ev)~group, ...)
    #
    # or
    #
    # mysurvfit(cbind(ti, ev*evty)~group, ...)
    #
    # for competing risks
    #
    # finds nelson aalen hazard estimates separately within each level of 'group'
    # can accomodate competing risks (as in second call above)
    #
    # evty = integer, taking values 0, 1, ..., #event types.
    # group = either a factor variable or an integer taking no more than a reasonable number of levels
    #
    m <- .call. <- match.call()
    m[[1]] <- as.name("model.frame")
    m <- eval(m, parent.frame())

    R <- model.extract(m, "response")
    ind.too.small <- (R[, 1] < 1e-10)
    n.too.small <- sum(ind.too.small)
    TOS <- R[, 1]
    Event.0 <- R[, 2]
    isfctr <- is.factor(Event.0)
    if(isfctr) Event.levs <- sort(unique(as.character(Event.0)))[-1]
    Event <- as.numeric(Event.0)
    Event.f <- as.factor(as.character(Event))
    if(!isfctr) Event.levs <- levels(Event.f)[-1]
    nevty <- length(Event.levs)
    Event <- as.numeric(Event.f) - 1

    ntimes <- length(unique(TOS[Event != 0]))
    n <- length(Event)
    Arm.vnm <- as.character(.call.$formula[[3]])
    Arm <- 0 * Event
    Arm.levs <- 0
    nb <- 1
    
    int.only <- (length(m) == 1)
    if (int.only) {
        Arm.vnm <- ""
        Arm.levs <- ""
    }
    if(!int.only)
    {
      isfctr <- is.factor(m[[names(m)[2]]])
      if(isfctr) Arm.levs <- sort(unique(as.character(m[[names(m)[2]]])))
      Arm <- as.numeric(m[[names(m)[2]]])
      Arm.f <- as.factor(as.character(Arm))
      if(!isfctr) Arm.levs <- levels(Arm.f)
      nb <- length(Arm.levs)
      Arm <- as.numeric(Arm.f) - 1
    }
    
    ans <- .C("mysurvfit", TOS = as.double(TOS), Event = as.integer(Event), 
        Arm = as.integer(Arm), pn = as.integer(n), time = double(ntimes), 
        nrisk = integer(ntimes * nb), nevent = integer(ntimes * 
            nevty * nb), pntimes = as.integer(ntimes), pnevtypes = as.integer(nevty), 
        pnblocks = as.integer(nb), PACKAGE = "PwrGSD")
    if (nevty > 1) {
        tbl <- as.data.frame(cbind(ans$time, t(matrix(ans$nrisk, 
            nb, ntimes)), matrix(aperm(array(ans$nevent, c(nevty, 
            nb, ntimes)), c(3, 1, 2)), ntimes, nb * nevty)))
        ty.bl.levs <- c(outer("Ev" %,% Event.levs %,% ".", Arm.vnm %,% 
            Arm.levs, FUN = "%,%"))
        nms <- c("time", "nrisk." %,% Arm.vnm %,% Arm.levs, "nevent." %,% 
            ty.bl.levs)
    }
    if (nevty == 1) {
        tbl <- as.data.frame(cbind(ans$time, t(matrix(ans$nrisk, 
            nb, ntimes)), t(matrix(ans$nevent, nb, ntimes))))
        nms <- c("time", "nrisk." %,% Arm.vnm %,% Arm.levs, "nevent." %,%  Arm.vnm %,% Arm.levs)
    }
    names(tbl) <- nms
    out <- list(call = .call., Table = tbl)
    class(out) <- "blkdcp"
    out
}

"print.blkdcp" <-
function(x, ...)
{
    cat("Call:\n",deparse(x$call),"\n\n")
    cat("Table:\n")
    print(x$Table)
    invisible(x$Table)
}

"summary.blkdcp" <-
function (object, event.name = "Incidence", colors = NULL, ...) 
{
    M <- as.matrix(object$Table)

    ri.ind <- which(regexpr("nrisk", dimnames(M)[[2]])>0)
    ev.ind <- which(regexpr("nevent", dimnames(M)[[2]])>0)

    nb <- length(ri.ind)
    nevty <- length(ev.ind)/nb
    
    ti <- M[, 1]
    n.ti <- length(ti)
    A <- matrix(0, n.ti, nb*nevty)
    for(k in 1:nevty)
      A[,nevty*(0:(nb-1))+k] <- apply(M[, ev.ind[nevty*(0:(nb-1))+k], drop=FALSE]/M[, ri.ind, drop=FALSE], 2, FUN = cumsum)
    cumhaz <- data.frame(ti, A)
    names(cumhaz) <- c("times", dimnames(A)[[2]])
    object$cumhaz <- cumhaz

    pu <- apply(M[, ri.ind, drop=FALSE] * DX(ti), 2, FUN=sum)
    Ev <- apply(M[, ev.ind, drop=FALSE], 2, FUN=sum)
    object$crudehaz <- Ev/rep(pu, each=nevty)

    object$events <- apply(M[, ev.ind, drop=FALSE], 2, FUN=sum)
    object$pu <- apply(M[, ri.ind, drop=FALSE] * DX(ti), 2, FUN=sum)
   
    object
}

"plot.blkdcp" <-
function(x, event.name="Incidence", colors=NULL, lty=NULL,...)
{
    M <- as.matrix(x$Table)
    
    ri.ind <- which(regexpr("nrisk", dimnames(M)[[2]])>0)
    ev.ind <- which(regexpr("nevent", dimnames(M)[[2]])>0)

    nb <- length(ri.ind)
    nevty <- length(ev.ind)/nb
    
    if (missing(colors)) colors <- 1:nevty
    if (missing(lty)) lty <- ((1:nevty)-1)%%6 + 1

    ti <- M[,1]
    n.ti <- length(ti)
    
    A <- matrix(0, n.ti, nb*nevty)
    for(k in 1:nevty)
      A[,nb*(k-1)+(1:nb)] <- apply(M[, ev.ind[nb*(k-1)+(1:nb)], drop=FALSE]/M[, ri.ind, drop=FALSE], 2, FUN = cumsum)

    rng.ti <- range(ti)
    rng.A <- range(quantile(A, 0.01*(0:100))[-101])
    plot(rng.ti, rng.A, type="n", xlab="Time on Study", ylab="Cummulative " %,% event.name)
    for(k in 1:nb)
      for(j in 1:nevty) lines(ti, A[,nb*(k-1) + j], type="s", col=colors[k], lty=lty[j])
    invisible(x)
}
