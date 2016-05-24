plot.Surv <- function (x, strata = NULL,
                       fn = c("cum", "surv", "log", "loglog"), 
                       limits = TRUE, conf = 0.95, main = NULL, xlab = NULL,
                       ylab = NULL, xlim = NULL, ylim = NULL, lty = NULL,
                       col = NULL, lty.con = NULL, 
                       col.con = NULL, x.axis = TRUE,
                       printLegend = TRUE, ...) 
{
    fn <- fn[1]
    if (is.logical(printLegend)){
        if (printLegend) {
            if (fn == "surv"){
                where <- "bottomright"
            }else{
                where <- "bottomleft"
            }
        }
    }else{
        where <- printLegend
        printLegend <- TRUE
        if (!(where %in% c("bottomleft", "bottomright", "topleft", "topright",
                           "left", "right", "top", "bottom", "center")))
            stop(paste(where, " is not allowed as a value of 'printLegend'"))
    }
        
    if (!inherits(x, "Surv")) 
        stop("First argument must be of type 'Surv'")
    if (ncol(x) == 3) {
        enter <- x[, 1]
        exit <- x[, 2]
        event <- x[, 3]
        n <- length(exit)
    }
    else {
        exit <- x[, 1]
        n <- length(exit)
        enter <- rep(0, n)
        event <- x[, 2]
    }
    yVal <- function(x) {
        if (fn == "cum") 
            return(cumsum(x))
        if (fn %in% c("log", "loglog")) 
            return(log(cumsum(x)))
        n <- length(x)
        s <- numeric(n)
        s[1] <- 1 - x[1]
        if (n > 1) {
            for (rs in 2:n) {
                s[rs] <- s[rs - 1] * (1 - x[rs])
            }
        }
        return(s)
    }
    n <- length(exit)
    if (is.null(strata)) 
        group <- rep(1, n)
    else group <- strata
    if (is.factor(group)) {
        strata <- levels(group)
    }
    else {
        group <- as.character(group)
        strata <- sort(unique(group))
    }
    out <- is.na(group) # Added 2011-11-30 (if missing values in strata)
    if (sum(out)){
        warning("Missing values in stratum variable. Removed cases.")
        group <- group[!out]
        enter <- enter[!out]
        exit <- exit[!out]
        event <- event[!out]
    }
    noOfGroups <- length(strata)
    if (noOfGroups > 1) 
        limits <- FALSE
    if (is.null(col)) {
        col <- rep(1, noOfGroups)
    }
    if (is.null(lty)) {
        lty <- 1:noOfGroups
    }
    fn <- fn[1]
    times <- list()
    atoms <- list()
    i <- 0
    x.min <- min(enter)
    x.max <- max(exit)
    if (fn == "loglog") {
        x.min <- log(min(exit)/2)
        x.max <- log(x.max)
    }
    if (is.null(xlim)) 
        xlim <- c(x.min, x.max)
    if (is.null(ylim)) {
        y.max <- -1e+103
        y.min <- 1e+103
    }
    for (stratum in strata) {
        i <- i + 1
        atom <- table.events(enter[group == stratum], exit[group == 
            stratum], event[group == stratum])
        times[[i]] <- atom$times
        atoms[[i]] <- atom$events/atom$riskset.sizes
        if (fn %in% c("surv", "cum")) {
            times[[i]] <- c(x.min, times[[i]])
            atoms[[i]] <- c(0, atoms[[i]])
        }
        nat <- length(times[[i]])
        mnat <- max(exit[group == stratum])
        if (times[[i]][nat] < mnat) {
            times[[i]] <- c(times[[i]], mnat)
            atoms[[i]] <- c(atoms[[i]], 0)
        }
        if (fn == "loglog") 
            times[[i]] <- log(times[[i]])
        atoms[[i]] <- yVal(atoms[[i]])
        if (is.null(ylim) && (fn != "surv")) {
            y.max <- max(y.max, atoms[[i]][nat])
            y.min <- min(y.min, atoms[[i]][i])
        }
    }
    if (is.null(ylim)) {
        if (fn == "surv") {
            y.min <- 0
            y.max <- 1
            ylim <- c(0, 1)
        }
        else {
            if (fn == "cum") 
                y.min <- 0
            ylim <- c(y.min, y.max)
        }
    }
    else {
        if (length(ylim) != 2) 
            stop("ylim must have length 2.")
        y.min <- min(ylim)
        y.max <- max(ylim)
    }
    if (is.null(ylab)) 
        ylab <- ""
    if (is.null(xlab)) 
        xlab <- "Duration"
    if (is.null(main)) {
        if (fn == "cum") 
            main <- "Cumulative hazard function"
        else if (fn == "surv") 
            main <- "Survivor function"
        else {
            main <- "Log cumulative hazard function"
            if (fn == "loglog") 
                xlab <- "Log(duration)"
        }
    }
    plot(times[[1]], atoms[[1]], main = main, xlab = xlab, ylab = ylab, 
        xlim = xlim, ylim = ylim, type = "s", lty = lty[1], col = col[1], 
        ...)
    if (noOfGroups > 1) {
        for (st in 2:noOfGroups) {
            lines(times[[st]], atoms[[st]], type = "s", lty = lty[st], 
                col = col[st])
        }
        x <- 0.7 * (x.max - x.min) + x.min
        if (fn == "surv") {
            y <- 0.9 * (y.max - y.min) + y.min
        }
        else if (fn == "cumhaz") {
            x <- x.min + 0.01 * (x.max - x.min)
            y <- 0.9 * (y.max - y.min) + y.min
        }
        else {
            y <- 0.4 * (y.max - y.min) + y.min
        }
        if (printLegend){
            if (fn == "surv"){
                where <- "bottomleft"
            }else{
                where <- "bottomright"
            }
            legend(where, legend = strata,
                   lty = lty[1:noOfGroups], 
                   col = col[1:noOfGroups])
        }
    }
    if (fn %in% c("surv", "cum") & x.axis) 
        abline(h = 0)
    if (limits && (fn == "surv")) {
        q.alpha <- abs(qnorm((1 - conf)/2))
        survived <- (atom$riskset.size - atom$events)
        se <- sqrt(cumsum(atom$events/(atom$riskset.sizes * survived))) /
            cumsum(-log(survived/atom$riskset.sizes))
        se <- c(0, se)
        surv <- atoms[[1]]
        n <- length(se)
        if (length(surv) > n) 
            se <- c(se, se[n])
        upper <- surv^exp(q.alpha * se)
        lower <- surv^exp(-q.alpha * se)
        n <- length(upper)
        if (length(times[[1]]) > n) {
            upper <- c(upper, upper[n])
            lower <- c(lower, upper[n])
        }
        if (is.null(lty.con)) {
            lty.con = 2
        }
        if (is.null(col.con)) {
            col.con = 2
        }
        lines(times[[1]], upper, type = "s", lty = lty.con, col = col.con)
        lines(times[[1]], lower, type = "s", lty = lty.con, col = col.con)
    }
}
