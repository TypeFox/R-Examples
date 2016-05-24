p.ts <-
    function(x, nrplots = max(1, min(8, n%/%400)), overlap = nk %/% 16,
             date.x = NULL, do.x.axis = !is.null(date.x), do.x.rug = FALSE,
             ax.format, main.tit = NULL, ylim = NULL, ylab = "", xlab = "Time",
             quiet = FALSE, mgp = c(1.25, .5, 0), ...)
{
    ## Purpose: plot.ts with multi-plots + Auto-Title -- currently all on 1 page
    ## -------------------------------------------------------------------------
    ## Arguments: x      : timeseries [ts,rts,its,cts] or numeric vector
    ##            nrplots: number of sub-plots    [DEFAULT: in {1..8}, ~= n/400]
    ##            overlap: how much should subsequent plots overlap [DEFAULT:..]
    ##
    ## Depends on   mult.fig()
    ##
    ## ---> help page  ?p.ts
    ##
    ## -------------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 1 Jul 1994; 18 Dec 1998.

    if(is.null(main.tit)) main.tit <- paste(deparse(substitute(x)))
    isMat <- is.matrix(x)
    n <- if(isMat) nrow(x) else length(x)
    has.date.x <- !is.null(date.x)
    if(do.x.axis && !has.date.x)
        stop("'do.x.axis' is true, but 'date.x' is NULL")
    if(has.date.x) {
        if(n != length(date.x))
            stop("'date.x' must be date vector of the same length as series")
        if(do.x.axis)
            date.x <- as.POSIXct(date.x)    # work, or give error now
        if(is.unsorted(date.x, na.rm=TRUE)) {
            i <- order(date.x)
            x <- if(isMat) x[i,] else x[i]
            date.x <- date.x[i]
        }
        xaxt <- "n"
    } else xaxt <- par("xaxt")
    if(nrplots == 1) {
        if(has.date.x) {
            plot(date.x, x, ..., ylim = ylim, type = 'l',
                 main = main.tit, xlab = xlab, ylab = ylab, xaxt = xaxt)
            if(do.x.axis) axis.POSIXct(1, x = date.x, format = ax.format)
        }
        else
            plot.ts(x, ..., ylim = ylim,
                    main = main.tit, xlab = xlab, ylab = ylab, xaxt = xaxt)
    }
    else if(nrplots <= 0)
        return(nrplots)
    else {                              # nrplots >= 2 :
        if(n <= 1) stop("`x' must have at least two points!")
        if(!is.ts(x)) x <- as.ts(x)
        ##-  do.dates <- !is.null(class(x)) && class(x) == "cts"
        ##-  if(do.dates) x <- as.rts(x)# dates() as below fails [S+ 3.4]
        ## NB: end() and start() are of length 1 _or_ 2 (!)
        scal <- (end(x) - (t1 <- start(x))) / (n-1)
        nk <- n %/% nrplots
        if(is.null(ylim))
            ylim <- range(pretty(range(x, na.rm = TRUE)))
        ##    --------
        if(!quiet)
            Form <- function(x)
                paste("(",paste(formatC(x, digits=6, width=1), collapse=", "),
                      ")",sep='')
        pp <- mult.fig(mfrow=c(nrplots,1), main = main.tit, quiet= TRUE,
                       mgp = mgp, marP = c(-1,-1,-2,0))
        on.exit(par(pp $ old.par))
        for(i in 1:nrplots) {
            i0  <- as.integer(max(0, (-overlap + (i-1)*nk)-1) )
            in1 <- as.integer(min(n, i*nk + overlap)-1 )
            st <- t1 + scal*i0 ##; if(do.dates) st <- dates(st)
            en <- t1 + scal*in1 ##; if(do.dates) en <- dates(en)
            if(!quiet)
                cat(sprintf("%2d -- start{%d}= %s; end{%d}= %s\n",
                            i, i0,Form(st), in1, Form(en)))
        if(has.date.x) {
            plot(date.x[1+ i0:in1], window(x, start= st, end = en),
                 ..., ylim = ylim, type = 'l',
                 xlab = xlab, ylab = ylab, xaxt = xaxt)
            if(do.x.axis) {
                if(!quiet) {
                    cat("summary(date.x):\n"); print(summary(date.x[1+ i0:in1]))
                }
                axis.POSIXct(1, x = date.x[1+ i0:in1], format = ax.format)
                ## (I've lost my improved version of this which had 'nYrs = 12'

                if(do.x.rug) ## this can be ugly
                    rug(date.x[1+ i0:in1])
            }
        }
        else
            plot(window(x, start= st, end = en), ylim = ylim,
                 xlab = xlab, ylab = ylab, xaxt = xaxt, ...,
                 plot.type= "single")# plot.type : for plot.mts only,
        }
    }
}
