###  demo(tkdensity) ## is at
### /u/maechler/R/D/r-devel/Linux-inst/library/tcltk/demo/tkdensity.R

tkdensity <- function(y, n = 1024, log.bw = TRUE, showvalue = TRUE,
                      xlim = NULL, do.rug = size < 1000, kernels = NULL,
                      from.f = if(log.bw) -2   else 1/1000,
                      to.f   = if(log.bw) +2.2 else 2, col = 2)
{
    ## Purpose: as density() but with  scrollbar - bandwidth selection
    ## -----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 8 Nov 2000, 19:00

    requireNamespace("tcltk") || stop("tcltk support is absent")
    tclVar <- tcltk::tclVar
    tclvalue <- tcltk::tclvalue
    tkframe <- tcltk::tkframe
    tkpack <- tcltk::tkpack
    tklabel <- tcltk::tklabel
    tkscale <- tcltk::tkscale
    nbw <- xZ <- xM <- NA_real_ # so '<<-' keeps them here

    dFun <- density.default
    all.kerns <- eval(formals(dFun)$kernel)
    kernels <-
        if(is.null(kernels)) all.kerns
        else match.arg(kernels, all.kerns, several.ok = TRUE)
    ynam <- deparse(substitute(y))
    size <- length(y)
    sd.y <- sqrt(var(y))

    ## Use Silverman's  Rule of Thumb initially :
    hi <- sd.y
    if ((lo <- min(hi, IQR(y)/1.34)) == 0)
        (lo <- hi) || (lo <- abs(y[1])) || (lo <- 1)
    bw <- bw0 <- 0.9 * lo * size^(-0.2)
    if(log.bw) lbw <- lbw0 <- log10(bw0)

    ry <- range(y)
    xlim <- if(is.null(xlim)) ry + c(-2,2)* bw0 else as.numeric(xlim)
    xlmid <- xm0 <- mean(xlim)
    xr0 <- diff(xlim)

    ## Initialize Tcl variables:

    xZoom <- tclVar(100)# %
    xlmid <- tclVar(xlmid)

    if(log.bw)
        Lbw <- tclVar(log10(bw))
    else
        bw <- tclVar(bw)

    kernel <- tclVar("gaussian")

    ## Tvar <- function(v) as.numeric(tclvalue(substitute(v)))

    replot <- function(...) {
        if (is.null(y)) return() # too early...

        b <- if(log.bw) 10 ^ (lbw <<- as.numeric(tclvalue(Lbw))) else
                              nbw <<- as.numeric(tclvalue(bw))
        ##Dbg cat("b = ", formatC(b),"\n")
        k <- tclvalue(kernel) # *is* char
        ##Dbg cat("tclvalue(kernel)"); str(k)

        xZ <<- as.numeric(tclvalue(xZoom))
        xM <<- as.numeric(tclvalue(xlmid))

	## "codetools, please do believe that we do use 'b', 'k', 'xlim' !":
	if(0 > 1)
	    b <- xlim + b + k
        xr.half <- (xr0 / 2) * 100 / xZ
        xlim <- xM + c(-xr.half, xr.half)
        eval(substitute(plot(density(y, bw = b, kernel = k, n = n),
                             main =  paste("density(",ynam,
                             ", bw = ",format(b, dig = 3),
                             ", kernel = \"", k, "\")", sep=""),
                             xlim = xlim, col = col)))
        if(do.rug) rug(y) ## points(y,rep(0,size), col = 3)
    }

    replot.maybe <- function(...)
        if ((log.bw  && !identical(lbw, as.numeric(tclvalue(Lbw)))) ||
	    (!log.bw && !identical(nbw, as.numeric(tclvalue(bw)))) ||
	    !identical(xZ, as.numeric(tclvalue(xZoom))) ||
	    !identical(xM, as.numeric(tclvalue(xlmid)))
	    )
            replot()

    base <- tcltk::tktoplevel()
    tcltk::tkwm.title(base, paste("Tk Density(",ynam,")"))

    base.frame <- tkframe(base, borderwidth = 2)
    bw.frame   <- tkframe(base.frame, relief = "groove", borderwidth = 3)
    kern.frame <- tkframe(base.frame, relief = "groove", borderwidth = 2)

    x.frame    <- tkframe(base.frame)
    xr.frame   <- tkframe(x.frame)
    xmid.frame <- tkframe(x.frame)
    tkpack(xr.frame, xmid.frame, side = "left", anchor = "s")

    q.but <- tcltk::tkbutton(base, text = "Quit", command = function() {
	par(op) ## see par() below !
	tcltk::tkdestroy(base) })
    tkpack(base.frame,
           bw.frame, kern.frame,
           x.frame,
           q.but)

    ## Bandwith Frame :
    tkpack(tklabel (bw.frame,
                    text = if(log.bw)"log10(Bandwidth)" else "Bandwidth"))
    tkpack(tkscale (bw.frame, command = replot.maybe,
                    from = if(log.bw) lbw0 + (from.f) else bw0 * from.f,
                    to   = if(log.bw) lbw0 + (to.f)   else bw0 * to.f,
                    showvalue = showvalue,
                    variable = if(log.bw) Lbw else bw,
                    resolution = if(log.bw) lbw0/20 else bw0/4 * from.f,
                    length = 200,
                    orient = "horiz"))

    ## Kernel Frame :
    tkpack(tklabel(kern.frame, text = "Kernel"))
    for (k.name in kernels)
        tkpack(tcltk::tkradiobutton(kern.frame, command = replot,
                                    text = k.name, value = k.name, variable=kernel),
               anchor = "w")

    ## [x Zoom] Frame :
    tkpack(tklabel (xr.frame, text = "x zoom [%]"))
    tkpack(tkscale (xr.frame, command = replot.maybe,
                    from = 5,# = 1/20
                    to   = 500,# = * 5
                    showvalue = TRUE, variable = xZoom,
                    length = 80, orient = "horiz"))

    ## [x Pan] Frame :
    tkpack(tklabel (xmid.frame, text = "x pan"))
    tkpack(tkscale (xmid.frame, command = replot.maybe,
                    from = xm0 - xr0,
                    to   = xm0 + xr0,
                    showvalue = FALSE, variable = xlmid,
                    resolution = xr0/2000,
                    length = 80, orient = "horiz"))


    if((op <- par("ask")) || prod(par("mfrow")) > 1)
        op <- par(ask = FALSE, mfrow = c(1,1))
    ## on.exit(par(op)) is *NOT* sufficient; do it only when quitting tk !!

    ##Dbg cat("Before calling `replot()' : \n")
    replot()

    ## Returning doesn't work!!
    ##return(tclvar[c("bw", "kernel")])
}

###---

## tkpack() :
##- .Tcl(.Tcl.args(...)) :
##-  [tcl] unknown or ambiguous option  "":  must be \
## 	-after, -anchor, -before, -expand, -fill, -in,
##      -ipadx, -ipady, -padx, -pady, or -side.
