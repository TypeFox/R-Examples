normplot <- function() {
### Interactive density plots. Based on TCL version by Guido Masarotto

require(tcltk) || stop("tcltk support is absent")
local({
    y <- NULL
    x <- seq(-5,8,0.1)
    xlim <-NULL
    size  <- tclVar(0.05)
    dist  <- tclVar(1)
    kernel<- tclVar("gaussian")
    bw    <- tclVar(1)
    bw.sav <- 1 # in case replot.maybe is called too early

    replot <- function(...) {
        if (is.null(y)) return() # too early...
        bw.sav <<- b <- as.numeric(tclvalue(bw))
        k <- tclvalue(kernel)
        sz <- as.numeric(tclvalue(size))
        thresh = qnorm(sz,lower.tail=F)
        pow = pnorm(thresh,lower.tail=F,mean=b)
        plot(x,dnorm(x, mean=b),type="l",main=paste("Power = ",round(pow,2)," at the ",sz," level; Mean =",b))
        points(x,dnorm(x, mean=0),type="l",col="red")
	lines(x=c(thresh,thresh),y=c(0,0.4),col="blue")
    }

    replot.maybe <- function(...)
    {
        if (as.numeric(tclvalue(bw)) != bw.sav) replot()
    }

    regen <- function(...) {
        if (tclvalue(dist)==1) y<<-rnorm(as.numeric(tclvalue(size)))
        else y<<-rexp(as.numeric(tclvalue(size)))
        xlim <<- range(y) + c(-2,2)
        replot()
    }



    base <- tktoplevel()
    tkwm.title(base, "Options")

    spec.frm <- tkframe(base,borderwidth=2)
    left.frm <- tkframe(spec.frm)
    right.frm <- tkframe(spec.frm)

    ## left frames:
    frame1 <-tkframe(left.frm, relief="groove", borderwidth=2)
    tkpack(tklabel(frame1, text="Significance Level"))
    for ( i in c(0.05,0.01,0.001,0.0001) ) {
        tmp <- tkradiobutton(frame1, command=regen,
                             text=i,value=i,variable=size)
        tkpack(tmp, anchor="w")

    }

    ## right frames:
    frame4 <-tkframe(right.frm, relief="groove", borderwidth=2)
    tkpack(tklabel (frame4, text="Mean under H.alt"))
    tkpack(tkscale(frame4, command=replot.maybe, from=0.0, to=4.00,
                   showvalue=F, variable=bw,
                   resolution=0.05, orient="horiz"))

    tkpack(frame1, fill="x")
    tkpack(frame4, fill="x")
    tkpack(left.frm, right.frm,side="left", anchor="n")

    ## `Bottom frame' (on base):
    q.but <- tkbutton(base,text="Quit",
                      command=function()tkdestroy(base))

    tkpack(spec.frm, q.but)

    regen()
})
}
