see.smooth.tck <- function(){



local({
    have_ttk <- as.character(tcl("info", "tclversion")) >= "8.5"
    if(have_ttk) {
        tkbutton <- ttkbutton
        tkframe <- ttkframe
        tklabel <- ttklabel
        tkradiobutton <- ttkradiobutton
    }

plot_lines <- function(x, y, xlab, ylab, degree, span, kernel, bandwidth, nknots, spar){
    dev.hold()
    par(cex=1.2)
    plot(x, y, xlab=xlab, ylab=ylab)
    legend("topleft", lty = c(1,2,3), lwd = 2, legend = c("lowess", "kernel", "spline"), bty = "n")
    ##################
    o <- order(x)
    lines(x[o],fitted(loess(y[o] ~ x[o], degree = degree, span = span)[o]), lty = 1, lwd = 2)
    
    ##################
    lines(ksmooth(x, y, kernel = kernel, bandwidth = bandwidth), lty = 2, lwd = 2)
    
    ##################    
    lines(smooth.spline(x, y, cv = NA, nknots = nknots, spar = spar, keep.data = FALSE), lty = 3, lwd = 2)
    dev.flush()
    }
        
     
dialog.ci <- function(){
    grDevices::devAskNewPage(FALSE) 
    tclServiceMode(FALSE)
    base <- tktoplevel()
    tkwm.title(base, "Smoothers")

    spec.frm <- tkframe(base, borderwidth=2)
    left.frm <- tkframe(spec.frm)
    
    X.entry <- tkentry(spec.frm, textvariable = X, width = 15)
    Y.entry <- tkentry(spec.frm, textvariable = Y, width = 15)
    Xlab.entry <- tkentry(spec.frm, textvariable = Xlab, width = 15)
    Ylab.entry <- tkentry(spec.frm, textvariable = Ylab, width = 15)
    
    degree <- tclVar(1)
    span <- tclVar(0.5)
    kernel<- tclVar("normal")
    bw    <- tclVar(5)
    spar <- tclVar(0.5)
    knots <- tclVar(8)
    
    replot <- function(...) {
    X <- parse(text = tclvalue(X))[[1]];X <- eval(X)
    Y <- parse(text = tclvalue(Y))[[1]];Y <- eval(Y)
    Xlab <- tclvalue(Xlab)
    Ylab <- tclvalue(Ylab)
    degree <- as.numeric(tclvalue(degree)); degree <- as.numeric(degree)
    span <- tclvalue(span); span <- as.numeric(span)
    
    kernel <- tclvalue(kernel)
    bandwidth <- tclvalue(bw); bandwidth <- as.numeric(bandwidth)
    spar <- tclvalue(spar); spar <- as.numeric(spar)
    knots <- tclvalue(knots); nknots <-as.numeric(knots)
    plot_lines(x = X, y = Y, xlab = Xlab, ylab = Ylab, degree = degree, span=span, kernel=kernel, bandwidth=bandwidth, nknots = nknots, spar=spar)                            
    
}  
     
    
   
    ## Three left frames:
    #LOWESS
    frame1 <-tkframe(left.frm, relief="groove", borderwidth=2)
    tkpack(tklabel(frame1, text="LOWESS"))
    tkpack(tklabel(frame1, text="Degree"))
    for ( i in c(1,2) ) {
        tmp <- tkradiobutton(frame1, command=replot,
                             text=i,value=i,variable = degree) 
                                     tkpack(tmp, anchor="w")
                                        }
    tkpack(tklabel (frame1, text="Span"))
    tkpack(tkscale(frame1, command = replot, from = 0.1, to = 1.0,
                   showvalue = TRUE, variable = span,
                   resolution=0.05, orient="horiz"))

    
    #Kernel
    frame2 <- tkframe(left.frm, relief="groove", borderwidth=2)
    tkpack(tklabel(frame2, text="Kernel smoother"))
    tkpack(tklabel(frame2, text="Kernel"))
    for (i in c("normal", "box")) {
        tmp <- tkradiobutton(frame2, command=replot,
                             text=i, value=i, variable=kernel)
        tkpack(tmp, anchor="w")
    }
    tkpack(tklabel (frame2, text="Bandwidth"))
    tkpack(tkscale(frame2, command=replot, from = 0.5, to = 500,
                   showvalue = TRUE, variable=bw,
                   resolution=0.5, orient="horiz"))
   
   
    #spline
    frame3 <-tkframe(left.frm, relief="groove", borderwidth=2)
    tkpack(tklabel(frame3, text="Spline"))
     
    tkpack(tklabel (frame3, text="Knots"))
    tkpack(tkscale(frame3, command = replot, from = 1, to = 28,
                   showvalue = TRUE, variable = knots,
                   resolution=1, orient="horiz"))
    tkpack(tklabel (frame3, text="Smoothing parameter"))
    tkpack(tkscale(frame3, command = replot, from = 0.01, to = 1.0,                       
                   showvalue = TRUE, variable = spar,
                   resolution=0.05, orient="horiz"))
   
#######################################################################
   
    tkpack(tklabel(spec.frm, text = "X"), X.entry, anchor = "w", fill = "x")
    tkpack(tklabel(spec.frm, text = "Y"), Y.entry, anchor = "w", fill = "x")
    tkpack(tklabel(spec.frm, text = "X-axis label"), Xlab.entry, anchor = "w", fill = "x")
    tkpack(tklabel(spec.frm, text = "Y-axis label"), Ylab.entry, anchor = "w", fill = "x") 
    tkpack(tklabel(spec.frm, text = ""))
        
       
    tkpack(frame1, frame2, frame3)
    tkpack(left.frm, side = "left", anchor = "n")
    q.but <- tkbutton(base, text = " Exit", command=function() tkdestroy(base))
    
    
    tkpack(spec.frm, q.but)
    tclServiceMode(TRUE)
}
        
        X <- tclVar("c(7.1, 4.0, 5.3, 4.5, 4.1, 3.1, 9.0, 5.6, 3.7, 5.8, 9.5, 9.3, 5.7, 3.0, 4.9, 0.8, 5.9, 3.1, 8.3, 8.6, 2.7, 8.2, 5.2, 1.0, 3.8, 9.6, 5, 1.6, 3.6, 2.5)")
        Y <- tclVar("X + X^2 + X^3")
        Xlab <- tclVar("X")
        Ylab <- tclVar("Y")
        Xl <- parse(text = tclvalue(X))[[1]]; n <- length(Xl)
        dialog.ci()
})
  
}