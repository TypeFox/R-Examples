# ---------------------------------------
# Author: Andreas Alfons
#         Vienna University of Technology
# ---------------------------------------

TKRdevice <- function(title) {
    tt <- tktoplevel()
    # set title
    if(missing(title)) title <- " "
    tkwm.title(tt,title)
    # get resolution
    screenwidth <- as.numeric(tclvalue(tkwinfo("screenwidth", tt)))
    screenheight <- as.numeric(tclvalue(tkwinfo("screenheight", tt)))
    # define canvas and frame
    canvasFrame <- tkframe(tt)
    canvas <- tkcanvas(canvasFrame, highlightthickness=0)
    frame <- tkframe(canvas)
    # return object
    res <- list(top=tt, canvasFrame=canvasFrame, canvas=canvas, frame=frame, 
        resolution=c(screenwidth, screenheight))
    class(res) <- "TKRdevice"
    res
}

finish <- function(window, ...) UseMethod("finish")

finish.TKRdevice <- function(window, img, label = NULL, ...) {
    # a label for the device
    if(length(label)) tkpack(tklabel(window$top, text=label), side="top")
    # display image in canvas
    tkcreate(window$canvas, "window", 0, 0, anchor="nw", window=window$frame)
    tkpack(img, fill="both", expand=TRUE)
    .Tcl("update idletasks")
    # horizontal scrollbar
    hscr <- tkscrollbar(window$canvasFrame, repeatinterval=5, 
        orient="horizontal", command=function(...)tkxview(window$canvas, ...))
    # vertical scrollbar
    vscr <- tkscrollbar(window$canvasFrame, repeatinterval=5, 
        orient="vertical", command=function(...) tkyview(window$canvas, ...))
    .Tcl("update idletasks")
    # pack elements
    tkgrid(window$canvas, row=0, column=0, sticky="news")
    tkgrid(hscr, row=1, column=0, sticky="ew")
    tkgrid(vscr, row=0, column=1, sticky="ns")
    .Tcl("update idletasks")
    # make sure canvas and scrollbars expand or shrink when device is resized
    tkgrid.rowconfigure(window$canvasFrame, 0, "-weight", 1)
    tkgrid.columnconfigure(window$canvasFrame, 0, "-weight", 1)
    tkpack(window$canvasFrame, side="bottom", fill="both", expand=TRUE)
    .Tcl("update idletasks")
    # set canvas size
    framewidth <- as.numeric(tclvalue(tkwinfo("reqwidth", window$frame)))
    frameheight <- as.numeric(tclvalue(tkwinfo("reqheight", window$frame)))
    tkconfigure(window$canvas, width=min(framewidth, window$resolution[1]-100),
        height=min(frameheight, window$resolution[2]-100))
    .Tcl("update idletasks")
    # set scroll region
    chscrreg <- paste("0 0 ", framewidth, frameheight)
    tkconfigure(window$canvas, 
        xscrollcommand=function(...) tkset(hscr,...), scrollregion=chscrreg)
    tkconfigure(window$canvas, 
        yscrollcommand=function(...) tkset(vscr,...), scrollregion=chscrreg)
    .Tcl("update idletasks")
    # limit device size
    devwidth <- as.numeric(tclvalue(tkwinfo("reqwidth", window$top)))
    devheight <- as.numeric(tclvalue(tkwinfo("reqheight", window$top)))
    tkwm.maxsize(window$top, min(devwidth, window$resolution[1]-100), 
        min(devheight, window$resolution-100))
    tkwm.minsize(window$top, 50, 50)
}

