### to do: add tangent line, parabola, cubic

TkSpline <- function(x, y, method='natural', snap.to.x=FALSE, digits=4,
                     col=c('blue','#009900','red','black'),
                     xlab=deparse(substitute(x)), ylab=deparse(substitute(y)),
                     hscale=1.5, vscale=1.5, wait=TRUE,
                     ...) {
    if( !requireNamespace('tkrplot', quietly = TRUE) ) stop('This function depends on the tkrplot package being available')

    snap.x <- tcltk::tclVar()
    tcltk::tclvalue(snap.x) <- ifelse(snap.to.x,"T","F")
    d1 <- tcltk::tclVar()
    d2 <- tcltk::tclVar()
    d3 <- tcltk::tclVar()
    tcltk::tclvalue(d1) <- 'F'
    tcltk::tclvalue(d2) <- 'F'
    tcltk::tclvalue(d3) <- 'F'

    xxx <- as.numeric(x)
    ax <- (min(x)+max(x))/2

    sf <- splinefun(x, y, method=method)

    ay <- sf(ax)
    ay1 <- sf(ax, 1)
    ay2 <- sf(ax, 2)
    ay3 <- sf(ax, 3)

    yy <- c(ay,ay1,ay2,ay3)

    txtvar <- tcltk::tclVar()
    tcltk::tclvalue(txtvar) <- " \n \n \n "

    first <- TRUE
    ul <- ur <- 0

    replot <- function() {
        par(mar=c(5,4,2,2)+0.1)
        plot(x,y, xlab=xlab, ylab=ylab, ...)
        u <- par('usr')
        curve(sf(x), from=u[1], to=u[2],add=TRUE)
        lines( c(ax,ax,u[1]), c(u[3], ay, ay), col=col )
        mtext( format( ax, digits=digits), side=3, at=ax, line=1, col=col[1])
        mtext( format( ay, digits=digits), side=4, at=ay, line=1, col=col[1])

        if(as.logical(tcltk::tclvalue(d1))) {
            curve( ay+(x-ax)*yy[2], from=u[1], to=u[2], add=TRUE, col=col[2])
        }
        if(as.logical(tcltk::tclvalue(d2))) {
            curve( ay+(x-ax)*yy[2]+((x-ax)^2)*yy[3],
                  from=u[1], to=u[2], add=TRUE, col=col[3])
        }
        if(as.logical(tcltk::tclvalue(d3))) {
            curve( ay+(x-ax)*yy[2]+((x-ax)^2)*yy[3]+((x-ax)^3)*yy[4],
                  from=u[1], to=u[2], add=TRUE, col=col[4])
        }

        tcltk::tclvalue(txtvar) <<- paste( c('y:   ','d1:   ','d2:   ','d3:   '),
                                   format( yy, digits=digits ), collapse='\n')

        if(first) {
            first <<- FALSE
#            tmp <- cnvrt.coords(c(0,1),c(0,1), input='dev')$usr
            tmpx <- grconvertX(c(0,1), from='ndc')
            ul <<- tmpx[1]
            ur <<- tmpx[2]

        }
    }

    tt <- tcltk::tktoplevel()
    tcltk::tkwm.title(tt, "TkSpline")

    img <- tkrplot::tkrplot(tt, replot, vscale=vscale, hscale=hscale)
    tcltk::tkpack(img, side='top')

    tcltk::tkpack(fr <- tcltk::tkframe(tt), side='left')
    tcltk::tkpack(tcltk::tkcheckbutton(fr,variable=d1, onvalue="T", offvalue="F",
                         text="Show d1", command=function()tkrplot::tkrreplot(img)),
           side='top')
    tcltk::tkpack(tcltk::tkcheckbutton(fr,variable=d2, onvalue="T", offvalue="F",
                         text="Show d2", command=function()tkrplot::tkrreplot(img)),
           side='top')
    tcltk::tkpack(tcltk::tkcheckbutton(fr,variable=d3, onvalue="T", offvalue="F",
                         text="Show d3", command=function()tkrplot::tkrreplot(img)),
           side='top')



    tcltk::tkpack(tcltk::tklabel(tt, textvariable=txtvar), side='top')

    tcltk::tkpack(tcltk::tkcheckbutton(fr,variable=snap.x, onvalue="T", offvalue="F",
                         text="Snap to points"),
           side='top')

    tcltk::tkpack(tcltk::tkbutton(tt,text='Quit', command=function() tcltk::tkdestroy(tt)),
           side='right')

    md <- FALSE
    iw <- as.numeric(tcltk::tcl('image','width',tcltk::tkcget(img,'-image')))
    ih <- as.numeric(tcltk::tcl('image','height',tcltk::tkcget(img,'-image')))

    ccx <- ccy <- 0
    ci <- 0

    mouse.move <- function(x,y) {
        if(md) {
            tx <- (as.numeric(x)-1)/iw
            ccx <<- tx*ur + (1-tx)*ul
            if(as.logical(tcltk::tclvalue(snap.x))) {
                ccx <<- xxx[ which.min( abs(ccx-xxx) ) ]
            }

            ax <<- ccx
            ccy <<- sf(ccx)
            yy <<- c( ccy, sf(ccx,1), sf(ccx,2), sf(ccx,3) )
            ay <<- ccy

            tkrplot::tkrreplot(img)
        }
    }

    mouse.down <- function(x,y) {
        md <<- TRUE
        mouse.move(x,y)
    }

    mouse.up <- function(x,y) {
        md <<- FALSE
    }

    tcltk::tkbind(img, '<Motion>', mouse.move)
    tcltk::tkbind(img, '<ButtonPress-1>', mouse.down)
    tcltk::tkbind(img, '<ButtonRelease-1>', mouse.up)

    if(wait) {
        tcltk::tkwait.window(tt)
        out <- list( x=ccx, y=yy )
    } else {
        out <- NULL
    }

    invisible(out)
}
