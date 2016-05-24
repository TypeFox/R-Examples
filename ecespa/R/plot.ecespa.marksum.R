 ##PLOT.ecespa.marksum
 ## Modificaciones respecto a ecespa 1.1-0
 ## ribbon=TRUE, controla si pone leyenda de valores
 ## col=NULL, permite seleccionar la gama de colores
 ## main=NULL, permite cambiar el titulo
 ## xlab, ylab, permite cambiar las etiquetas de los ejes
 
 
 plot.ecespa.marksum <-
function (x, what = "normalized", contour = FALSE, grid = FALSE, ribbon=TRUE,col=NULL ,main=NULL,xlab="",ylab="",...) 
{
    
    if (what == "normalized") {
        cosa <- x$normalized
        what = "normalized mark-sum"
    }
    if (what == "pointsum") {
        cosa <- x$pointsum
        what = "point-sum"
    }
    if (what == "marksum") {
        cosa <- x$marksum
        what = "mark-sum"
    }
    plot(Smooth(setmarks(x$grid.ppp, cosa), ...), main = "", col=col, ribbon=ribbon, xlab=xlab,ylab=ylab)
    if(is.null(main)) title(main = paste(x$dataname, "\n", noquote(what), "measure; R=", x$R)) else title(main=main)
    if (contour == TRUE) 
        contour(Smooth(setmarks(x$grid.ppp, cosa), ...), 
            add = TRUE)
    if (grid == TRUE) 
        plot(setmarks(x$grid.ppp, cosa), add = TRUE)
}
