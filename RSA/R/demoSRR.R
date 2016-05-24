#' @export
demoSRR <- function(a4=0, C=0, bm=0, type="3d", zlim=c(-2, 2), xlim=c(-2, 2), ylim=c(-2, 2), xlab="X", ylab="Y", zlab="Z", points = TRUE, model="full", project=c("PA1", "PA2"), legend=FALSE, coefs=TRUE, ...) {

	type <- match.arg(type, c("interactive", "3d", "contour"))
	type2 <- type
	if (type2 == "interactive") stop("demoRSA only works with type == '3d' or 'contour'!")
		
	if (!requireNamespace("tkrplot", quietly = TRUE)) {
		stop('`tkrplot` package needed for demos. Please install with install.packages("tkrplot")', call. = FALSE)
	}	

	fit <- NULL

	# Define all input components	
	I_a4 <- tclVar(); tclvalue(I_a4) <- a4
	I_C <- tclVar(); tclvalue(I_C) <- C
	I_bm <- tclVar(); tclvalue(I_bm) <- bm
	
	update <- function(...) {
        tkrreplot(img, hscale=1.5, vscale=1.5)
    }

    replot <- function() {
		# read parameters from sliders
		a4 <- as.numeric(tclvalue(I_a4))
		C <- as.numeric(tclvalue(I_C))
		bm <- as.numeric(tclvalue(I_bm))

		# derive regression weights from surface parameters
		x2 <- a4/4
		y2 <- x2*1.001	# slight inequality, so that PA1 and PA2 are shown
		xy <- -a4/2
		
		x <- (bm - C*4*x2)/2
		y <- bm - x
		
		plot(plotRSA(x=x, y=y, x2=x2, y2=y2, xy=xy, zlim=zlim, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, zlab=zlab, demo=TRUE, project=project, legend=legend, coefs=coefs, ...))
    }

	# define framework
    tt <- tktoplevel()
    tkwm.title(tt, "Response surface plot - Shifted Rising Ridge")

    img <- tkrplot(tt, replot, vscale=2, hscale=2)
    tkpack(img, side='left')
	

	# define sliders: polynomial model
	tkpack(tfr <- tkframe(tt, relief='groove', borderwidth=3), side='left')
	tkpack(fr1 <- tkframe(tfr), side='top',fill='x')
	tkpack(fr2 <- tkframe(tfr), side='top',fill='x')
	tkpack(fr4 <- tkframe(tfr), side='top',fill='x')
	I_a4.lab <- tklabel(fr1,text='a4: ')
	I_C.lab <- tklabel(fr2,text='C: ')
	I_bm.lab <- tklabel(fr4,text='bm: ')
		
    tkpack(I_a4.lab, side='left',anchor='s')
	tkpack(tkscale(fr1, variable=I_a4, orient='vertical', command=update, from=2, to=-2, resolution=0.01), side='left')

    tkpack(I_C.lab, side='left',anchor='s')
	tkpack(tkscale(fr2, variable=I_C, orient='horizontal', command=update, from=-1, to=1, resolution=0.01), side='left')
	
    tkpack(I_bm.lab, side='left',anchor='s')
	tkpack(tkscale(fr4, variable=I_bm, orient='vertical', command=update, from=2, to=-2, resolution=0.01), side='left')

    return(invisible(NULL))
}

# hack to please CRAN ...
if(getRversion() >= "2.15.1")  {utils::globalVariables('tclvalue')}