#' @title Plots a response surface of a polynomial equation of second degree with interactive controls
#'
#' @description
#' Plots an RSA object, or a response surface with specified parameters, with interactive controls for coefficients.
#'
#' @details
#' No details so far. Just play around with the interface!
#'
#' @aliases demoSRR demoSRRR
#'
#' @export
#' @import tkrplot
#' @import tcltk
#' @param x Either an RSA object (returned by the \code{RSA} function), or the coefficient for the X predictor
#' @param y Y coefficient
#' @param x2 X^2 coefficient
#' @param y2 Y^2 coefficient
#' @param xy XY interaction coefficient
#' @param y3 Y^3 coefficient
#' @param x3 X^3 coefficient
#' @param xy2 XY^2 coefficient
#' @param x2y X^2Y coefficient
#' @param w W coefficient (for (un)constrained absolute difference model)
#' @param wx WX coefficient (for (un)constrained absolute difference model)
#' @param wy WY coefficient (for (un)constrained absolute difference model)
#' @param b0 Intercept
#' @param xlim Limits of the x axis
#' @param ylim Limits of the y axis
#' @param zlim Limits of the z axis
#' @param xlab Label of the x axis
#' @param ylab Label of the y axis
#' @param zlab Label of the z axis

#' @param type \code{3d} for 3d surface plot, \code{contour} for 2d contour plot. Shortcuts (i.e., first letter of string) are sufficient; be careful: "contour" is very slow at the moment
#' @param points A list of parameters which define the appearance of the raw scatter points: show = TRUE: Should the original data points be overplotted? value="raw": Plot the original z value, "predicted": plot the predicted z value. jitter=0: Amount of jitter for the raw data points. cex = .5: multiplication factor for point size. See ?plotRSA for details.
#' @param project Which geatures should be projected on the floor? See ?plotRSA for details.
#' @param model If x is an RSA object: from which model should the response surface be computed?
#' @param extended Show additional controls (not implemented yet)
#' @param ... Other parameters passed through to plot.RSA (e.g., xlab, ylab, zlab, cex, legend)
#'
#'
#' @seealso \code{\link{plotRSA}}, \code{\link{RSA}}
#'
#' @examples
#' # Plot response surfaces from known parameters
#' # example of Edwards (2002), Figure 3
#' \dontrun{
#' demoRSA(x=.314, y=-.118, x2=-.145, y2=-.102, xy=.299, b0=5.628, type="3d")
#' demoRSA(x=.314, y=-.118, x2=-.145, y2=-.102, xy=.299, b0=5.628, legend=FALSE, type="c")
#' }
#'
#' # Plot response surface from an RSA object
#' \dontrun{
#' set.seed(0xBEEF)
#' n <- 300
#' err <- 2
#' x <- rnorm(n, 0, 5)
#' y <- rnorm(n, 0, 5)
#' df <- data.frame(x, y)
#' df <- within(df, {
#' 	diff <- x-y
#' 	absdiff <- abs(x-y)
#' 	SD <- (x-y)^2
#' 	z.diff <- diff + rnorm(n, 0, err)
#' 	z.abs <- absdiff + rnorm(n, 0, err)
#' 	z.sq <- SD + rnorm(n, 0, err)
#' 	z.add <- diff + 0.4*x + rnorm(n, 0, err)
#' 	z.complex <- 0.4*x + - 0.2*x*y + + 0.1*x^2 - 0.03*y^2 + rnorm(n, 0, err)
#' })
#' 
#' r1 <- RSA(z.sq~x*y, df)
#' demoRSA(r1)
#' demoRSA(r1, points=TRUE, model="SQD")
#' }

## TODO: Convert to Shiny app

demoRSA <- function(x=NULL, y=0, x2=0, y2=0, xy=0, w=0, wx=0, wy=0, x3=0, xy2=0, x2y=0, y3=0, b0=0, type="3d", zlim=c(-2, 2), xlim=c(-2, 2), ylim=c(-2, 2), xlab=NULL, ylab=NULL, zlab=NULL, points = TRUE, model="full", project=c("PA1", "PA2"), extended=FALSE, ...) {

	type <- match.arg(type, c("interactive", "3d", "contour"))
	type2 <- type
	if (type2 == "interactive") stop("demoRSA only works with type == '3d' or 'contour'!")
		
	if (!requireNamespace("tkrplot", quietly = TRUE)) {
		stop('`tkrplot` package needed for modeltrees. Please install with install.packages("tkrplot")', call. = FALSE)
	}	
	

	# if model is provided: take its parameters as starting values
	if (is.null(x)) {
		x <- 0
		fit <- NULL
		points <- FALSE
		if (is.null(xlab)) {xlab <- "X"}
		if (is.null(ylab)) {ylab <- "Y"}
		if (is.null(zlab)) {zlab <- "Z"}
			
	} else if (!is.null(x) & !is.null(attr(x, "class"))) {
		if (attr(x, "class") == "RSA") {
			fit <- x
			C <- coef(fit$models[[model]])
			b0.0 <- b0 <- as.numeric(ifelse(is.na(C[paste0(fit$DV, "~1")]), b0, C[paste0(fit$DV, "~1")]))
			x.0 <- x <- as.numeric(ifelse(is.na(C["b1"]), 0, C["b1"]))
			y.0 <- y <- as.numeric(ifelse(is.na(C["b2"]), y, C["b2"]))
			x2.0 <- x2 <- as.numeric(ifelse(is.na(C["b3"]), x2, C["b3"]))
			y2.0 <- y2 <- as.numeric(ifelse(is.na(C["b5"]), y2, C["b5"]))
			xy.0 <- xy <- as.numeric(ifelse(is.na(C["b4"]), xy, C["b4"]))
			w.0 <- w <- as.numeric(ifelse(is.na(C["b6"]), w, C["b6"]))
			wx.0 <- wx <- as.numeric(ifelse(is.na(C["b7"]), wx, C["b7"]))
			wy.0 <- wy <- as.numeric(ifelse(is.na(C["b8"]), wy, C["b8"]))
			
			x3.0 <- x3 <- as.numeric(ifelse(is.na(C["b9"]), x3, C["b9"]))
			xy2.0 <- xy2 <- as.numeric(ifelse(is.na(C["b10"]), xy2, C["b10"]))
			x2y.0 <- x2y <- as.numeric(ifelse(is.na(C["b11"]), x2y, C["b11"]))
			y3.0 <- y3 <- as.numeric(ifelse(is.na(C["b12"]), y3, C["b12"]))
		
			xlim <- c(min(fit$data[, fit$IV1], na.rm=TRUE), max(fit$data[, fit$IV1], na.rm=TRUE))
			ylim <- c(min(fit$data[, fit$IV2], na.rm=TRUE), max(fit$data[, fit$IV2], na.rm=TRUE))
			
			if (is.null(xlab)) xlab <- fit$IV1
			if (is.null(ylab)) ylab <- fit$IV2
			if (is.null(zlab)) zlab <- fit$DV
				
			# expand range by 20% at each end
			xlim[1] <- xlim[1]*ifelse(xlim[1]<0, 1.1, 0.9)
			xlim[2] <- xlim[2]*ifelse(xlim[2]<0, 0.9, 1.1)
			ylim[1] <- ylim[1]*ifelse(ylim[1]<0, 1.1, 0.9)
			ylim[2] <- ylim[2]*ifelse(ylim[2]<0, 0.9, 1.1)
				
			# for the correct visual diagonal: same range for X and Y
			xlim[1] <- ylim[1] <- min(xlim[1], ylim[1])
			xlim[2] <- ylim[2] <- max(xlim[2], ylim[2])
		
			zlim <- c(min(fit$data[, fit$DV], na.rm=TRUE)*0.8, max(fit$data[, fit$DV], na.rm=TRUE)*1.2)
			
			# define the defaults
			if (is.null(points) || (typeof(points) == "logical" && points == TRUE)) {
				points <- list(show=TRUE, value="raw", jitter=0, color="black", cex=.5, out.mark=FALSE)
			}
			if (is.null(points) || (typeof(points) == "logical" && points == FALSE)) {
				points <- list(show=FALSE, value="raw", jitter=0, color="black", cex=.5, out.mark=FALSE)
			}
			if (is.null(points$out.mark)) points$out.mark <- FALSE

			if (points$out.mark == FALSE) {data.used <- fit$data[fit$data$out==FALSE, ]}
			if (points$out.mark == TRUE) {data.used <- fit$data}
			
			points$data <- data.used[, c(fit$IV1, fit$IV2, fit$DV, colnames(fit$data)[which(!colnames(fit$data) %in% c(fit$IV1, fit$IV2, fit$DV))])]
		}
	} else {
		fit <- NULL
		points <- FALSE
		if (is.null(xlab)) {xlab <- "X"}
		if (is.null(ylab)) {ylab <- "Y"}
		if (is.null(zlab)) {zlab <- "Z"}
	}
	

    TYPE <- tclVar(); tclvalue(TYPE) <- "full"
	
	B0 <- tclVar(); tclvalue(B0) <- b0
	X <- tclVar(); tclvalue(X) <- x
	Y <- tclVar(); tclvalue(Y) <- y
	X2 <- tclVar(); tclvalue(X2) <- x2
	Y2 <- tclVar(); tclvalue(Y2) <- y2
	XY <- tclVar(); tclvalue(XY) <- xy
	
	W <- tclVar(); tclvalue(W) <- w
	WX <- tclVar(); tclvalue(WX) <- wx
	WY <- tclVar(); tclvalue(WY) <- wy
	
	# rotation of the 3d-frame
	RX <- tclVar(); tclvalue(RX) <- -45
	RY <- tclVar(); tclvalue(RY) <- 45
	RZ <- tclVar(); tclvalue(RZ) <- 35
	
	# Dummy variables: Shift and rotation
	C <- tclVar(); tclvalue(C) <- 0
	S <- tclVar(); tclvalue(S) <- 1
	
	if (extended==TRUE) {
		X.Y2 <- tclVar(); tclvalue(X.Y2) <- 0
	}
	
	
	setAllBlack <- function() {
		sapply(list(X.lab, Y.lab, X2.lab, Y2.lab, XY.lab, W.lab, WX.lab, WY.lab), tkconfigure, foreground="black")
	}

	update <- function(...) {
		
		# hack to please CRAN ...
		#if(getRversion() >= "2.15.1")  {utils::globalVariables('tclvalue')}
				
		# read parameters from sliders
        type <- as.character(tclvalue(TYPE))
		b0 <- as.numeric(tclvalue(B0))
		x <- as.numeric(tclvalue(X))
		y <- as.numeric(tclvalue(Y))
		x2 <- as.numeric(tclvalue(X2))
		y2 <- as.numeric(tclvalue(Y2))
		xy <- as.numeric(tclvalue(XY))
		w <- as.numeric(tclvalue(W))
		wx <- as.numeric(tclvalue(WX))
		wy <- as.numeric(tclvalue(WY))
		rx <- as.numeric(tclvalue(RX))
		ry <- as.numeric(tclvalue(RY))
		rz <- as.numeric(tclvalue(RZ))
		c <- as.numeric(tclvalue(C))
		s <- as.numeric(tclvalue(S))
		if (extended==TRUE) {x.y2 <- as.numeric(tclvalue(X.Y2))}
        
		setAllBlack()
		
		# set constraints
		if (type == "all") {
		}
		if (type == "poly") {
			tclvalue(W) <- tclvalue(WX) <- tclvalue(WY) <- 0
			sapply(list(W.lab, WX.lab, WY.lab), tkconfigure, foreground="grey40")
		}
		
		if (type == "diff") {
			tclvalue(Y) <- -x
			tclvalue(X2) <- tclvalue(Y2) <- tclvalue(XY) <- tclvalue(W) <- tclvalue(WX) <- tclvalue(WY) <- 0
			sapply(list(X2.lab, Y2.lab, XY.lab, W.lab, WX.lab, WY.lab), tkconfigure, foreground="grey40")
		}
		if (type == "SQD") {
			tclvalue(Y) <- 0
			tclvalue(X) <- 0
			tclvalue(Y2) <- x2
			tclvalue(XY) <- -2*x2
			tclvalue(W) <- tclvalue(WX) <- tclvalue(WY) <- 0
			sapply(list(X.lab, Y.lab, Y2.lab, XY.lab, W.lab, WX.lab, WY.lab), tkconfigure, foreground="grey40")
		}
		if (type == "sq.shift") {
			tclvalue(Y) <- -x
			tclvalue(Y2) <- x2
			tclvalue(XY) <- -2*x2
			if (x2 != 0) {
				tclvalue(B0) <- x^2 / (4*x2)				
				tclvalue(C) <- x/(2*x2)
			}
			tclvalue(W) <- tclvalue(WX) <- tclvalue(WY) <- 0
			sapply(list(Y.lab, Y2.lab, XY.lab, W.lab, WX.lab, WY.lab), tkconfigure, foreground="grey40")
		}
		if (type == "sq.rot") {
			#tclvalue(X) <- 2*c*s*y2
			#tclvalue(Y) <- -2*c*y2
			#tclvalue(X2) <- s^2*y2
			#tclvalue(XY) <- -2*s*y2
			
			if (y2 != 0) {
				tclvalue(X2) <- (xy^2) / (4*y2)
				x <- tclvalue(X) <- (y*xy)/(2*y2)
			}
			
			if (y2 != 0 & y != 0) {
				tclvalue(B0) <- y^2 / (4*y2)
				tclvalue(C) <- -0.5*(y/y2)
				tclvalue(S) <- -(x/y)
			}
			tclvalue(W) <- tclvalue(WX) <- tclvalue(WY) <- 0
			sapply(list(X.lab, X2.lab, W.lab, WX.lab, WY.lab), tkconfigure, foreground="grey40")
		}
		if (type == "IA") {
			tclvalue(X2) <- 0
			tclvalue(Y2) <- 0
			tclvalue(W) <- tclvalue(WX) <- tclvalue(WY) <- 0
			sapply(list(X2.lab, Y2.lab, W.lab, WX.lab, WY.lab), tkconfigure, foreground="grey40")
		}
		if (type == "absunc") {
			tclvalue(X2) <- 0
			tclvalue(Y2) <- 0
			tclvalue(XY) <- 0
			sapply(list(X2.lab, Y2.lab, XY.lab), tkconfigure, foreground="grey40")
		}
		if (type == "absdiff") {
			tclvalue(X) <- 0
			tclvalue(Y) <- 0
			tclvalue(X2) <- 0
			tclvalue(Y2) <- 0
			tclvalue(XY) <- 0
			tclvalue(W) <- 0
			tclvalue(WY) <- -wx
			sapply(list(X.lab, Y.lab, X2.lab, Y2.lab, XY.lab, W.lab, WY.lab), tkconfigure, foreground="grey40")
		}

        tkrreplot(img, hscale=1.5, vscale=1.5)
    }

    replot <- function() {
		# read parameters from sliders
        type <- as.character(tclvalue(TYPE))
		b0 <- as.numeric(tclvalue(B0))
		x <- as.numeric(tclvalue(X))
		y <- as.numeric(tclvalue(Y))
		x2 <- as.numeric(tclvalue(X2))
		y2 <- as.numeric(tclvalue(Y2))
		xy <- as.numeric(tclvalue(XY))
		w <- as.numeric(tclvalue(W))
		wx <- as.numeric(tclvalue(WX))
		wy <- as.numeric(tclvalue(WY))
		rx <- as.numeric(tclvalue(RX))
		ry <- as.numeric(tclvalue(RY))
		rz <- as.numeric(tclvalue(RZ))
		
		plot(plotRSA(x=x, y=y, x2=x2, y2=y2, xy=xy, w=w, wx=wx, wy=wy, b0=b0, rotation=list(x=rx, y=ry, z=rz), zlim=zlim, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, zlab=zlab, points=points, demo=TRUE, type=type2, fit=fit, project=project, ...))
    }

	# define framework
    tt <- tktoplevel()
    tkwm.title(tt, "Response surface plot - polynomial model")

    img <- tkrplot(tt, replot, vscale=1.5, hscale=1.5)
    tkpack(img, side='left')
	
	# define radiobuttons
	tkpack(tfr <- tkframe(tt, relief='groove', borderwidth=3), side='top')
	
	tkpack(typebox <- tkframe(tfr), side='top', fill='x')
    tkpack(tklabel(typebox,text='Constraints: '), side='left',anchor='s')
	tkpack(tkradiobutton(typebox, variable=TYPE, command=update, value="all", text="All parameters"))
	tkpack(tkradiobutton(typebox, variable=TYPE, command=update, value="poly", text="Full polynomial"))
	tkpack(tkradiobutton(typebox, variable=TYPE, command=update, value="IA", text="Interaction"))
	tkpack(tkradiobutton(typebox, variable=TYPE, command=update, value="SQD", text="Squared difference"))
	tkpack(tkradiobutton(typebox, variable=TYPE, command=update, value="sq.shift", text="Shifted squared difference"))
	tkpack(tkradiobutton(typebox, variable=TYPE, command=update, value="sq.rot", text="Shifted and rotated squared difference"))
	tkpack(tkradiobutton(typebox, variable=TYPE, command=update, value="diff", text="Difference score X-Y"))
	
	tkpack(tkradiobutton(typebox, variable=TYPE, command=update, value="absunc", text="Unconstrained absolute difference"))
	tkpack(tkradiobutton(typebox, variable=TYPE, command=update, value="absdiff", text="Absolute difference"))

	

	# define sliders: polynomial model
	tkpack(tfr <- tkframe(tt, relief='groove', borderwidth=3), side='left')
	tkpack(fr0 <- tkframe(tfr), side='top',fill='x')
	tkpack(fr1 <- tkframe(tfr), side='top',fill='x')
	tkpack(fr2 <- tkframe(tfr), side='top',fill='x')
	tkpack(fr3 <- tkframe(tfr), side='top',fill='x')
	tkpack(fr4 <- tkframe(tfr), side='top',fill='x')
	tkpack(fr5 <- tkframe(tfr), side='top',fill='x')
	tkpack(fr6 <- tkframe(tfr), side='top',fill='x')
	tkpack(fr7 <- tkframe(tfr), side='top',fill='x')
	tkpack(fr8 <- tkframe(tfr), side='top',fill='x')
	B0.lab <- tklabel(fr0,text='Intercept: ')
	X.lab <- tklabel(fr1,text='x: ')
	Y.lab <- tklabel(fr2,text='y: ')
	XY.lab <- tklabel(fr3,text='xy: ')
	X2.lab <- tklabel(fr4,text='x2: ')
	Y2.lab <- tklabel(fr5,text='y2: ')
	W.lab <- tklabel(fr6,text='w: ')
	WX.lab <- tklabel(fr7,text='wx: ')
	WY.lab <- tklabel(fr8,text='wy: ')
	
    tkpack(B0.lab, side='left',anchor='s')
	tkpack(tkscale(fr0, variable=B0, orient='horizontal', command=update, from=-5, to=5, resolution=.1), side='left')
	
    tkpack(X.lab, side='left',anchor='s')
	tkpack(tkscale(fr1, variable=X, orient='horizontal', command=update, from=ifelse(is.null(fit), -5, -abs(x.0)*2), to=ifelse(is.null(fit), 5, abs(x.0)*2), resolution=0.01), side='left')

    tkpack(Y.lab, side='left',anchor='s')
	tkpack(tkscale(fr2, variable=Y, orient='horizontal', command=update, from=ifelse(is.null(fit), -5, -abs(y.0)*2), to=ifelse(is.null(fit), 5, abs(y.0)*2), resolution=0.01), side='left')
	
    tkpack(XY.lab, side='left',anchor='s')
	tkpack(tkscale(fr3, variable=XY, orient='horizontal', command=update, from=ifelse(is.null(fit), -3, -abs(xy.0)*2), to=ifelse(is.null(fit), 3, abs(xy.0)*2), resolution=0.01), side='left')
	
    tkpack(X2.lab, side='left',anchor='s')
	tkpack(tkscale(fr4, variable=X2, orient='horizontal', command=update, from=ifelse(is.null(fit), -3, -abs(x2.0)*2), to=ifelse(is.null(fit), 3, abs(x2.0)*2), resolution=0.01), side='left')

    tkpack(Y2.lab, side='left',anchor='s')
	tkpack(tkscale(fr5, variable=Y2, orient='horizontal', command=update, from=ifelse(is.null(fit), -3, -abs(y2.0)*2), to=ifelse(is.null(fit), 3, abs(y2.0)*2), resolution=0.01), side='left')
	
	# define sliders: absdiff model
	tkpack(tfr <- tkframe(tt, relief='groove', borderwidth=3), side='right')
	
    tkpack(W.lab, side='left',anchor='s')
	tkpack(tkscale(fr6, variable=W, orient='horizontal', command=update, from=ifelse(is.null(fit), -5, -abs(w.0)*2), to=ifelse(is.null(fit), 5, abs(w.0)*2), resolution=0.01), side='left')

    tkpack(WX.lab, side='left',anchor='s')
	tkpack(tkscale(fr7, variable=WX, orient='horizontal', command=update, from=ifelse(is.null(fit), -1, -abs(wx.0)*2), to=ifelse(is.null(fit), 1, abs(wx.0)*2), resolution=0.01), side='left')
	
    tkpack(WY.lab, side='left',anchor='s')
	tkpack(tkscale(fr8, variable=WY, orient='horizontal', command=update, from=ifelse(is.null(fit), -1, -abs(wy.0)*2), to=ifelse(is.null(fit), 1, abs(wy.0)*2), resolution=0.01), side='left')
	
	
	## Rotation of display
	tkpack(tfr3d <- tkframe(tt, relief='groove', borderwidth=3), side='right')
	tkpack(fr3.1 <- tkframe(tfr3d), side='top',fill='x')
	tkpack(fr3.2 <- tkframe(tfr3d), side='top',fill='x')
	tkpack(fr3.3 <- tkframe(tfr3d), side='top',fill='x')
	X3.lab <- tklabel(fr3.1,text='x rotation: ')
	Y3.lab <- tklabel(fr3.2,text='y rotation: ')
	Z3.lab <- tklabel(fr3.3,text='z rotation: ')
    tkpack(X3.lab, side='left',anchor='s')
	tkpack(tkscale(fr3.1, variable=RX, orient='horizontal', command=update, from=-90, to=90, resolution=1), side='left')

    tkpack(Y3.lab, side='left',anchor='s')
	tkpack(tkscale(fr3.2, variable=RY, orient='horizontal', command=update, from=-90, to=90, resolution=1), side='left')
	
    tkpack(Z3.lab, side='left',anchor='s')
	tkpack(tkscale(fr3.3, variable=RZ, orient='horizontal', command=update, from=-90, to=90, resolution=1), side='left')
	
	
	
	## Extra (dummy) parameters
	tkpack(frROT1 <- tkframe(tfr3d), side='top',fill='x')
	tkpack(frROT2 <- tkframe(tfr3d), side='top',fill='x')
	
    tkpack(tklabel(frROT1,text='Shift (C): '), side='left',anchor='s')
	tkpack(tkscale(frROT1, variable=C, orient='horizontal', command=update, from=-20, to=20, resolution=0.1), side='left')

    tkpack(tklabel(frROT2,text='Rotation (S): '), side='left',anchor='s')
	tkpack(tkscale(frROT2, variable=S, orient='horizontal', command=update, from=0, to=3, resolution=0.1), side='left')
	
    return(invisible(NULL))
}


#demoRSA()
#demoRSA(fit=r1, points=TRUE)

# Hack to please CRAN:
if(getRversion() >= "2.15.1")  {
	utils::globalVariables(c('tclVar', 'tclvalue', 'tkconfigure' , 'tkframe', 'tklabel', 'tkpack', 'tkradiobutton', 'tkscale', 'tktoplevel', 'tkwm.title'))
}