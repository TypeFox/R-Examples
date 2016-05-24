#' @title Plots a response surface of a polynomial equation of second degree
#'
#' @description
#' Plots an RSA object, or a response surface with specified parameters
#'
#' @details
#' Each plot type has its distinctive advantages. The two-dimensional contour plot gives a clear view of the position of the principal axes and the stationary point. The 3d plot gives a three dimensional impression of the surface, allows overplotting of the original data points (in case an RSA object is provided), and allows the interactive adjustment of regression weights in the \code{\link{RSA}} function. The interactive plot allows rotating and exploring a three-dimensional surface with the mouse (nice for demonstration purposes).
#' If you want to export publication-ready plots, it is recommended to export it with following commands:
#' \code{p1 <- plot(r1, bw=TRUE)
#' trellis.device(device="cairo_pdf", filename="RSA_plot.pdf")
#' print(p1)
#' dev.off()}
#'
#' @aliases plotRSA
#'
#' @importFrom aplpack compute.bagplot
#'
#' @export
#' @param x Either an RSA object (returned by the \code{RSA} function), or the coefficient for the X predictor
#' @param y Y coefficient
#' @param x2 X^2 coefficient
#' @param y2 Y^2 coefficient
#' @param xy XY interaction coefficient
#' @param w W coefficient (for (un)constrained absolute difference model)
#' @param wx WX coefficient (for (un)constrained absolute difference model)
#' @param wy WY coefficient (for (un)constrained absolute difference model)
#' @param y3 Y^3 coefficient
#' @param x3 X^3 coefficient
#' @param xy2 XY^2 coefficient
#' @param x2y X^2Y coefficient
#' @param b0 Intercept
#' @param xlim Limits of the x axis
#' @param ylim Limits of the y axis
#' @param zlim Limits of the z axis
#' @param xlab Label for x axis
#' @param ylab Label for y axis
#' @param zlab Label for z axis
#' @param main the main title of the plot
#' @param cex.main Factor for main title size
#' @param surface Method for the calculation of the surface z values. "predict" takes the predicted values from the model, "smooth" uses a thin plate smoother (function \code{Tps} from the \code{fields} package) of the raw data
#' @param lambda lambda parameter for the smoother. Default (NULL) means that it is estimated by the smoother function. Small lambdas around 1 lead to rugged surfaces, big lambdas to very smooth surfaces.
#' @param rotation Rotation of the 3d surface plot (when type == "3d")
#' @param label.rotation Rotation of the axis labls (when type == "3d")
#' @param gridsize Number of grid nodes in each dimension
#' @param bw Print surface in black and white instead of colors?
#' @param legend Print color legend for z values?
#' @param cex Font size factor for axes labels and axes titles
#' @param type \code{3d} for 3d surface plot, \code{contour} for 2d contour plot, "interactive" for interactive rotatable plot. Shortcuts (i.e., first letter of string) are sufficient
#' @param points A list of parameters which define the appearance of the raw scatter points: 
#'	\itemize{
#'		\item data: Data frame which contains the coordinates of the raw data points. First column = x, second = y, third = z. This data frame is automatically generated when the plot is based on a fitted RSA-object
#'		\item show = TRUE: Should the original data points be overplotted?
#'		\item color = "black": Color of the points
#' 		\item value="raw": Plot the original z value, "predicted": plot the predicted z value
#'		\item jitter = 0: Amount of jitter for the raw data points. For z values, a value of 0.005 is reasonable
#'		\item cex = .5: multiplication factor for point size
#' 		\item out.mark = FALSE: If set to TRUE, outliers according to Bollen & Jackman (1980) are printed as red X symbols, but only when they have been removed in the RSA function: \code{RSA(..., out.rm=TRUE)}.
#'			\itemize{
#'				\item If out.rm == TRUE (in RSA()) and out.mark == FALSE (in plotRSA()), the outlier is removed from the model and *not plotted* in plotRSA.
#'				\item If out.rm == TRUE (in RSA()) and out.mark == TRUE (in plotRSA()), the outlier is removed from the model but plotted and marked in plotRSA.
#'				\item If out.rm == FALSE (in RSA()): Outliers are not removed and cannot be plotted.
#'				\item Example syntax: \code{plotRSA(r1, points=list(show=TRUE, out.mark=TRUE))}
#'		}
#'	}
#' As a shortcut, you can also set \code{points=TRUE} to set the defaults.

#' @param model If x is an RSA object: from which model should the response surface be computed?
#' @param demo Do not change that parameter (internal use only)
#' @param fit Do not change that parameter (internal use only)
#' @param param Should the surface parameters a1 to a4 be shown on the plot? In case of a 3d plot a1 to a4 are printed on top of the plot; in case of a contour plot the principal axes are plotted.
#' @param coefs Should the regression coefficients b1 to b5 be shown on the plot? (Only for 3d plot)
#' @param axes A vector of strings specifying the axes that should be plotted. Can be any combination of c("LOC", "LOIC", "PA1", "PA2"). LOC = line of congruence, LOIC = line of incongruence, PA1 = first principal axis, PA2 = second principal axis
#' @param project A vector of graphic elements that should be projected on the floor of the cube. Can include any combination of c("LOC", "LOIC", "PA1", "PA2", "contour", "points")
#' @param maxlines Should the maximum lines be plotted? (red: maximum X for a given Y, blue: maximum Y for a given X). Works only in type="3d"
#' @param link Link function to transform the z axes. Implemented are "identity" (no transformation; default), "probit", and "logit"
#' @param suppress.surface Should the surface be suppressed (only for \code{type="3d"})? Useful for only showing the data points, or for didactic purposes (e.g., first show the cube, then fade in the surface).
#' @param suppress.box Should the surrounding box be suppressed (only for \code{type="3d"})?
#' @param border Should a thicker border around the surface be plotted? Sometimes this border leaves the surrounding box, which does not look good. In this case the border can be suppressed by setting \code{border=FALSE}.
#' @param contour A list defining the appearance of contour lines (aka. height lines). show=TRUE: Should the contour lines be plotted on the 3d wireframe plot? (Parameter only relevant for \code{type="3d"}). color = "grey40": Color of the contour lines. highlight = c(): A vector of heights which should be highlighted (i.e., printed in bold). Be careful: the highlighted line is not necessarily exactly at the specified height; instead the nearest height line is selected.
#' @param hull Plot a bag plot on the surface (This is a bivariate extension of the boxplot. 50\% of points are in the inner bag, 50\% in the outer region). See Rousseeuw, Ruts, & Tukey (1999).
#' @param showSP Plot the stationary point? (only relevant for \code{type="contour"})
#' @param showSP.CI Plot the CI of the stationary point? (only relevant for \code{type="contour"})
#' @param distance A vector of three values defining the distance of labels to the axes
#' @param tck A vector of three values defining the position of labels to the axes (see ?wireframe)
#' @param pal A palette for shading. You can use \code{\link{colorRampPalette}} to construct a color ramp, e.g. \code{plot(r.m, pal=colorRampPalette(c("darkgreen", "yellow", "darkred"))(20))}. If \code{pal="flip"}, the default palette is used, but reversed (so that red is on top and green on the bottom).
#' @param pal.range Should the color range be scaled to the box (\code{pal.range = "box"}, default), or to the min and max of the surface (\code{pal.range = "surface"})? If set to "box", different surface plots can be compared along their color, as long as the zlim is the same for both.
#' @param pad Pad controls the margin around the figure (positive numbers: larger margin, negative numbers: smaller margin)
#' @param ... Additional parameters passed to the plotting function (e.g., sub="Title"). A useful title might be the R squared of the plotted model: \code{sub = as.expression(bquote(R^2==.(round(getPar(x, "r2", model="full"), 3))))}
#'
#' @references
#' Rousseeuw, P. J., Ruts, I., & Tukey, J. W. (1999). The Bagplot: A Bivariate Boxplot. The American Statistician, 53(4), 382-387. doi:10.1080/00031305.1999.10474494
#' @seealso \code{\link{demoRSA}}, \code{\link{RSA}}
#'
#' @examples
#' # Plot response surfaces from known parameters
#' # example of Edwards (2002), Figure 3
#' # Default: 3d plot:
#' plotRSA(x=.314, y=-.118, x2=-.145, y2=-.102, xy=.299, b0=5.628)
#' # Contour plot:
#' plotRSA(x=.314, y=-.118, x2=-.145, y2=-.102, xy=.299, b0=5.628, type="c")
#' \dontrun{
#' # Interactive plot (try the mouse!):
#' plotRSA(x=.314, y=-.118, x2=-.145, y2=-.102, xy=.299, b0=5.628, type="i")
#' }
#'
#' # Plot response surface from an RSA object
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
#' r1 <- RSA(z.sq~x*y, df, models=c("SQD", "full", "IA"))
#' plot(r1)	# default: model = "full"
#' plot(r1, model="SQD", points=list(show=TRUE, value="predicted"))



#b0=0; x=0; y=0; x2=0; y2=0; xy=0; w=0; wx=0; wy=0;  zlim=NULL; xlim=c(-2, 2); ylim=c(-2, 2); rotation=list(x=-45, y=45, z=35); legend=TRUE; cex=1.2; type="3d"; points=TRUE; demo=FALSE; model="full"; 

#b0=-9; x=0; y=0; x2=0; y2=0; xy=0; w=0; wx=1; wy=-1;  zlim=NULL; xlim=c(-2, 2); ylim=c(-2, 2); rotation=list(x=-45, y=45, z=35); legend=TRUE; cex=1.2; type="3d"; points=TRUE; demo=FALSE; model="full"; fit=NULL; link="identity"; param=TRUE; gridsize=21;bw=FALSE; pal=NULL; axes=c("LOC", "LOIC", "PA1", "PA2"); distance=c(1, 1, 1); tck=c(1, 1, 1); xlab="X"; ylab="Y"; zlab="Z"; border=TRUE;

## old rotation
# rotation=list(x=-45, y=45, z=35), label.rotation=list(x=45, y=-25, z=94)
# distance=c(1, 1, 1), tck=c(1, 1, 1)

plotRSA <- function(x=0, y=0, x2=0, y2=0, xy=0, w=0, wx=0, wy=0, x3=0, xy2=0, x2y=0, y3=0, b0=0, 
	type="3d", model="full", 
	xlim=NULL, ylim=NULL, zlim=NULL, 
	xlab=NULL, ylab=NULL, zlab=NULL, main="",
	surface="predict", lambda=NULL, 
	suppress.surface=FALSE, suppress.box = FALSE,
	rotation=list(x=-63, y=32, z=15), label.rotation=list(x=19, y=-40, z=92), 
	gridsize=21, bw=FALSE, legend=TRUE, param=TRUE, coefs=FALSE,
	axes=c("LOC", "LOIC", "PA1", "PA2"), 
	project=c("contour"), maxlines=FALSE,
	cex=1, cex.main=1, 
	points = list(data=NULL, show=NA, value="raw", jitter=0, color="black", cex=.5, out.mark=FALSE),
	fit=NULL, link="identity", 
	tck=c(1.5, 1.5, 1.5), distance=c(1.3, 1.3, 1.4), border=FALSE, 
	contour = list(show=FALSE, color="grey40", highlight = c()),
	hull=NA, showSP=FALSE, showSP.CI=FALSE, 
	pal=NULL, pal.range="box", 
	pad=0, demo=FALSE, ...) {
	
	
	# ---------------------------------------------------------------------
	# Warnings and error handling ...
	if (!identical(xlim, ylim)) {print("Note: Axes dimensions are not equal. The visual diagonal is *not* the line of numerical congruence! Consider choosing identical values for xlim and ylim.")}
		
	if (class(x) == "RSA") {
		stop("If you want to plot an RSA object, please use plot(...); plotRSA should be only used when you directly provide the regression coefficients.")
	}
	
	# remove LOC, LOIC etc. when they do not make sense.
	if (any(c(x3, xy2, x2y, y3, w, wx, wy) != 0)) {
		axes <- ""
		project <- project[!project %in% c("PA1", "PA2", "LOC", "LOIC")]
	}
	
	
	# define the defaults
	if (is.null(points) || (typeof(points) == "logical" && points == TRUE)) {
		points <- list(show=TRUE, value="raw", jitter=0, color="black", cex=.5, out.mark=FALSE)
	}
	if (is.null(points) || (typeof(points) == "logical" && points == FALSE)) {
		points <- list(show=FALSE)
	}
	
	if (is.null(points$show)) points$show <- TRUE
	if (is.na(points$show)) {
		if (is.null(points$data)) {
			points$show <- FALSE	
		} else {
			points$show <- TRUE
		}
	}
	if (is.null(points$value)) points$value <- "raw"
	if (is.null(points$color)) points$color <- "black"
	if (is.null(points$jitter)) points$jitter <- 0
	if (is.null(points$cex)) points$cex <- 0.5
	if (is.null(points$out.mark)) points$out.mark <- FALSE
	if (points$show==TRUE & is.null(points$data)) {
		warning("You must provide a data frame with the coordinates of the raw data points (points = list(show = TRUE, data = ???)). Points are not plotted.")
		points$show <- FALSE
	}
	if (!is.null(fit)) {
		if (points$out.mark==TRUE & fit$out.rm==FALSE) {
			warning("Outliers can only be marked in the plot when they were removed in the RSA function: RSA(..., out.rm=TRUE).")
			points$out.mark <- FALSE
		}
	}
		
	if (is.null(contour$show)) contour$show <- TRUE
	if (is.null(contour$color)) contour$color <- "grey40"
	if (is.null(contour$highlight)) contour$highlight <- c()
	
	# define default behavior for the "hull" parameter
	if (is.na(hull)) {
		if (!is.null(fit) | !is.null(points$data)) {
			hull <- TRUE
		} else {
			hull <- FALSE
		}
	}
	
	type <- match.arg(type, c("interactive", "3d", "contour"))
	surface <- match.arg(surface, c("predict", "smooth"))
	points[["value"]] <- match.arg(points[["value"]], c("raw", "predicted"))

	if (demo == FALSE) {
			if (is.null(xlab)) {
				if (!is.null(points$data)) {
					xlab <- colnames(points$data)[1]
				} else {
					xlab <- "X"
				}
			}
			if (is.null(ylab)) {
				if (!is.null(points$data)) {
					ylab <- colnames(points$data)[2]
				} else {
					ylab <- "Y"
				}
			}
			if (is.null(zlab)) {
				if (!is.null(points$data)) {
					zlab <- colnames(points$data)[3]
				} else {
					zlab <- "Z"
				}
			}
			
			if (is.null(xlim)) {xlim <- c(-2.1, 2.1)}
			if (is.null(ylim)) {ylim <- c(-2.1, 2.1)}
	}
	
	
	if (is.null(points$data) & surface == "smooth") {
		warning("Smoothing only works if data points are provided (points=list(data=???))! Reverting to surface = 'predict'")
		surface <- "predict"
	}
	
	C <- c(x, y, x2, y2, xy, w, wx, wy, x3, xy2, x2y, y3)
	
	if (!model %in% c("absunc", "absdiff")) {
		if (!is.null(fit) & model != "cubic" & model != "null") {
			SP <- RSA.ST(fit, model=model)
			PAR <- getPar(fit, "coef", model=model)
			
			SP.text <- paste0("a", 1:4, ": ", f2(SP$SP$estimate, 2), p2star(SP$SP$p.value), collapse="   ")
			
			a4rs_par <- PAR[PAR$label == "a4.rescaled", ]
			if (nrow(a4rs_par) == 1) {
				SP.text <- paste0(SP.text, "/ a4(rescaled): ", f2(a4rs_par$est, 2), p2star(a4rs_par$pvalue))
			}
			SP.text <- paste0(SP.text, "\n")
			
			
			meanlevel <- PAR[PAR$label == "meaneffect", ]
			if (nrow(meanlevel) == 1) {
				SP.text <- paste0(SP.text, "mean-level effect = ", f2(meanlevel$est, 2), p2star(meanlevel$pvalue), "    ")
			}
			
			C_par <- PAR[PAR$label == "C", ]
			if (nrow(C_par) == 1) {
				SP.text <- paste0(SP.text, "C = ", f2(C_par$est, 2), p2star(C_par$pvalue), "    ")
			}
			
			S_par <- PAR[PAR$label == "S", ]
			if (nrow(S_par) == 1) {
				SP.text <- paste0(SP.text, "S = ", f2(S_par$est, 2), p2star(S_par$pvalue), "    ")
			}
			
		} else {
			SP <- RSA.ST(x=x, y=y, xy=xy, x2=x2, y2=y2)
			SP.text <- paste0("a", 1:4, ": ", f2(SP$SP$estimate, 2), p2star(SP$SP$p.value), collapse="    ")			
		}		
	} else {
		SP <- NULL
		param <- FALSE
		SP.text <- ""
	}
	
	# Print coefs in 3d plot?
	COEFS <- ""
	if (coefs == TRUE) {
		COEFS <- paste0("b1 = ", f2(x, 3), "\n", "b2 = ", f2(y, 3), "\n", "b3 = ", f2(x2, 3), "\n", "b4 = ", f2(xy, 3), "\n", "b5 = ", f2(y2, 3), "\n")
	}
	
	
	# ---------------------------------------------------------------------
	# Calculate positions of raw points

	xpoints <- ypoints <- zpoints <- NA
	if (points$out.mark == TRUE & is.null(fit)) {
		warning("Outliers can only be marked if an RSA-object is provided. Points are not plotted.")
		points$show <- FALSE
	}
	if (points$show == TRUE | hull==TRUE) {
		if (is.null(points$data)) stop("You must provide raw data if you want to plot raw points or the bagplot.")
		if (points$out.mark == FALSE) {
			data.used <- points$data
		}
		if (points$out.mark == TRUE & !is.null(fit)) {
			data.used <- fit$data[, c(fit$IV1, fit$IV2, fit$DV)]	# this includes all data points
		}
	
		if (points$jitter > 0) {
			data.used[, 1] <- data.used[, 1] + rnorm(length(data.used[, 1]), 0, points$jitter)
			data.used[, 2] <- data.used[, 2] + rnorm(length(data.used[, 2]), 0, points$jitter)
		}
	
		if (points$value == "raw") {
			zpoints <- as.vector(data.used[, 3])
		} else if (points$value == "predicted") {
			N <- colnames(data.used)
			data.used2 <- add.variables(formula(paste0(N[3], " ~ ", N[1], "*", N[2])), data.used)
			
			# calculate predicted values
			zpoints <- b0 + colSums(C*t(data.used2[, c(
				N[1],
				N[2],
				paste0(N[1], "2"),
				paste0(N[2], "2"),
				paste0(N[1], "_", N[2]),
				"W",
				paste0("W_", N[1]),
				paste0("W_", N[2]),
				paste0(N[1], "3"),
				paste0(N[1], "_", N[2], "2"),
				paste0(N[1], "2", "_", N[2]),
				paste0(N[2], "3")
			)]))
			
			zpoints <- as.vector(zpoints)
		}

		xpoints <- as.vector(data.used[, 1])
		ypoints <- as.vector(data.used[, 2])	
	}
	
	
		
	# build data set
	grid <- gridsize
	new <- data.frame(x = rep(seq(xlim[1], xlim[2], length.out=grid), grid), y = rep(seq(ylim[1], ylim[2], length.out=grid), each=grid))
	new2 <- add.variables(z~x+y, new)
		
	# calculate z values
	if (surface == "predict") {
		new2$z <- b0 + colSums(C*t(new2[, c(1:5, 9:11, 15:18)]))
	}
	if (surface == "smooth") {
		
		if (!requireNamespace("fields", quietly = TRUE)) {
			stop('`fields` package needed for smooth surfaces to work. Please install it with install.packages("fields")', call. = FALSE)
		}
		
		tpsfit <- fields::Tps(points$data[, 1:2], points$data[, 3], scale.type="unscaled", lambda=lambda)
		new2$z <- fields::predict.Krig(tpsfit, new[, c("x", "y")])
		param <- FALSE
		axes <- ""
	}
	
	# impose link functions
	logit <- function (x) {log(x/(1-x))}
	invlogit <- function (x) {1/(1+exp(-x))}
	link <- match.arg(link, c("identity", "logit", "probit"))
	if (link == "probit") {
		z.trans <- 1.7 * new2$z
		new2$z <- invlogit(z.trans)
	}
	if (link == "logit") {
		new2$z <- invlogit(new2$z)
	}


	# determine zlim
	if (!is.null(points$data) & demo==FALSE & is.null(zlim)) {
		# old: set zlim according to fitted surface
		#zlim <- c(min(min(new2$z, na.rm=TRUE), min(points$data[, 3], na.rm=TRUE)), max(max(new2$z, na.rm=TRUE), max(points$data[, 3], na.rm=TRUE)))
		
		# new: set zlim according to actual data range
		zlim <- c(min(points$data[, 3], na.rm=TRUE), max(points$data[, 3], na.rm=TRUE))
	} else {
		if (is.null(zlim)) zlim <- c(min(new2$z), max(new2$z))
	}
	zlim.final <- zlim
	
	
	# Catch border case: completely flat surface. Remove contour, redefine zlim
	if (var(new2$z) == 0) {
		contour$show <- FALSE
		project <- project[-which(project == "contour")]
		if (zlim[1] == 0 && zlim[2] == 0) zlim <- xlim
	}
	
	
	## Define colors
	
	# flip palette?
	flip <- FALSE
	if (!is.null(pal) && pal=="flip") {
		flip <- TRUE
		pal <- NULL
	}
	if (bw == FALSE) {
		
		# RdYlGn palette
		if (is.null(pal)) {
			pal <- c("#A50026","#D73027","#F46D43","#FDAE61","#FEE08B","#FFFFBF","#D9EF8B","#A6D96A","#66BD63","#1A9850","#006837")
			if (flip==TRUE) {pal <- rev(pal)}
		}
		
		gridCol <- ifelse(contour$show == TRUE, "grey60", "grey30")
		
		LOC.col <- LOIC.col <- "blue"
		PA1.col <- PA2.col <- "grey30"
	} else {
		if (is.null(pal)) {
			pal <- colorRampPalette(c("#FFFFFF", "#AAAAAA", "#030303"), bias=2)(11)
			if (flip==TRUE) {pal <- rev(pal)}
		}
		gridCol <- ifelse(contour$show == TRUE, "grey10", "grey10")
		
		LOC.col <- LOIC.col <- "black"
		PA1.col <- PA2.col <- "black"
	}
	if (length(pal) < 2) {legend <- FALSE}
	
	
	
	# ---------------------------------------------------------------------
	#  calculate bag plot: bag = outer, loop = inner
	
	if (hull==TRUE) {
		BAG <- aplpack::compute.bagplot(xpoints, ypoints)
		
		# close the polygon
		h1 <- rbind(BAG$hull.bag, BAG$hull.bag[1, ])
		h2 <- rbind(BAG$hull.loop, BAG$hull.loop[1, ])

		# approx: interpolate the points of the bag (in order to get a more smooth fitting line on the z-axis)
		minDist <- min(diff(xlim)/gridsize, diff(ylim)/gridsize)/2

		h1 <- interpolatePolygon(h1[, 1], h1[, 2], minDist=minDist)	
		h2 <- interpolatePolygon(h2[, 1], h2[, 2], minDist=minDist)	
		
		# calculate predicted values
		bagpoints <- add.variables(z~x+y, data.frame(x=h1$x, y=h1$y))
		bagpoints$z <- b0 + colSums(C*t(bagpoints[, c(1:5, 9:11, 15:18)]))
		bag <- data.frame(X  = bagpoints$x, Y  = bagpoints$y, Z = bagpoints$z)
		
		looppoints <- add.variables(z~x+y, data.frame(x=h2$x, y=h2$y))
		looppoints$z <- b0 + colSums(C*t(looppoints[, c(1:5, 9:11, 15:18)]))
		loop <- data.frame(X  = looppoints$x, Y  = looppoints$y, Z = looppoints$z)			
	}
	
		
	## ======================================================================
	## Interactive plot
	## ======================================================================
	
	if (type == "interactive") {
		if (!requireNamespace("rgl", quietly = TRUE)) {
			stop("`rgl` package needed for interactive plots. Please install it with install.packages('rgl').", call. = FALSE)
		}
		
		P <- list(x=seq(xlim[1], xlim[2], length.out=grid), y=seq(ylim[1], ylim[2], length.out=grid))
		DV2 <- matrix(new2$z, nrow=grid, ncol=grid, byrow=FALSE)
		R <- range(DV2)
		col2 <- as.character(cut(1:(R[2] - R[1] + 1), breaks=length(pal), labels=pal))
		
		rgl::open3d(cex=cex)
		rgl::rgl.viewpoint(-30, -90, fov=0)
		rgl::rgl.light(theta = 0, phi = 90, viewpoint.rel = TRUE, ambient = "#FF0000", diffuse = "#FFFFFF", specular = "#FFFFFF")
		rgl::persp3d(P$x, P$y, DV2, xlab = xlab, ylab = ylab, zlab = zlab, color=col2[DV2 - R[1] + 1], main=main, ...)

		if (contour$show == TRUE) {
		    contours <- contourLines(P, z=DV2)
		     for (i in 1:length(contours)) {
				 with(contours[[i]], rgl::lines3d(x, y, level, col=contour$color))
			 }
		 }
		
		if (points$show == TRUE) {
			if (points$out.mark == FALSE) {
				rgl::points3d(data.frame(xpoints, ypoints, zpoints), col=points$color)
			}
			if (points$out.mark == TRUE) {
				if (!is.null(fit)) {
					colvec <- rep(points$color, nrow(fit$data))
					colvec[fit$outliers] <- "red"
					rgl::points3d(fit$data[, c(fit$IV1, fit$IV2, fit$DV)], col=colvec)
					rgl::text3d(fit$data[fit$outliers, c(fit$IV1, fit$IV2, fit$DV)], col="red", texts="X")
				} else {
					warning("Please provide an RSA-object to mark outliers.")
				}
			}
		}
		
		p1 <- NULL	# no plot object is returned	
	}	
	
	
	
	## ======================================================================
	## Wireframe plot
	## ======================================================================
	
		
	if (type == "3d") {
				
			mypanel2 <- function(x, y, z, xlim, ylim, zlim, xlim.scaled, ylim.scaled, zlim.scaled, axes, axesList, x.points=NULL, y.points=NULL, z.points=NULL, SPs="", ...) {
				
				
			   # rescale absolute x, y, and z values so that they fit into the box
			   RESCALE.Z <- function(z1) {
	              Z2 <- zlim.scaled[1] + diff(zlim.scaled) * (z1 - zlim[1]) / diff(zlim)
				  return(Z2)
			   }
			   RESCALE <- function(n) {
				  X2 <- xlim.scaled[1] + diff(xlim.scaled) * (n$X - xlim[1]) / diff(xlim)
	              Y2 <- ylim.scaled[1] + diff(ylim.scaled) * (n$Y - ylim[1]) / diff(ylim)
	              Z2 <- zlim.scaled[1] + diff(zlim.scaled) * (n$Z - zlim[1]) / diff(zlim)
				  df <- data.frame(X=X2, Y=Y2, Z=Z2)
				  df <- df[df$X >= min(xlim.scaled) & df$X <= max(xlim.scaled) & df$Y >= min(ylim.scaled) & df$Y <= max(ylim.scaled) &  df$Z >= min(zlim.scaled) & df$Z <= max(zlim.scaled), ]
				  return(df)
			   }
				
			
			# ---------------------------------------------------------------------
			# 1. Projection on bottom of cube
			  if (length(project) > 0) {
				  for (p in project) {
					  if (p %in% c("LOC", "LOIC", "PA1", "PA2")) {
						  if (is.null(axesList[[p]])) break;
						  a0 <- RESCALE(getIntersect2(p0=axesList[[p]]$p0, p1=axesList[[p]]$p1))
						  if (nrow(a0) <= 1) break;
							  panel.3dscatter(x = a0$X, y = a0$Y, z = rep(RESCALE.Z(min(zlim.final) + .01), nrow(a0)), 
						  				xlim = xlim, ylim = ylim, zlim = zlim,
			                            xlim.scaled = xlim.scaled, ylim.scaled = ylim.scaled, zlim.scaled = zlim.scaled, 
										type="l", col.line=axesList[[p]]$col, lty=axesList[[p]]$lty, lwd=2, ...)
					  }
					  
					  if (p == "hull") {
	  					  bag.rescale <- RESCALE(bag)
						  panel.3dscatter(x = bag.rescale$X, y = bag.rescale$Y, z = rep(RESCALE.Z(min(zlim.final) + .01), nrow(bag.rescale)), 
					  				xlim = xlim, ylim = ylim, zlim = zlim,
		                            xlim.scaled = xlim.scaled, ylim.scaled = ylim.scaled, zlim.scaled = zlim.scaled, 
									type="l", col.line="grey30", lty="dashed", lwd=2, ...)
									
  	  					  loop.rescale <- RESCALE(loop)
  						  panel.3dscatter(x = loop.rescale$X, y = loop.rescale$Y, z = rep(RESCALE.Z(min(zlim.final) + .01), nrow(loop.rescale)), 
  					  				xlim = xlim, ylim = ylim, zlim = zlim,
  		                            xlim.scaled = xlim.scaled, ylim.scaled = ylim.scaled, zlim.scaled = zlim.scaled, 
  									type="l", col.line="black", lty="dashed", lwd=2, ...)			
										
					  }
					  
					  if (p == "points") {
  			              x2 <- xlim.scaled[1] + diff(xlim.scaled) * (x.points - xlim[1]) / diff(xlim)
  			              y2 <- ylim.scaled[1] + diff(ylim.scaled) * (y.points - ylim[1]) / diff(ylim)
						  z2 <- rep(RESCALE.Z(min(zlim.final) + .01), length(x2))
						
  			              panel.3dscatter(x = x2, y = y2, z = z2, 
							  			xlim = xlim, ylim = ylim, zlim = zlim,
  			                              xlim.scaled = xlim.scaled, ylim.scaled = ylim.scaled, zlim.scaled = zlim.scaled,
  										  pch=20, col=points$color, cex=points$cex, ...)
					  }
				  }
			  }
				  
				  
				  # project contour on bottom
				if (contour$show == TRUE | "contour" %in% project) {
					# "abuse" ggplot to compute the contour lines
					cs <- ggplot(new2, aes_string(x="x", y="y", fill="z", z="z")) + stat_contour(bins=ifelse(length(pal)>1, length(pal)+1, 8))
					cLines <- ggplot_build(cs)
					C0 <- cLines$data[[1]][, c("x", "y", "level", "group")]
					colnames(C0) <- c("X", "Y", "Z", "group")	# C0 keeps the contour lines
				
				
					if ("contour" %in% project) {	
						for (cL in C0$group) {
						  C1 <- RESCALE(C0[C0$group==cL, c("X", "Y", "Z")])
						  panel.3dscatter(x = C1$X, y = C1$Y, z = rep(RESCALE.Z(min(zlim.final) + .01), nrow(C1)), 
							 		xlim = xlim, ylim = ylim, zlim = zlim,
				                    xlim.scaled = xlim.scaled, ylim.scaled = ylim.scaled, zlim.scaled = zlim.scaled,
									type="l", col.line=contour$color, lty="solid", lwd=1, ...)
						}	
					}
				}
				  
		  	
				
			   # ---------------------------------------------------------------------
			   # 2. Borders, back part
			   
			   
					if (border==TRUE & suppress.surface==FALSE) {
						  # Make boundary of grid a bit thicker
						  box1 <- new2[new2$y == max(new2$y), ]
						  box2 <- new2[new2$y == min(new2$y), ]
						  box3 <- new2[new2$x == max(new2$x), ]
						  box4 <- new2[new2$x == min(new2$x), ]
						  box <- rbind(data.frame(box1, side=1), data.frame(box2, side=2), data.frame(box3, side=3), data.frame(box4, side=4))
						  
			              x.box <- xlim.scaled[1] + diff(xlim.scaled) * (box$x - xlim[1]) / diff(xlim)
			              y.box <- ylim.scaled[1] + diff(ylim.scaled) * (box$y - ylim[1]) / diff(ylim)
			              z.box <- zlim.scaled[1] + diff(zlim.scaled) * (box$z - zlim[1]) / diff(zlim)
						  
						  # plot the back lines of the border
			              panel.3dscatter(x = x.box[box$side==1], y = y.box[box$side==1], z = z.box[box$side==1], xlim = xlim, ylim = ylim, zlim = zlim, xlim.scaled = xlim.scaled, ylim.scaled = ylim.scaled, zlim.scaled = zlim.scaled, type="l", col.line=gridCol, lwd=4, ...)
						  panel.3dscatter(x = x.box[box$side==3], y = y.box[box$side==3], z = z.box[box$side==3], xlim = xlim, ylim = ylim, zlim = zlim, xlim.scaled = xlim.scaled, ylim.scaled = ylim.scaled, zlim.scaled = zlim.scaled, type="l", col.line=gridCol, lwd=4, ...)
					  }
				
					  # ---------------------------------------------------------------------
					  # 3. the surface
					  if (suppress.surface==FALSE) {
						  panel.3dwire(x = x, y = y, z = z, xlim = xlim, ylim = ylim, zlim = zlim,
		                           xlim.scaled = xlim.scaled, ylim.scaled = ylim.scaled, zlim.scaled = zlim.scaled,
								   col=gridCol, lwd=0.5, ...)
								   
					   }
								   
								   
					# ---------------------------------------------------------------------
					# 4. plot of LOC and LOIC, and other axes
					if (suppress.surface==FALSE) {
						  for (a in axes) {
							  if (!is.null(axesList[[a]])) {
								  a0 <- RESCALE(getIntersect2(p0=axesList[[a]]$p0, p1=axesList[[a]]$p1))
								  if (nrow(a0) <= 1) break;
					              panel.3dscatter(x = a0$X, y = a0$Y, z = a0$Z, xlim = xlim, ylim = ylim, zlim = zlim,
					                      	xlim.scaled = xlim.scaled, ylim.scaled = ylim.scaled, zlim.scaled = zlim.scaled, 
											type="l", col.line=axesList[[a]]$col, lty=axesList[[a]]$lty, lwd=2, ...)
							  }
						  }   
					  }
						  
  					# ---------------------------------------------------------------------
  					# 4b. plot of maximum lines
					if (maxlines == TRUE & suppress.surface==FALSE) {
						# maximum X for a given Y
  							  a0 <- RESCALE(getIntersect2(p0=-(C[1]/C[5]), p1=-((2*C[3])/C[5])))
  				              panel.3dscatter(x = a0$X, y = a0$Y, z = a0$Z, xlim = xlim, ylim = ylim, zlim = zlim,
  				                              xlim.scaled = xlim.scaled, ylim.scaled = ylim.scaled, zlim.scaled = zlim.scaled,
  											  type="l", col.line="red", lty="dashed", lwd=2, ...)
											  
											  
											  
  							  a0 <- RESCALE(getIntersect2(p0=-(C[2]/(2*C[4])), p1=-((C[5])/(2*C[4]))))
							  #a0 <- RESCALE(getIntersect2(p0=-(C[2]/C[5]), p1=-((2*C[4])/C[5])))
  				              panel.3dscatter(x = a0$X, y = a0$Y, z = a0$Z, xlim = xlim, ylim = ylim, zlim = zlim,
  				                              xlim.scaled = xlim.scaled, ylim.scaled = ylim.scaled, zlim.scaled = zlim.scaled,
  											  type="l", col.line="blue", lty="dashed", lwd=2, ...)
											  
					}
						  
					

					# ---------------------------------------------------------------------
					# 5. Borders, front part	  
									   
					   if (border==TRUE & suppress.surface==FALSE) {
 						  # plot the front boundary lines
 			              panel.3dscatter(x = x.box[box$side==2], y = y.box[box$side==2], z = z.box[box$side==2], xlim = xlim, ylim = ylim, zlim = zlim, xlim.scaled = xlim.scaled, ylim.scaled = ylim.scaled, zlim.scaled = zlim.scaled, type="l", col.line=gridCol, lwd=3, ...)
						  panel.3dscatter(x = x.box[box$side==4], y = y.box[box$side==4], z = z.box[box$side==4], xlim = xlim, ylim = ylim, zlim = zlim, xlim.scaled = xlim.scaled, ylim.scaled = ylim.scaled, zlim.scaled = zlim.scaled, type="l", col.line=gridCol, lwd=3, ...)
					   }
										  
						  if (param == TRUE) {
							  grid::grid.text(SPs, .02, .95, just="left", gp=grid::gpar(cex=cex*0.8))
						  }  
						  
						  if (coefs == TRUE) {
							  grid::grid.text(COEFS, .80, .87, just="left", gp=grid::gpar(cex=cex*0.8))
						  }  
						  
						
						
						# ---------------------------------------------------------------------
						# 6a: The bag plot, if requested
					
						if (hull==TRUE & suppress.surface==FALSE) {	
							
							# bag (= inner bag)
							if (any(bag$X < xlim[1] | bag$X > xlim[2] | bag$Y < ylim[1] | bag$Y > ylim[2])) {
								warning("The bag is partly outside the plotting region. Bag is not displayed, please adjust xlim and ylim to include the full range of raw data.")
							} else {
		  					  bag.rescale <- RESCALE(bag)
		  		              panel.3dscatter(x = bag.rescale$X, y = bag.rescale$Y, z = bag.rescale$Z, xlim = xlim, ylim = ylim, zlim = zlim, xlim.scaled = xlim.scaled, ylim.scaled = ylim.scaled, zlim.scaled = zlim.scaled, type="l", col.line="grey30", lty="dashed", lwd=2, ...)
							
							}
							
							# loop (= outer bag)
							if (any(loop$X < xlim[1] | loop$X > xlim[2] | loop$Y < ylim[1] | loop$Y > ylim[2])) {
								warning("The loop is partly outside the plotting region. Loop is not displayed, please adjust xlim and ylim to include the full range of raw data.")
							} else {
		  					  loop.rescale <- RESCALE(loop)
		  		              panel.3dscatter(x = loop.rescale$X, y = loop.rescale$Y, z = loop.rescale$Z, xlim = xlim, ylim = ylim, zlim = zlim, xlim.scaled = xlim.scaled, ylim.scaled = ylim.scaled, zlim.scaled = zlim.scaled, type="l", col.line="black", lty="dashed", lwd=2, ...)
							
							}
						}	  	
						
						
						# ---------------------------------------------------------------------
						# 6b. Raw data points scatter plot	  
					  		if (points$show == TRUE) {
							
	  			              x2 <- xlim.scaled[1] + diff(xlim.scaled) * (x.points - xlim[1]) / diff(xlim)
	  			              y2 <- ylim.scaled[1] + diff(ylim.scaled) * (y.points - ylim[1]) / diff(ylim)
							  z2 <- zlim.scaled[1] + diff(zlim.scaled) * (z.points - zlim[1]) / diff(zlim)
							
	  			              panel.3dscatter(x = x2, y = y2, z = z2, xlim = xlim, ylim = ylim, zlim = zlim,
	  			                              xlim.scaled = xlim.scaled, ylim.scaled = ylim.scaled, zlim.scaled = zlim.scaled,
	  										  pch=20, col=points$color, cex=points$cex, ...)
							  # plot outliers
							  if (points$out.mark==TRUE) {
	  			              	panel.3dscatter(x = x2[fit$outliers], y = y2[fit$outliers], z = z2[fit$outliers], xlim = xlim, ylim = ylim, zlim = zlim,
	  			                              xlim.scaled = xlim.scaled, ylim.scaled = ylim.scaled, zlim.scaled = zlim.scaled,
	  										  pch=4, col="red", cex=points$cex, ...)
							  }
								
					  		}	
						
					# ---------------------------------------------------------------------
					# 7. plot contour lines on surface:
					
					if (contour$show == TRUE & suppress.surface==FALSE) {
						# C0 keeps the contour lines and has been computed before
							for (cL in C0$group) {
							  C1 <- RESCALE(C0[C0$group==cL, c("X", "Y", "Z")])
							  
							  if (contour$show == TRUE) {
								  	panel.3dscatter(x = C1$X, y = C1$Y, z = C1$Z, 
										xlim = xlim, ylim = ylim, zlim = zlim,
				                        xlim.scaled = xlim.scaled, ylim.scaled = ylim.scaled, zlim.scaled = zlim.scaled,
										type="l", col.line=contour$color, lty="solid", lwd=1, ...)
								}								
							}
							
							# highlight specific contour lines?
							if (length(contour$highlight) > 0) {
								C2 <- C0[C0$Z %in% f0(unique(C0$Z), contour$highlight), ]
								for (cL in C2$group) {
								  C3 <- RESCALE(C2[C2$group==cL, c("X", "Y", "Z")])
					              panel.3dscatter(x = C3$X, y = C3$Y, z = C3$Z, xlim = xlim, ylim = ylim, zlim = zlim,
					                              xlim.scaled = xlim.scaled, ylim.scaled = ylim.scaled, zlim.scaled = zlim.scaled,
												  type="l", col.line=contour$color, lty="solid", lwd=2, ...)
								}
								
							}	
						}
						
				}  # of mypanel2
			


				# local function: compute the surface line, defined by a line on the X-Y plane (p0 = intercept, p1=slope)
				getIntersect2 <- function(p0, p1, Z=NULL) {
					X <- seq(min(xlim), max(xlim), length.out=grid*2)
					Y <- p0 + p1*X
					n <- data.frame(X, Y)
					n2 <- add.variables(z~X+Y, n)
					n2$Z <- b0 + colSums(c(x, y, x2, y2, xy)*t(n2[, c(1:5)]))
					if (!is.null(Z)) n2$Z <- Z
					return(n2[, c("X", "Y", "Z")])
				}
				
				axesList <- list()
				axesList[["LOC"]]  <- list(p0=0, p1=1, lty="solid", col=LOC.col)
				axesList[["LOIC"]] <- list(p0=0, p1=-1, lty="solid", col=LOIC.col)
				if (x2 != y2) {
					axesList[["PA1"]] <- list(p0=SP$p10, p1=SP$p11, lty="dotted", col=PA1.col)
					axesList[["PA2"]] <- list(p0=SP$p20, p1=SP$p21, lty="dotted", col=PA2.col)	
				}			
				
				
				# Define color range: Relative to surface min/max, or relative to box (zlim)?
				if (pal.range == "box") {
					at <- seq(zlim[1], zlim[2], length.out=length(pal)-1)
				} else if (pal.range == "surface") {
					at <- seq(min(new2$z), max(new2$z), length.out=length(pal)-1)
				}
				
				# define the appearance of the color legend
				CK <- FALSE
				if (legend == TRUE) {
					CK <- list(labels=list(cex=cex))
				}
				
				# Define appearance of the surrounding box
				axesCol <- "black"
				boxCol <- "black"
				
				if (suppress.box == TRUE) {
					axesCol <- "transparent"
					boxCol <- NA
				}
				
			p1 <- wireframe(z ~ x*y, new2,  drape=TRUE, 
					scales 	= list(arrows = FALSE, cex=cex, col = axesCol, font = 1, tck=tck, distance=distance), 
					xlab	= list(cex=cex, label=xlab, rot=label.rotation[["x"]]), 
					ylab	= list(cex=cex, label=ylab, rot=label.rotation[["y"]]), 
					zlab	= list(cex=cex, label=zlab, rot=label.rotation[["z"]]), zlim=zlim, 
					main	= list(cex=cex.main, label=main),
					screen	= rotation,
					at		= at, col.regions=pal, colorkey=CK, 
					par.settings = list(
						axis.line = list(col = "transparent"), 
						layout.heights = list(top.padding=pad, bottom.padding=pad), 
						layout.widths=list(left.padding=pad, right.padding=pad),
						box.3d = list(col=boxCol)), 
					axes	= axes,
					axesList= axesList, 
					SPs		= SP.text,
					COEFS	= COEFS, 
					panel.3d.wireframe = mypanel2,
					x.points=xpoints, y.points=ypoints, z.points=zpoints)
				
	}  # of type == "3d"
	
	
	
	## ======================================================================
	## Contour plot
	## ======================================================================
	if (type == "contour") {
		if (!all(C == 0)) {
			
			# Define color range: Relative to surface min/max, or relative to box (zlim)?
			if (pal.range == "box") {
				limits <- c(zlim[1], zlim[2])
			} else if (pal.range == "surface") {
				limits <- c(min(new2$z), max(new2$z))
			}
			
			p1 <- ggplot(new2, aes_string(x="x", y="y", z="z")) + geom_tile(aes_string(fill="z")) + scale_fill_gradientn(zlab, colours=pal, limits=limits) + theme_bw() + theme(aspect.ratio=1) + xlab(xlab) + ylab(ylab)
			
			if (legend==FALSE) {
				p1 <- p1 + guides(fill=FALSE)
			}
	
			p1 <- p1 + stat_contour(bins=40, alpha=.4)
			
			# highlight specific contour lines?
			if (length(contour$highlight) > 0) {
				cLines <- ggplot_build(p1)
				C0 <- cLines$data[[2]][, c("x", "y", "level", "group")]

				# Find closest values in contours
				C1 <- C0[C0$level %in% f0(unique(C0$level), contour$highlight), ]
				p1 <- p1 + geom_path(data=C1, aes_string(x="x", y="y", group="group", z="level"), size=1.1)
			}
				
			# (in)congruence lines
			if ("LOC" %in% axes) {
				p1 <- p1 + geom_abline(aes(intercept=0, slope=1), color="grey20")
			}
			if ("LOIC" %in% axes) {
				p1 <- p1 + geom_abline(aes(intercept=0, slope=-1), linetype="dotted", size=1, color="grey20")
			}
			if (("PA1" %in% axes) & !any(is.na(SP[c("p10", "p11")]))) {
				p1 <- p1 + geom_abline(data=data.frame(SP[c("p10", "p11")]), aes_string(intercept="p10", slope="p11"), color="grey20")
			}
			if (("PA2" %in% axes) & !any(is.na(SP[c("p20", "p21")]))) {
				p1 <- p1+ geom_abline(data=data.frame(SP[c("p20", "p21")]), aes_string(intercept="p20", slope="p21"), linetype="dotted", color="grey20")
			}
			
			if (showSP==TRUE & !any(is.na(SP[c("X0", "Y0")])) & !model %in% c("RR", "SQD", "SSQD", "SRSQD", "SRR", "SRRR")) {
				p1 <- p1 + annotate("point", x=SP$X0, y=SP$Y0, z=max(new2$z))
			}
				
				
			if (points$show == TRUE) {
				if (points$out.mark==FALSE) {
					p1 <- p1 + annotate("point", x=xpoints, y=ypoints, color=points$color, size=3*points$cex)
				}
				if (points$out.mark==TRUE) {
					colvec <- rep(points$color, nrow(fit$data))
					colvec[fit$outliers] <- "red"
					shapevec <- rep(19, nrow(fit$data))
					shapevec[fit$outliers] <- 4
					p1 <- p1 + annotate("point", x=fit$data[, fit$IV1], y=fit$data[, fit$IV2], color=colvec, size=3*points$cex, shape=shapevec)
				}
			}
			
			if (hull==TRUE & !is.null(points$data)) {
				p1 <- p1 + annotate("path", x=bag$X, y=bag$Y, linetype="solid", size=1, color="grey10")
				p1 <- p1 + annotate("path", x=loop$X, y=loop$Y, linetype="dashed", size=1, color="grey10")
			}
			
			# plot CI of SP
			if (showSP==TRUE & showSP.CI==TRUE & !is.null(fit)) {
				PAR <- getPar(fit, "coef", model=model)
				p1 <- p1 + annotate("errorbar", x=SP$X0, y=SP$Y0, ymin=PAR[PAR$label=="Y0", "ci.lower"], ymax=PAR[PAR$label=="Y0", "ci.upper"], z=max(new2$z), width=.3)
				p1 <- p1 + annotate("errorbarh", x=SP$X0, y=SP$Y0, xmin=PAR[PAR$label=="X0", "ci.lower"], xmax=PAR[PAR$label=="X0", "ci.upper"], z=max(new2$z), height=.3)
			}
			
			# Title
			if (main != "") p1 <- p1 + ggtitle(main)
				
			p1 <- p1 + coord_cartesian(xlim=xlim, ylim=ylim)
			
		}
	}
	
	
	return(p1)
}


#' @method plot RSA
#' @export

# Purpose: Extract the model parameters, xlim, xlab, etc. from the fitted object and give it to the plotRSA function
plot.RSA <- function(x, ...) {
	fit <- x
	
	extras <- match.call(expand.dots = FALSE)$...

	if (is.null(extras)) {extras <- list()}
	if (is.null(extras[["model"]])) {extras[["model"]] <- "full"}
	
	C <- coef(fit$models[[extras$model]])
	if (fit$models[[extras$model]]@Options$estimator != "DWLS") {
		extras[["b0"]] <- as.numeric(ifelse(is.na(C[paste0(fit$DV, "~1")]), 0, C[paste0(fit$DV, "~1")]))
	} else {			
		# the threshold is the negative of the intercept ...
		extras[["b0"]] <- -as.numeric(ifelse(is.na(C[paste0(fit$DV, "|t1")]), 0, C[paste0(fit$DV, "|t1")]))
	}
	extras$x <- as.numeric(ifelse(is.na(C["b1"]), 0, C["b1"]))
	extras$y <- as.numeric(ifelse(is.na(C["b2"]), 0, C["b2"]))
	extras$x2 <- as.numeric(ifelse(is.na(C["b3"]), 0, C["b3"]))
	extras$y2 <- as.numeric(ifelse(is.na(C["b5"]), 0, C["b5"]))
	extras$xy <- as.numeric(ifelse(is.na(C["b4"]), 0, C["b4"]))
	extras$w <- as.numeric(ifelse(is.na(C["b6"]), 0, C["b6"]))
	extras$wx <- as.numeric(ifelse(is.na(C["b7"]), 0, C["b7"]))
	extras$wy <- as.numeric(ifelse(is.na(C["b8"]), 0, C["b8"]))
	
	# cubic parameters
	extras$x3 <- as.numeric(ifelse(is.na(C["b9"]), 0, C["b9"]))
	extras$xy2 <- as.numeric(ifelse(is.na(C["b10"]), 0, C["b10"]))
	extras$x2y <- as.numeric(ifelse(is.na(C["b11"]), 0, C["b11"]))
	extras$y3 <- as.numeric(ifelse(is.na(C["b12"]), 0, C["b12"]))
	
	if (is.null(extras[["xlab"]])) {extras[["xlab"]] <- fit$IV1}
	if (is.null(extras[["ylab"]])) {extras[["ylab"]] <- fit$IV2}
	if (is.null(extras[["zlab"]])) {extras[["zlab"]] <- fit$DV}
		
	extras$fit <- fit
	
	# define the defaults
	if (is.null(extras$points) || (typeof(extras$points) == "logical" && extras$points == TRUE)) {
		extras$points <- list(show=TRUE, value="raw", jitter=0, color="black", cex=.5, out.mark=FALSE)
	}
	if (is.null(extras$points) || (typeof(extras$points) == "logical" && extras$points == FALSE)) {
		extras$points <- list(show=FALSE, value="raw", jitter=0, color="black", cex=.5, out.mark=FALSE)
	}
	if (is.null(extras$points$out.mark)) extras$points$out.mark <- FALSE

	if (extras$points$out.mark == FALSE) {
		data.used <- fit$data[fit$data$out==FALSE, ]
	}
	if (extras$points$out.mark == TRUE) {
		data.used <- fit$data
	}
	
	extras$points$data <- data.used[, c(fit$IV1, fit$IV2, fit$DV, colnames(fit$data)[which(!colnames(fit$data) %in% c(fit$IV1, fit$IV2, fit$DV))])]

	adjust <- FALSE
	if (is.null(extras$xlim)) {
		extras$xlim <- c(min(data.used[, fit$IV1], na.rm=TRUE), max(data.used[, fit$IV1], na.rm=TRUE))
		# expand range by 20% at each end
		extras$xlim[1] <- extras$xlim[1]*ifelse(extras$xlim[1]<0, 1.1, 0.9)
		extras$xlim[2] <- extras$xlim[2]*ifelse(extras$xlim[2]<0, 0.9, 1.1)
		adjust <- TRUE
	}
	
	if (is.null(extras$ylim)) {
		extras$ylim <- c(min(data.used[, fit$IV2], na.rm=TRUE), max(data.used[, fit$IV2], na.rm=TRUE))
		extras$ylim[1] <- extras$ylim[1]*ifelse(extras$ylim[1]<0, 1.1, 0.9)
		extras$ylim[2] <- extras$ylim[2]*ifelse(extras$ylim[2]<0, 0.9, 1.1)
		adjust <- TRUE
	}
	
	if (adjust == TRUE) {
		extras$xlim[1] <- extras$ylim[1] <- min(extras$xlim[1], extras$ylim[1])
		extras$xlim[2] <- extras$ylim[2] <- max(extras$xlim[2], extras$ylim[2])
	}
	
	do.call(plotRSA, as.list(extras), envir = parent.frame())
}

