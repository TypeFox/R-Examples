#' Empirical variogram
#' 
#' @description This program computes the empirical variogram of a selected variable.
#' If the coordinates are given in decimal degrees, set latlon = TRUE. The program
#' return a table with the mean class distances (d.mean) and the semivariances (obs)
#' for each class.
#' 
#' @param Z Vector for the analysis.
#' @param XY Data frame or matrix with individual's position (projected coordinates).
#' @param int Distance interval in the units of XY.
#' @param smin Minimum class distance in the units of XY.
#' @param smax Maximum class distance iin the units of XY.
#' @param nclass Number of classes.
#' @param seqvec Vector with breaks in the units of XY.
#' @param size Number of individuals per class.
#' @param bin Rule for constructing intervals when a partition parameter (int, 
#' nclass or size) is not given. Default is Sturge's rule (Sturges, 1926). Other
#' option is Freedman-Diaconis method (Freedman and Diaconis, 1981).
#' @param row.sd Logical. Should be row standardized the matrix? Default FALSE 
#' (binary weights).
#' @param latlon Are the coordinates in decimal degrees format? Defalut FALSE. If TRUE,
#' the coordinates must be in a matrix/data frame with the longitude in the first
#' column and latitude in the second. The position is projected onto a plane in
#' meters with the function \code{\link[SoDA]{geoXY}}.
#' 
#' @return The program returns an object of class "eco.correlog" 
#' with the following slots:
#' @return > OUT analysis output
#' @return > IN analysis input data
#' @return > BEAKS breaks
#' @return > CARDINAL number of elements in each class
#' @return > DISTMETHOD method used in the construction of breaks
#' 
#' 
#' \strong{ACCESS TO THE SLOTS}
#' The content of the slots can be accessed 
#' with the corresponding accessors, using
#' the generic notation of EcoGenetics 
#' (<ecoslot.> + <name of the slot> + <name of the object>).
#' See help("EcoGenetics accessors") and the Examples
#' section below
#' 
#' 
#' @examples
#' \dontrun{
#' 
#' data(eco.test)
#' variog <- eco.variogram(Z = eco[["P"]][, 2],XY =  eco[["XY"]])
#' plot(variog)
#' 
#' # variogram plots support the use of ggplot2 syntax
#' variogplot <- plot(variog) + theme_bw() + theme(legend.position="none")
#' variogplot
#' 
#' #-----------------------
#' # ACCESSORS USE EXAMPLE
#' #-----------------------
#' 
#' # the slots are accessed with the generic format 
#' # (ecoslot. + name of the slot + name of the object). 
#' # See help("EcoGenetics accessors")
#' 
#' ecoslot.OUT(variog)        # slot OUT
#' ecoslot.BREAKS(variog)     # slot BREAKS
#' 
#'}
#'
#' @references  
#' Borcard D., F. Gillet, and P. Legendre. 2011. Numerical ecology with R. 
#' Springer Science & Business Media.
#' 
#' Legendre P., and L. Legendre. 2012. Numerical ecology. Third English edition.
#' Elsevier Science, Amsterdam, Netherlands.
#' 
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' 
#' @export

setGeneric("eco.variogram",  
					 function(Z, XY, 
					         int = NULL,
					 				 smin = 0,
					 				 smax = NULL,
					 				 nclass = NULL,
					 				 seqvec = NULL,
					 				 size = NULL,
					 				 bin = c("sturges", "FD"),
					 				 row.sd = FALSE,
					 				 latlon = FALSE) {
	
	bin <- match.arg(bin)
	
	#CHECKING XY DATA
	
	if(ncol(XY) > 2) {
	  message("XY slot with > 2 columns. The first two are taken as X-Y coordinates")
	  XY <- XY[,1:2]
	} 
	
	if(latlon == TRUE) {
	  XY <- SoDA::geoXY(XY[,2], XY[,1], unit=1)
	} 
	
					 	
	####
	
	mat <- as.matrix(dist(Z))
	
listaw <- eco.lagweight(XY, 
							int = int, 
							smin = smin,
							smax = smax, 
							nclass = nclass,
							size = size,
							seqvec = seqvec,
							row.sd = row.sd,
							bin = bin)

	wg <- listaw@W
	
breakpoints <- listaw@BREAKS
d.max <- round(breakpoints[-1], 3)
d.min <- round(breakpoints[-length(breakpoints)], 3)
classint <- listaw@MEAN
classint <- round(classint, 3)
cardinal <- listaw@CARDINAL
	
	dist.dat<-paste("d=", d.min, "-", d.max)

  d.mean <- listaw@MEAN

	mat2 <- mat ^ 2
  wsub <- (2 * sapply(wg, sum))
	est <- sapply(wg, function(x) sum(x * mat2)) / wsub 
 

		tab <- data.frame(matrix(nrow = length(dist.dat), ncol = 2))
		rownames(tab) <- dist.dat
		tab[, 1] <- d.mean
		tab[, 2] <- est
		colnames(tab) <- c("d.mean","obs")

salida <- new("eco.correlog")

salida@OUT <- list(tab)
salida@IN <- list(XY = XY, Z = Z)
salida@BREAKS <- breakpoints
salida@CARDINAL <- cardinal
salida@METHOD <- "empirical variogram"
salida@DISTMETHOD <- listaw@METHOD
salida@TEST <- "none"

salida

})
		
