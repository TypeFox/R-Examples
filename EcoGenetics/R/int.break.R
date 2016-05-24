#' breaks obtention
#' 
#' @param XY data frame or matrix with individual's position.
#' @param int Distance interval in the units of the @@XY slot data.
#' @param smin Minimum class distance in the units of the @@XY slot data.
#' @param smax Maximum class distance in the units of the @@XY slot data.
#' @param kmax Number of nearest neighbors.
#' @param nclass Number of classes.
#' @param seqvec Vector with breaks in the units of the @@XY slot data.
#' @param size Number of individuals per class.
#' @param bin Rule for constructing intervals when a partition parameter (int, 
#' nclass or size) is not given. Default is Sturge's rule (Sturges, 1926). Other
#' option is Freedman-Diaconis method (Freedman and Diaconis, 1981).
#' @param latlon Are the coordinates in decimal degrees format? Defalut FALSE. If TRUE,
#' the coordinates must be in a matrix/data frame with the longitude in the first
#' column and latitude in the second. The position is projected onto a plane in
#' meters with the function \code{\link[SoDA]{geoXY}}.
#' 
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' 
#' @keywords internal


int.break <- function(XY, 
											int = NULL, 
											smin = 0,
											smax = NULL,
											seqvec = NULL,
										  nclass = NULL,
										  size = NULL,
										  bin = c("FD", "sturges"),
										  latlon = FALSE) {
	
	#-----some checking of the data-----------------------#
	
	#filtering.
	 control <- sum(!is.null(smax), !is.null(seqvec))
	 if(control > 1) {
	 	stop("only one of smax or seqvec must be selected")
	 }
	 
	 control <- sum(!is.null(int), !is.null(nclass), !is.null(seqvec))
	 if(control > 1) {
	 	stop("only one of int, nclass or seqvec must be selected")
	 }
	 
	 #matching and distance matrix obtention.
    bin <- match.arg(bin)

	
		if(latlon == FALSE) {
			distancia <- dist(XY)
		} else {
		  if(ncol(XY != 2)) {
		    message("more than 3 columns in coordinates. The first two will we
								taken as X-Y data for estimating distance intervals")
		    XY <- XY[, 1:2]
		  }
			distancia <- dist(SoDA::geoXY(XY[, 2], XY[, 1], unit = 1))
			}
	

		#controlling case defined by interval. We defined above "distancia".
		if(!is.null(int)) {
			if(is.null(smax)) {
				smax <- max(distancia)
			} 
			#we defined 5 individuals as the minimum size for initial and terminal classes
			hmuch <- sum(distancia > smin & distancia <= smin + int)
			if(hmuch < 5) {
				stop("Scale not apropiated, few individuals in initial classes.Increase int")
			}
			
			hlast <- sum(distancia > smax - int)
			if(hlast <5 ) {
				stop("Range not apropiated, few individuals in higher classes. Decrease smax 
						 or increase int")
			}
		}
		
	  #---end of data checking----------------------------------#
	 
		
		#----automatic bin selection methods----------------------#
	 
		#sturges rule
		sturges <- function(max.pair) {
		n <- sum(distancia > smin & distancia <= max.pair)
		nclases <- ceiling(1 + log(n, base = 2))
		min.pair <- smin
		ancho  <- (max.pair - smin) / nclases
		bp <- seq(min.pair, max.pair, ancho)
		bp
		}
		
		#freedman diaconis rule
		FD <- function(data) {
			nclasses <- nclass.FD(data)
			min.pair <- min(data[data != 0])
			max.pair <- max(data)
			ancho  <- (max.pair - min.pair) / nclasses
			bp <- seq(min.pair, max.pair, ancho)
			bp
		}
		
		#size partition function
	 sizefun <- function(maxarg) {
	   distarg <- as.matrix(distancia)
	   j <- 1
	   size2 <- 2 * size
	   method <- "equal.size"
	   vec.mat <- as.vector(distarg)
	   largo <- length(vec.mat)
	   names(vec.mat) <- 1:largo
	   vec.mat.nodiag <- vec.mat[vec.mat != 0]
	   vec.sort <- sort(vec.mat.nodiag)
	   cortes <- seq(size2, length(vec.mat.nodiag), size2)
	   cortes <- c(0, cortes)
	   bp <- 0
	   
	   for(i in 2:length(cortes)) {
	     temp <- as.numeric(names(vec.sort[(cortes[i - 1] +1):cortes[i]]))
	     distemp <- distarg[temp]
	     control.corte <- max(distemp)
	     
	     #control
	     if(control.corte > maxarg) {
	       break
	     }
	     
	     if(i >2) {
	       temp0 <- as.numeric(names(vec.sort[(cortes[i - 2] +1):cortes[i-1]]))
	       distemp0 <- distarg[temp0]
	       control.corte0 <- max(distemp0)
	       if(control.corte0 == control.corte) {
	         stop(paste("size not appropiated, max distance non different between consecutive breaks.
                 Increase size"))
	       }
	     }
	     ##
	     
	     bp[j] <- control.corte
	     j <- j+1
	   }
	   bp <- c(min(distarg), bp)
	   bp
	 }
		#-----------End of bin selection method--------#
		
		
		#---interval construction----------------------#
	 
		#case for smax not null
			if(!is.null(smax)) {
				
				if(!is.null(int)) {
					#case defined by interval + smax
					breakpoints <- seq(smin, smax, int)
					method <- "int - smax"
					
				# case n classes + smax
				} else if(!is.null(nclass)) {
					ancho  <- (smax - smin) / nclass
					breakpoints <- seq(smin, smax, ancho)
					method <- "n.classes-smax"
					
				#case defined by size
				} else if (!is.null(size)) {
					
				 breakpoints <- sizefun(smax)
				 method <- "equal.size-smax"
				 
					#case defined by FD or Sturges rule + smax
					} else {
					distancia2 <- distancia[distancia > smin & distancia <= smax]
					#case max == smax, interval == "nclass"
					if(!is.null(nclass)) {
						ancho  <- (smax - smin) / nclass
						breakpoints <- seq(smin, smax, ancho)
						method <- "n.classes-smax"
						
						#case max == smax, interval == "FD"
					} else  if(bin == "FD") {
					  breakpoints <- FD(distancia2)
					  method <- "FD-max"
						
						#case max == smax, interval == "sturges"
					} else {
						breakpoints <- sturges(smax)
						method <- "sturges-smax"
					}
				}
					
			#cases for smax null
			} else if(!is.null(seqvec)) {
				# first, breakpoints defined by predefined sequence
				breakpoints <- seqvec
				method <- "seqvec"
				
			} else if(!is.null(nclass)) {
				max.pair <- max(distancia)
				ancho  <- (max.pair - smin) / nclass
				breakpoints <- seq(smin, max.pair, ancho)
				method <- "n.classes-complete range"
				
			}	else if (!is.null(size)) {
				
				breakpoints <- sizefun(max(distancia))
				method <- "size-complete range"
				
			} else {
	    	
		#case defined by FD or Sturges rule 
				if(bin == "FD") {
					breakpoints <- FD(distancia)
					method <- "FD-complete range"
				} else {
				smax <- max(distancia)
				breakpoints <- sturges(smax)
				method <- "sturges-complete range"
			}
			}
	 #---end of interval construction------#
	 
		out <- list(breakpoints = round(breakpoints, 5), 
												method = method)
		out
	}
		