#' Obtention of a list of spatial weights for classes
#' defined by inter-individual distances or nearest-neighbors
#' 
#' @description This program returns a list of weights matrices (binary or
#' row-standardized), one for each spatial class. For a given maximum and minimum 
#' inter-individual distance (IID), the data can be partitioned in different ways. 
#' The program set as default the highest IID as the maximum, 
#' and the lowest as the minimum. 
#' These values can be changed with "smax" and "smin", respectively. 
#' Intervals may be generated with the parameters "int" (which divides the range 
#' each int distance units), "nclass" (which divides the range in n-classes) and "size" (a fixed
#' size of pairs included in each class). When a partition argument is not given (int, nclass or size) the program
#' determines the number of classes using the Sturge's rule (default) or 
#' the Freedman- Diaconis method. Two additional methods can be used:
#' a list with nearest-neighbors matrices, from 1 to k nearest-neighbors, may be generated
#' with the argument "kmax".
#' A custom vector with breaks may be provided by the user with the argument "seqvec".
#' See the examples. 
#' 
#' @param XY Matrix/data frame with projected coordinates.
#' @param int Distance interval in the units of XY.
#' @param smin Minimum class distance in the units of XY.
#' @param smax Maximum class distance in the units of XY.
#' @param kmax Number of nearest-neighbors.
#' @param nclass Number of classes.
#' @param seqvec Vector with breaks in the units of XY.
#' @param size Number of individuals per class.
#' @param bin Rule for constructing intervals when a partition parameter (int, 
#' nclass or size) is not given. Default is Sturge's rule (Sturges, 1926). Other
#' option is Freedman-Diaconis method (Freedman and Diaconis, 1981).
#' @param cummulative Logical. Should be created matrices considering cummulative
#' distances instead of discrete classes? Default FALSE.
#' @param row.sd Logical. Should be row standardized each matrix? Default FALSE 
#' (binary weights).
#' @param self Logical. Should be included a first matrix with ones in the diagonal and zeros
#' zeros in the other cells? Default FALSE.
#' @param latlon Are the coordinates in decimal degrees format? Defalut FALSE. If TRUE,
#' the coordinates must be in a matrix/data frame with the longitude in the first
#' column and latitude in the second. The position is projected onto a plane in
#' meters with the function \code{\link[SoDA]{geoXY}}.
#' 
#' @return The program returns an object of class "eco.lagweight" with the following
#' slots:
#' @return > PAR parameters used for the construction of breaks
#' @return > PAR.VAL values of the parameters used for the construction of breaks
#' @return > ROW.SD row standardization (logical)
#' @return > SELF data self-included (logical)
#' @return > W weights list
#' @return > XY coordinates
#' @return > MEAN mean class distances
#' @return > LOGMEAN mean of the class distances logarithm
#' @return > CARDINAL number of elements in each class
#' @return > BREAKS breaks
#' @return > METHOD breaks construction method
#' 
#' 
#' \strong{ACCESS TO THE SLOTS}
#' The content of the slots can be accessed 
#' with the corresponding accessors, using
#' the generic notation of EcoGenetics 
#' (<ecoslot.> + <name of the slot> + <name of the object>).
#' See help("EcoGenetics accessors") and the Examples
#' section below.
#' 
#' @examples
#' \dontrun{
#' data(eco.test)
#' 
#' # method sturges-smax: in this case, the program generates 
#' # classes using the Sturge's rule.  
#' # As smax and smin are undefined, the program uses the default
#' # options (smin = 0, and smax = maximum inter-individual distance)
#' classlist <- eco.lagweight(eco[["XY"]]) 
#' classlist
#' 
#' # method sturges-smax: idem, but smax = 16
#' classlist <- eco.lagweight(eco[["XY"]], smax=16) 
#' 
#' ## using smax <16 in this case generates empty classes, 
#' ## which is not allowed
#' 
#' # method sturges-smax: idem, but smin = 3
#' classlist <- eco.lagweight(eco[["XY"]], smin = 3, smax = 15)
#' 
#' #-----------------------
#' # ACCESSORS USE EXAMPLE
#' #-----------------------
#' 
#' # the slots are accessed with the generic format 
#' # (ecoslot. + name of the slot + name of the object). 
#' # See help("EcoGenetics accessors")
#' 
#' ecoslot.BREAKS(classlist) # information about breaks. It includes the upper and lower limits
#' 
#' # method sturges-smax: complete range, 
#' # and cummulative = TRUE (instead of using
#' # lower and upper limits for each class, only the upper is used in turn)
#' classlist <- eco.lagweight(eco[["XY"]], cummulative = TRUE)
#' 
#' # method n.classes-smax: complete range partitioned in 4 classes
#' classlist <- eco.lagweight(eco[["XY"]], nclass = 4)
#' 
#' # method n.classes-smax: idem, but smax =  15
#' classlist <- eco.lagweight(eco[["XY"]], nclass = 4, smax = 15)
#' 
#' # method int-smax: the complete range partitioned each <int> units
#' # of inter-individual distance
#' classlist <- eco.lagweight(eco[["XY"]], int = 2)
#' 
#' # method int-smax: idem, but smax = 15 and smin = 3
#' classlist <- eco.lagweight(eco[["XY"]], int = 2, smin = 3, smax = 15)
#' 
#' # method equal.size: n individuals in each class, 
#' # partitioning the complete range.
#' classlist <- eco.lagweight(eco[["XY"]], size = 1000)
#' 
#' ## In the latter example, as an inter-individual distance
#' ## appear more than one time (different individuals pairs, 
#' ## identical distances), with a size <700 the limits
#' ## of some classes cannot be defined, and this is not allowed
#' 
#' # method equal.size: n individuals in each class, 
#' # but smax = 15
#' classlist <- eco.lagweight(eco[["XY"]], size = 1000, smax = 15)
#' 
#' # method kmax: sequence from k = 1 to k = n, in this case, n = 3
#' classlist <- eco.lagweight(eco[["XY"]], kmax = 3)
#' 
#' # method kmax: idem, but elements self-included
#' # (i.e., the pairs i-i, for all individuals i, are included)
#' classlist <- eco.lagweight(eco[["XY"]], kmax = 3, self = TRUE)
#' 
#' # method seqvec: a vector with the breaks is used
#' vec <- seq(0, 10, 2)
#' classlist <- eco.lagweight(eco[["XY"]], seqvec = vec)
#' }
#' 
#' @references 
#' 
#' Freedman D., and P. Diaconis. 1981. On the histogram as a density estimator: 
#' L 2 theory. Probability theory and related fields, 57: 453-476.
#' 
#' Sturges  H. 1926. The choice of a class interval. Journal of the American 
#' Statistical Association, 21: 65-66.
#' 
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' 
#' @export

setGeneric("eco.lagweight", 
          function(XY, 
                   int = NULL, 
          				 smin = 0,
          				 smax = NULL,
          				 kmax = NULL,
          				 nclass = NULL,
          				 seqvec = NULL,
          				 size = NULL,
          				 bin = c("sturges", "FD"),
          				 cummulative = FALSE,
					 				 row.sd = FALSE,
					 				 self = FALSE,
					 				 latlon = FALSE) {
  
  bin <- match.arg(bin)				 	

  
	#variables definitions
  if(ncol(XY) != 2) {
    message("more than 3 columns in coordinates. The first two will we
								taken as X-Y data for estimating distance intervals")
    XY <- XY[,1:2]
  }
  
	  if(latlon == FALSE) {
	    distancia <- dist(XY)
	  } else {
	    distancia <- dist(SoDA::geoXY(XY[,2], XY[,1], unit=1))
	  }
  	distancia <- as.matrix(distancia)
  	logdistancia <- log(distancia)

#method control
match.control <- sum(!is.null(smax), !is.null(kmax), !is.null(seqvec))

if(match.control == 0) {
	smax <- max(distancia)
} else if (match.control != 1) {
	stop("Only one of smax, kmax, nclass or seqvec should be given")
}


	meandist <- vector()
  logdist <- vector()
	cardinal <- vector()
	laglw <- list()
	j <- 1
	
	
	#-----------computation of lag matrices
	
	###computation based in distance
	
	#based on distance between individuals, different size
	if(is.null(kmax) & is.null(size)) {
	  
	  
	  input <- int.break(XY = XY, 
	                     int = int, 
	                     smin = smin,
	                     smax = smax,
	                     nclass = nclass,
	                     seqvec = seqvec,
	                     latlon = latlon,
	                     bin = bin)
	  breakpoints <- input$breakpoints
	  method <- input$method
	  
	  #  intervals (a, b]
	  
	  for(i in 2:length(breakpoints)) {
	    #cummulative distance
	    if(cummulative) {
	      temp <- which(distancia <= breakpoints[i] & (distancia > smin))
	    } else {
	      temp <- which((distancia <= breakpoints[i]) & (distancia > breakpoints[i-1]))
	    }
	    
	    distemp <- distancia[temp]
	    logdistemp <- logdistancia[temp]
	    meandist[j] <- mean(distemp)
	    logdist[j] <- mean(logdistemp)
	    dummy <- distancia
	    dummy <- dummy - distancia
	    dummy[temp] <- 1
	    laglw[[j]] <- dummy
	    cardinal[j] <- sum(dummy) / 2
	    j <- j+1
	  }
	  
	  if(self) {
	    dummy.self <- distancia - distancia
	    diag(dummy.self) <- 1
	    cardinal <- c(length(diag(dummy.self)), cardinal)
	    dummy.self <- list(dummy.self)
	    laglw <- append(dummy.self, laglw)
	    meandist <- c(0, meandist)
	    
	  }
	} 
	
	#based on distance between individuals, equal size
	if(is.null(kmax) & !is.null(size)) {
	  size2 <- 2 * size
	  method <- "equal.size"
	  vec.mat <- as.vector(distancia)
	  largo <- length(vec.mat)
	  names(vec.mat) <- 1:largo
	  vec.mat.nodiag <- vec.mat[vec.mat != 0]
	  vec.sort <- sort(vec.mat.nodiag)
	  cortes <- seq(size2, length(vec.mat.nodiag), size2)
	  cortes <- c(0, cortes)
	  breakpoints <- 0
	  
	  for(i in 2:length(cortes)) {
	    #cummulative distance
	    if(cummulative) {
	      temp <- as.numeric(names(vec.sort[1:cortes[i]]))
	    } else {
	      temp <- as.numeric(names(vec.sort[(cortes[i - 1] +1):cortes[i]]))
	    }
	    
	    distemp <- distancia[temp]
	    logdistemp <- logdistancia[temp]
	    control.corte <- max(distemp)
	    
	    #control
	    if(control.corte > smax) {
	      break
	    }
	    
	    if(i >2) {
	      temp0 <- as.numeric(names(vec.sort[(cortes[i - 2] +1):cortes[i-1]]))
	      distemp0 <- distancia[temp0]
	      control.corte0 <- max(distemp0)
	      if(control.corte0 == control.corte) {
	        stop(paste("size not appropiated, max distance non different between consecutive breaks.
	                   Increase size"))
	      }
	      }
	    ##
	    
	    breakpoints[j] <- control.corte
	    meandist[j] <- mean(distemp)
	    logdist[j] <- mean(logdistemp)
	    dummy <- distancia
	    dummy <- dummy - distancia
	    dummy[temp] <- 1
	    laglw[[j]] <- dummy
	    cardinal[j] <- sum(dummy) / 2
	    j <- j+1
	    }
	  breakpoints <- c(min(distancia), breakpoints)
	  
	  if(self) {
	    dummy.self <- distancia - distancia
	    diag(dummy.self) <- 1
	    cardinal <- c(length(diag(dummy.self)), cardinal)
	    dummy.self <- list(dummy.self)
	    laglw <- append(dummy.self, laglw)
	    meandist <- c(0, meandist)
	    
	  }
	}
	
	
	#case k neighbors
	
	if(!is.null(kmax)) {
	  #computation based in k max
	  method <- "kmax"
	  smin <- NULL
	  cummulative <- TRUE
	  if(class(XY) == "dist") {
	    stop("XY is a distance matrix. kmax require a matrix XY with coordinates")
	  }
	  for(i in 1:kmax) {
	    laglw[[i]] <- (eco.weight(XY, method = "knearest", k=i))@W
	    npair <- sum(laglw[[i]])
	    meandist[i] <- sum(distancia * laglw[[i]]) / npair
	    logdist[i] <- sum(logdistancia * laglw[[i]]) / npair
	    cardinal[i] <- sum(laglw[[i]]) / 2
	  }
	}
	
	#control
	if(any(cardinal == 0)) {
	  stop("empty classes. Change parameters setting")
	}
	  
	#conversion to row standardized weights
	if(row.sd) {
	  laglist <- list()
	  laglist <- lapply(laglw, function(y) y/apply(y, 1, sum))
	  for(i in 1:length(laglist)) {
	    laglist[[i]][is.na(laglist[[i]])] <- 0
	  }
	} else {
	  laglist <-laglw
	}
	
	if(method != "kmax") {
	  if(!self){
	    breaks <- breakpoints 
	  } else {
	    breaks <- c(0, breakpoints)
	  }
	} else {
	  breaks <- 1:kmax
	}

	#parameters for seqvec
	
	if(!is.null(seqvec)) {
	  smax <- max(seqvec)
	  smin <- min(seqvec)
	  nclass <- length(seqvec) - 1
	}
	  
	param <- c("int", "smin", "smax", "kmax", "nclass", "size")
	param.val <- c(is.null(int), is.null(smin), is.null(smax), 
	               is.null(kmax), is.null(nclass), is.null(size))
	cuales <- which(!param.val)
	PAR <- param[cuales]
	PAR.VAL <- c(int, smin, smax, kmax, nclass, size)
	
	
	res <- new("eco.lagweight")
	res@W <- laglist
	res@XY <- data.frame(XY)
	res@PAR <- PAR
	res@PAR.VAL <- PAR.VAL
	res@ROW.SD <- row.sd
	res@SELF <- self
	res@CUMMUL <- cummulative
	res@MEAN <- meandist
	res@LOGMEAN <- logdist
	res@CARDINAL <- cardinal
	res@BREAKS <- breaks
	res@METHOD <- method
	
	
	res
	
          })

