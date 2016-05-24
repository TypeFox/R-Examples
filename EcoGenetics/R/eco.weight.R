#' Spatial weights
#' 
#' @description Spatial weights for individuals with coordinates XY
#' @param XY Matrix/data frame with projected coordinates.
#' @param method Method of spatial weight matrix: "circle", "knearest", "inverse", 
#' "circle.inverse", "exponential", "circle.exponential".
#' @param k Number of neighbors for nearest neighbor distance. When equidistant
#' neighbors are present, the program select them randomly.
#' @param d1 Minimum distance for circle matrices.
#' @param d2 Maximum distance for circle matrices.
#' @param p Power for inverse distance. Default = 1.
#' @param alpha Alpha value for exponential distance. Default = 1.
#' @param dist.method Method of computing distance when XY is in metric
#' units. If latlon is TRUE, the method is euclidean. 
#' @param row.sd Logical. Should be row standardized the matrix? Default FALSE 
#' (binary weights).
#' @param self Should be the individuals self-included in circle or knearest
#' weights? Defalut FALSE.
#' @param latlon Are the coordinates in decimal degrees format? Defalut FALSE. If TRUE,
#' the coordinates must be in a matrix/data frame with the longitude in the first
#' column and latitude in the second. The position is projected onto a plane in
#' meters with the function \code{\link[SoDA]{geoXY}}.
#' @details This program computes a weights matrix (square matrix with individuals
#' in rows and columns, and weights wij in cells (i and j, individuals)) 
#' under the following available methodologies: 
#'   
#' - circle: all the connection between individuals i and j, included 
#' in a distance radius, higher than d1 and lower than d2, with center in 
#' the individual i, have a value of 1 for binary weights.
#' This distance requires the parameters d1 and d2 (default d1 = 0).
#'  
#' - knearest: the connections between an individual and its 
#' nearest neighbors of each individual i have a value of 1  for binary weights.
#' This distance requires the parameter k.
#' 
#' - inverse: inverse distance with exponent p (distance = 1/dij^p, with
#' dij the distance between individuals i and j).
#' This distance requires the parameter p (default p = 1).
#'   
#' - circle inverse: combination of "circle" and "inverse". 
#'  It is the matrix obtained by multiplying each element in a "circle" 
#'  binary matrix, and an "inverse" matrix.
#'  This distance requires the parameters p, d1 and d2 (default p = 1, d1 = 0).
#'  
#'  - exponential: inverse exponential distance with parameter alpha
#'   (distance = 1/e^(alpha *dij), with dij the distance between individuals i and j).
#'   This distance requires the parameter alpha (default alpha = 1).
#'  
#' - circle exponential: combination of "circle" and "exponential". 
#'  It is the matrix obtained by multiplying each element in a "circle" 
#'  binary matrix, and an "exponential" matrix.
#'  This distance requires the parameters alpha, d1 and d2 (default alpha = 1, d1 = 0).
#'  
#'  In row standardization, each weight wij for the individual i, is divided by the
#'  sum of the row weights (i.e., wij / sum(wij), where sum(wij) is computed over an
#'  individual i and all individuals j). 
#'  
#'  When self is TRUE, the connection j = i is also included.
#'  
#' @return An object of class eco.weight with the following slots:
#' @return > W weights matrix
#' @return > XY input coordinates
#' @return > METHOD weights construction method
#' @return > PAR parameters used for the construction of weights
#' @return > PAR.VAL values of the parameters used for the construction of weights
#' @return > ROW.SD row standardization (logical)
#' @return > SELF data self-included (logical)
#' @return > NONZERO percentage of non-zero connections
#' @return > NONZEROIND percentage of individuals
#' with non-zero connections 
#' @return > AVERAGE average number of connection per individual
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
#' 
#' \dontrun{
#' 
#' data(eco.test)
#' 
#' # "circle" method
#' con <- eco.weight(eco[["XY"]], method = "circle", d1 = 0, d2 = 2)
#' con
#' 
#' # "knearest" method
#' con <- eco.weight(eco[["XY"]], method = "knearest", k=3)
#' con
#' 
#' # "inverse" method
#' con <- eco.weight(eco[["XY"]], method = "inverse")
#' con
#' 
#' # "circle.inverse" method
#' con <- eco.weight(eco[["XY"]], method = "circle.inverse", d2 = 2)
#' con
#' 
#' # "exponential" method
#' con <- eco.weight(eco[["XY"]], method = "exponential")
#' con
#' 
#' # "circle.exponential" method
#' con <- eco.weight(eco[["XY"]], method = "circle.exponential", d2 = 2)
#' con
#' 
#' 
#' #-----------------------
#' # ACCESSORS USE EXAMPLE
#' #-----------------------
#' 
#' # the slots are accessed with the generic format 
#' # (ecoslot. + name of the slot + name of the object). 
#' # See help("EcoGenetics accessors")
#' 
#' ecoslot.METHOD(con)        # slot METHOD
#' ecoslot.PAR(con)           # slot PAR
#' ecoslot.PAR.VAL(con)       # slot PAR.VAL
#' 
#' }
#' 
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @export

setGeneric("eco.weight", function(XY,
											 method = c("circle", "knearest", "inverse",  
											            "circle.inverse", "exponential", 
											            "circle.exponential"),
											 d1 = 0, d2 = NULL,  k = NULL,  p = 1, alpha = 1, 
											 dist.method = "euclidean",
											 row.sd = FALSE, self = FALSE, latlon = FALSE) {
	
	method <- match.arg(method)
	dist.method <- match.arg(dist.method)
	
	#distance configuration
if(latlon == FALSE) {
		distancia <- as.matrix(dist(XY), upper = T, method = dist.method)
	} else {
	  distancia <- dist(SoDA::geoXY(XY[,2], XY[,1], unit=1))
	  distancia <- as.matrix(distancia, upper = T, method = dist.method)
	}
	#####
	
	if(method == "inverse") {
		y <- as.matrix(dist(XY, upper = T, method = dist.method))
		y <- 1 / (y ^ p)
		diag(y) <- 0
		
		param <- "p"
		param.values <- p
		
	} else if(method == "circle") {
	
	   if(is.null(d2)) {
	   stop("a d2 argument must be given")
	 }
		temp <- which((distancia <= d2) & (distancia > d1))
		y <- distancia
		y <- y - distancia
		y[temp] <- 1
		if(self) {
		  diag(y) <- 1
		}
		
		param <- c("d1", "d2")
		param.values <- c(d1, d2)
	
} else if(method == "circle.inverse") {
		if(is.null(d2)) {
      stop("the argument d2 is missing")
		}
    
		temp <- which((distancia <= d2) & (distancia >= d1))
		dummy <- distancia - distancia
		dummy[temp] <- 1
		y <- 1 / (distancia ^ p)
		y <- dummy * y
		diag(y) <- 0
		
		param <- c("d1", "d2", "p")
		param.values <- c(d1, d2, p)
		
	} else if(method == "exponential") {
		y <- 1 / exp(alpha * distancia)
		diag(y) <- 0
		
		param <- "alpha"
		param.values <- alpha
		
	} else if(method == "circle.exponential") {
	 
	  if(is.null(d2)) {
	    stop("the argument d2 is missing")
	  }
		
		temp <- which((distancia <= d2) & (distancia > d1))
		dummy <- distancia - distancia
		dummy[temp] <- 1
		y <- 1 / exp(alpha * distancia)
		y <- dummy * y
		diag(y) <- 0
		
		param <- c("d1", "d2", "alpha")
		param.values <- c(d1, d2, alpha)
		
	} else if(method == "knearest") {
	  
	  if(is.null(k)) {
	    stop("the argument k is missing")
	  }
		
		y <- t(apply(distancia, 1, function(x) {
			x[x != 0] <- rank(x[x != 0], ties.method = "random")
			return(x)}))
		
		y <- t(apply(y, 1, function(x) as.numeric(x <= k & x != 0)))
		if(self) {
		  diag(y) <- 1
		}
		
		param <- "k"
		param.values <- k
		
	} 
	
	#row standardization
	if(row.sd) {
	  y <- y/apply(y, 1, sum)
	  y[is.na(y)] <- 0
	} 

	#output construction
	out <- new("eco.weight")
	out@METHOD <- method
	out@ROW.SD <- row.sd
	if(method == "circle" | method =="knearest") {
	out@SELF <- self
	} else {
	out@SELF <- NA
	}
	out@W <- y
	out@XY <- XY
	
	y2 <- y
	diag(y2) <- 0
	out@NONZERO <- round(100* sum(y2 != 0) / (nrow(y2)^2 - nrow(y2)), 1)
	out@NONZEROIND <- round(100 * sum(apply(y2, 1, sum) != 0) / nrow(y2), 1)
	out@AVG <- round(sum(apply(y2, 1, sum))  / nrow(y2), 1)
	out@PAR <- param
	out@PAR.VAL <- param.values
	
	out

})
