#' Exporting an ecogen genetic data frame into SPAGeDI format
#' 
#' @description This function converts the genetic data of an ecogen object 
#' in a SPAGeDI input file.  
#' When distance classes are required, they can be constructed by combining 
#' the parameters "int", "smin", "smax", "nclass", "seqvec" and "size", as
#' described in the function \code{\link{eco.lagweight}}.
#' A distance matrix can also be included using the "distmat" parameter.
#' Missing data must be coded as a single "NA" in the G data frame. 
#' @param eco Object of class "ecogen".
#' @param pop The name of the S slot column with the groups 
#' for the output data. The default option includes all the individuals into 
#' a single group.
#' @param ndig Number of digits coding each allele
#'  (e.g., 1: x, 2: xx, or 3: xxx). 
#' @param name The name of the output file.
#' @param int Distance interval in the units of the XY slot data.
#' @param smin Minimum class distance in the units of the XY slot data.
#' @param smax Maximum class distance in the units of the XY slot data.
#' @param nclass Number of classes.
#' @param seqvec Vector with breaks in the units of the XY slot data.
#' @param size Number of individuals per class.
#' @param bin Rule for constructing intervals when a partition parameter (int, 
#' nclass or size) is not given. Default is Sturge's rule (Sturges, 1926). Other
#' option is Freedman-Diaconis method (Freedman and Diaconis, 1981).
#' @param distmat Distance matrix to include (optional).
#' @param latlon Are the coordinates in decimal degrees format? Defalut FALSE. If TRUE,
#' the coordinates must be in a matrix/data frame with the longitude in the first
#' column and latitude in the second. The position is projected onto a plane in
#' meters with the function \code{\link[SoDA]{geoXY}}.
#' @examples
#' 
#' \dontrun{
#' 
#' data(eco.test)
#' eco.2spagedi(eco, "pop", ndig = 1,int=2, smax=6, name="infile.spagedi.txt")
#' 
#' }
#' 
#' @references 
#' 
#' Freedman D., and P. Diaconis. 1981. On the histogram as a density estimator: 
#' L 2 theory. Probability theory and related fields, 57: 453-476.
#' 
#' Hardy O. and X Vekemans. 2002. SPAGeDi: a versatile computer program 
#' to analyse spatial genetic structure at the individual or population levels. 
#' Molecular ecology notes, 2: 18-620.
#' 
#' Sturges  H. 1926. The choice of a class interval. Journal of the American 
#' Statistical Association, 21: 65-66.
#' 
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @export


eco.2spagedi <- function(eco, 
												 pop = NULL, 
												 ndig, 
												 name = "infile.spagedi.txt", 
												 smin = 0,
												 smax= NULL,
												 int = NULL, 
												 nclass = NULL,
												 seqvec = NULL,
												 size = NULL,
												 bin = c("sturges", "FD"),
												 distmat = NULL,
												 latlon = FALSE) {
	
	bin <- match.arg(bin)
	if(is.null(pop)) {
		pop <- rep(1, nrow(eco@XY))
	} else {
	pop <- as.numeric(eco@S[, which(colnames(eco@S)== pop)])
	}
	
	if(sum(pop == 0)) {
		stop("non matching S column name")
	}
	
  gmat <- as.matrix(eco@G)
	ploidy <- eco@INT@ploidy
  gmat[gmat == "NA" || is.na(gmat)] <- paste(rep("0", ploidy), collapse="")
  
	matriz <- data.frame(rownames(eco@P), pop, eco@XY, gmat)
	matriz <- as.matrix(matriz)
	colnames(matriz) <- c("Individual", "Population",
												colnames(eco@XY), colnames(eco@G))
	
	arriba <- rep("", ncol(matriz))
	arriba[1] <- nrow(matriz)
	arriba[2] <- max(pop)
	arriba[3] <- ncol(eco@XY)
	arriba[4] <- ncol(eco@G)
	arriba[5] <- ndig
	arriba[6] <- eco@INT@ploidy
	
	if(!is.null(smax) | !is.null(seqvec)) {
	 
  	xy <- eco@XY[,1:2]
	  if(latlon) {
	    dist(SoDA::geoXY(xy[,2], xy[,1], unit=1))
	  }
	   
  	input <- int.break(XY = xy, 
  	                   int = int, 
  	                   smin = smin,
  	                   smax = smax,
                       size = size,
  	                   nclass = nclass,
  	                   seqvec = seqvec,
  	                   latlon = latlon,
  	                   bin = bin)
	
	breakpoints <- input$breakpoints
	
	} else  {
		breakpoints <- NULL
	}
		
	final <- rep("", ncol(matriz))
	final[1] <- "END"

	
		sink(name)
		cat(arriba, sep = "\t")
  	cat("\n")
		if(!is.null(breakpoints)) {
			cat(length(breakpoints), breakpoints, sep = "\t")
			cat("\n")
		} else {
			cat(0)
		}
	  cat("\n")
		cat(colnames(matriz), sep = "\t")
	  cat("\n")
		write.table(matriz, sep = "\t", quote = FALSE, row.names = FALSE, 
								col.names = FALSE)
		cat(final, sep = "\t")
	
	  if(!is.null(distmat)) {
	  	if(class(distmat) == "dist") {
	  		distmat <- as.matrix(distmat)
	  	} 
	  	if(class(distmat) != "matrix" & class(distmat) != "data.frame") {
	  		stop("invalid distance matrix format (It should be of class: dist, matrix or data.frame")
	  	}
	  	distnames <- rownames(distmat)
	  	distmat<-data.frame(distnames, distmat)
	  	colnames(distmat)[1]<-paste("M", nrow(eco@XY), sep="")
	  	cat("\n")
	  	write.table(distmat, sep = "\t", quote = FALSE, row.names = FALSE, 
	  							col.names = TRUE)
	  	cat("END")
	  }
	 
		sink()

	 }
	
