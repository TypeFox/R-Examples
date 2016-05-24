################################################################################
# AUXILIARY FUNCTIONS
################################################################################

#------------------------------------------------------------------------------#
# Union of classes

setClassUnion("dataframeORmatrix",
							c("data.frame", "matrix"))
setClassUnion("characterORnull",
							c("character", "NULL"))
setClassUnion("characterORmissing",
              c("character", "missing"))
setClassUnion("listORnull", 
              c("list","NULL"))
setClassUnion("factorORnull",
              c("factor","NULL"))
setClassUnion("callORnull",
              c("call","NULL"))
setClassUnion("intORnumeric", 
              c("integer","numeric","NULL"))
setClassUnion("intORmissing", 
              c("integer","missing","NULL"))
setClassUnion("intORnull", 
              c("integer","NULL"))


#------------------------------------------------------------------------------#
#' Scaling a data frame or matrix to 0 - 1 range
#' 
#' @param dfm Data frame, matrix or vector to scale.
#' @description This program scales each column of a data frame or a matrix 
#' to 0-1 range, computing ((X)\emph{ij} - (Xmin)\emph{i}) / range((X)\emph{i}) 
#' for each individual \emph{j} of the variable \emph{i}.
#' @examples
#' 
#' \dontrun{
#' 
#' data(eco.test)
#' require(adegenet)
#' pc <- dudi.pca(eco[["P"]], scannf = FALSE, nf = 3)
#' pc <- pc$li
#' pc <- aue.rescale(pc)
#' plot(eco[["XY"]][, 1], eco[["XY"]][, 2], col = rgb(pc), pch = 16, 
#' cex = 1.5, xlab ="X", ylab= "Y")
#' 
#' }
#' 
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar} 
#' @export

aue.rescale  <- function(dfm) {
  dfm <- as.data.frame(dfm)
  col <- apply(dfm, 2, function(X) { 
    (X - min(X, na.rm = TRUE)) / diff(range(X, na.rm = TRUE))
  })
  return(col)
}


#------------------------------------------------------------------------------#
#' Detection of metacharacters
#' @param X character string
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar} 
#' @keywords internal

is.meta <- function(X) {
  meta <- c("\\.", "\\\\", "\\|", "\\(", "\\)", "\\[", "\\]", "\\{", "\\}", 
            "\\(", "\\)", "\\^", "\\*", "\\?", "\\+", "\\$")
  any(meta %in% paste("\\", X, sep = ""))
}


#------------------------------------------------------------------------------#
#' Metachacter to character
#' @param X character string
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @keywords internal

meta2char <- function(X) {
  if(!is.character(X)) {
  X <- deparse(substitute(X))
  }
  if(is.meta(X)) {
    X <- paste("\\", X, sep = "")
  }
  X
}


#------------------------------------------------------------------------------#
#' Ordering the content of cells in a matrix. Ordering alleles in a genetic matrix.
#' @description This program takes a matrix and orders
#' the content of each cell. It was specially designed 
#' for genetic data, but can be used with any data 
#' that can be rearrenged by the function \code{\link{order}}.
#' The arguments ploidy and ncode determine the mode of
#' ordering the data. 
#' The cells corresponding to each individual \emph{i} and 
#' loci \emph{j} are ordered in ascending order in default option
#' (it can be passed decreasing = TRUE as argument, if descending order is desired). 
#' For example, a locus with ploidy = 2 and ncod =1,  coded as 51, 
#' for an individual, will be recoded as 15. A locus with ploidy = 3
#' and coded as 143645453, will be recoded as 143453645 (alleles 143, 454 and 645).
#' 
#' 
#' @param X Any matrix with content to order.
#' @param ncod Number of digits coding each allele
#'  (e.g., 1: x, 2: xx, 3: xxx, etc.). If NULL, ncode will we 
#'  obtained from the ploidy and the maximum number of characters
#'  in the data cells.
#' @param ploidy Ploidy of the data.
#' @param sep.loc Character string separating alleles.
#' @param chk.plocod  Defalult TRUE. The function checks coherence 
#' in ploidy and number of digits coding alleles.
#' @param ... Additional arguments passed to \code{\link{order}}
#' @examples
#' 
#' \dontrun{
#' 
#' # Example 1----------------------
#' 
#' geno <- c(12, 52, 62, 45, 54, 21)
#' geno <- matrix(geno, 3, 2)
#' 
#' # ordering the data
#' aue.sort(geno, ploidy = 2)
#' 
#' # decreasing sort order
#' aue.sort(geno, ploidy = 2, decreasing = TRUE)
#' 
#' 
#' # Example 2----------------------
#' 
#' geno2 <- c(123456, 524556, 629359, 459459, 543950, 219405)
#' geno2 <- matrix(geno2, 3, 2)
#' 
#' # ordering the data as diploid
#' aue.sort(geno2, ploidy = 2)  # the data is ordered using blocks of 3 characters
#' 
#' # ordering the data as triploid
#' aue.sort(geno2, ploidy = 3)  # the data is ordered using blocks of 2 characters
#' 
#' # error: the ploidy and the number of characters are not congruent
#' aue.sort(geno2, ploidy = 5) 
#' 
#' # error: the ploidy and the number of characters are not congruent
#' aue.sort(geno2, ploidy = 5)
#' 
#' 
#' # Example 3----------------------
#' 
#' # any character data
#' generic <- c("aldk", "kdbf", "ndnd", "ndkd")
#' generic <- matrix(generic, 2, 2)
#' aue.sort(generic, ploidy = 2) 
#' aue.sort(generic, ploidy = 4)
#' 
#' }
#' 
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @export

aue.sort <- function(X, ncod = NULL, ploidy = 2, 
                     sep.loc = "", chk.plocod = TRUE, ...)  {
  
  
  if(ploidy %% 1 != 0) {
    stop("ploidy argument must be non-fractional")
  }
    
  X <- as.matrix(X)
  nind <- nrow(X)
  nloc <- ncol(X)
  ploidy
  
  #if(check) {
  #X <- int.check.colnames(X)
  #X <- int.check.rownames(X)
  #}
  
  lista <- int.loc2listal(X, ncod = ncod, ploidy = ploidy,
                          sep.in = sep.loc, chk.plocod = chk.plocod,
                          chk.names = FALSE)
  
  # creating a list with individual locus 
  # and a list with ordered alleled by locus
 
  lista.orden <- lapply(lista, function(x) t(apply(x, 1, function(u) order(u, ...))))

  # ordering alleles
  for(i in 1:ncol(X)) {
    for(j in 1:nrow(X)) {
      ordlist <- lista.orden[[i]][j, ]
      lista[[i]][j, ] <- lista[[i]][j, ordlist]
    }
    lista[[i]] <- apply(lista[[i]], 1, function(x) paste(x, collapse = sep.loc))
  }
  
  # creating the output
  out <- do.call(cbind, lista)
  
  #replacing multiples "NA" (NANANA..) by NA
  out <- gsub("(NA)+", NA, out)
  
  rownames(out) <- rownames(X)
  colnames(out) <- colnames(X)
  
  out
  
  }


#------------------------------------------------------------------------------#
#' Identification of polymorphic loci
#' @param X allelic frequencies matrix
#' @param poly.level polymorphism threshold
#' @keywords internal

aue.is.poly <- function(X, poly.level) {
  
  if(poly.level < 0 || poly.level > 100) {
    stop("poly.level must be a number between 0 and 100")
  }
  
  temp <- (100 * apply(X, 2, mean, na.rm = TRUE)) 
  (temp >= poly.level) & (100 - temp >= poly.level) 
  
}
   

#------------------------------------------------------------------------------#
#' Remotion of non polymorphic loci
#' @param X allelic frequencies matrix
#' @param poly.level polymorphism threshold
#' @keywords internal

aue.rm.nonpoly <- function(X, poly.level) {
  
  isPoly <- aue.is.poly(X, poly.level)
  out <- X[,  isPoly,  drop = FALSE]
  out
}


#------------------------------------------------------------------------------#
#' Creation of a sequence of numbers in matrix or list format, for indexing.
#' @param from start of sequence
#' @param to end of sequence
#' @param by interval between elements
#' @param out.output format ("matrix2 or "list")
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar} 
#' @keywords internal

aue.seqlist <- function(from, to, by, out.format = c("matrix", "list")) {
  
  out.format <- match.arg(out.format)
  
  temp <- list()
  j <- 1
  k <- from
  for(i in seq(from = from + by - 1, to = to, by = by)) {
    temp[[j]] <- k:i
    j <- j+1
    k <- i+1
  }
  if(out.format == "matrix") {
    temp <- sapply(temp, function(x) x)
    if(is.null(nrow(temp))){
      temp <- matrix(temp, ncol = length(temp))
    }
  }
  
  temp
}


#------------------------------------------------------------------------------#
#' Remove spaces and tabs at the begining and the end of each element of charvec
#' @param charvec character vector
#' @author Thibaut Jombart
#' @keywords internal

aue.rmspaces <- function(charvec){
  charvec <- gsub("^([[:blank:]]*)([[:space:]]*)","",charvec)
  charvec <- gsub("([[:blank:]]*)([[:space:]]*)$","",charvec)
  return(charvec)
}


#------------------------------------------------------------------------------#
#' Generation of generic labels of constant length
#' @param base a character string
#' @param n number of labels
#' @author Thibaut Jombart
#' @keywords internal

aue.genlab <- function(base, n) {
  f1 <- function(cha, n){
    if(nchar(cha) < n){
      cha <- paste("0", cha, sep="")
      return(f1(cha, n))
    } else {return(cha)}
  }
  w <- as.character(1:n)
  max0 <- max(nchar(w))
  w <- sapply(w, function(cha) f1(cha, max0))
  return(paste(base, w, sep = ""))
}


#------------------------------------------------------------------------------#
#' Allelic frequencies 
#' @param eco Object of class ecogen.
#' @param grp Column in the slot S for summary by group.
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @keywords internal

aue.fqal <- function(eco, grp = NULL) {
  if(!is.null(grp)) {
    cual <- which(colnames(eco@S) == grp)
    grp.num <- as.numeric(levels(eco@S[,cual]))[eco@S[,cual]] 
    nfact <- max(grp.num)
  } else {
    dummy <- rep(1, nrow(eco@G))
    eco@S <- data.frame(dummy)
    cual <- which(colnames(eco@S) =="dummy")
    grp.num <- as.numeric(levels(as.factor(eco@S[,cual])))[eco@S[,cual]] 
    nfact <- 1
  }
  
  out <- list()
  for(i in 1:nfact) {
    eco2 <- eco[which(eco@S[,cual] == i)]
    clases<- as.numeric(eco2@INT@loc.fac)
    tabla <- eco2@A
    tabla <- 2 * tabla
    frecuencia <- apply(tabla, 2, sum)
    alelos.locus <- tapply(frecuencia, clases, sum)
    for( j in 1:length(clases)) {
      temp <- clases[j]
      clases[j] <- alelos.locus[temp]
    }
    frecuencia <- frecuencia / clases
  }
  
  round(frecuencia, 4)
}

#' EcoGenetics slot standard notation. Returns an accessor to the slot <ecoslot>
#' of the object <X>
#' @param ecoslot Slot of EcoGenetics object
#' @param X object name
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @keywords internal

aue.access <- function(ecoslot, X) {
  
  if(class(X) == "name") {
    class(X) <- "character"
  }
  
  if(class(ecoslot) == "name") {
    class(ecoslot) <- "character"
  }
  
  if(!is.character(X)) {
    X <- deparse(substitute(X))
  }
  if(!is.character(ecoslot)) {
    ecoslot <- deparse(substitute(ecoslot))
  }
  paste("ecoslot.",ecoslot,"(",X,")", sep = "")
}

#' .printaccess
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @keywords internal
#' 
# slot access message

.printaccess <- function() {
  cat("--------------------------------------------------------------------------\n")
  cat(" Access to slots:",
      "<ecoslot.> + <name of the slot> + <(name of the object)>","\n",
      "See help(\"EcoGenetics accessors\")")
}


################################################################################
## GRAPHICAL FUNCTONS-----------------------------------------------------------
################################################################################


### Graphical elements

#------------------------------------------------------------------------------#
#' Circle perimeter
#' @param mat Input raster matrix.
#' @param r Radius of the circle in pixels.
#' @param x0 Circle center- x position.
#' @param y0 Circle center- y position.
#' @param smooth Smoothing factor.
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @keywords internal

aue.circle <- function(mat, r, x0, y0, smooth = 100) { 
	
	for(theta in seq(0, 2*pi, pi/smooth)) {
		xd <- round(r*cos(theta)) 
		yd <- round(r*sin(theta)) 
		mat[x0 + xd, y0 + yd] <- 1
	}
	mat
}


#------------------------------------------------------------------------------#
#' Solid circle
#' @param mat Input raster matrix.
#' @param r Radius of the circle in pixels.
#' @param x0 Circle center- x position.
#' @param y0 Circle center- y position. 
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar} 
#' @keywords internal

aue.point <- function(mat, r, x0, y0) { 
	xmat <- col(mat)
	ymat <- row(mat)
	mat2 <- mat - mat
	mat2[sqrt((xmat - x0)^2 + (ymat -y0)^2) <= r] <- 1
	mat2
}


#------------------------------------------------------------------------------#
#' Solid diamond
#' @param mat Input raster matrix.
#' @param r Radius of the square in pixels.
#' @param x0 Square center- x position.
#' @param y0 Square center- y position.
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @keywords internal

aue.diamond<- function(mat, r, x0, y0) {
	xmat <- col(mat)
	ymat <- row(mat)
	mat2 <- mat - mat
	mat2[abs(xmat-x0) + abs(ymat -y0) <= r] <- 1
	mat2
}


#------------------------------------------------------------------------------#
#' Solid square
#' @param mat Input raster matrix.
#' @param r Radius of the square in pixels.
#' @param x0 Square center in the x direction.
#' @param y0 Square center in the y direction.
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @keywords internal
 
aue.square <- function(mat, d, x0, y0) {
	xmat <- col(mat)
	ymat <- row(mat)
	mat2 <- mat - mat
	mat2[pmax(abs(xmat-x0), abs(ymat -y0)) <= d] <- 1
	mat2
}


### Other image manipulation functions
#------------------------------------------------------------------------------#
#' Local filter
#' @description This program applies a function defined by the user, 
#' into a circle of radius r around each pixel of a raster, assigning
#' the result to the focal pixel.
#' @param mat Input raster matrix.
#' @param r Radius of the filter in pixels.
#' @param fun Function to apply.
#' @examples
#' 
#' \dontrun{
#' 
#' ras <- matrix(eco[["P"]][,1],15,15)
#' image(ras)
#' ras.mean <- aue.filter(ras, 3, mean)
#' image(ras.mean)
#' 
#' }
#' 
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @keywords internal
#' @export

aue.filter <- function(mat, r, fun) {
  mresamp <- mat - mat
  for(i in 1:nrow(mat)) {
    for(j in 1:ncol(mat)) {
      area <- aue.point(mresamp, r, i, j)
      tot <- sum(area != 0)
      mresamp[i, j] <- fun(mat * area)
    }
  }
  t(mresamp)
}


#------------------------------------------------------------------------------#
#' Radial distance to a point.
#' @param mat Input raster matrix.
#' @param x0 Circle center- x position.
#' @param y0 Circle center- y position.
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @keywords internal

aue.circle.w <- function(mat, x0, y0) { 
  xmat <- col(mat)
  ymat <- row(mat)
  mat2 <- mat - mat
  mat2 <- sqrt((xmat - x0)^2 + (ymat -y0)^2)
  mat2
}


#------------------------------------------------------------------------------#
#' Transforming a raster into a data frame with cartesian coordinates
#' 
#' @description This function returns a data frame with the column number (x),
#' row number (y) and cell value (z) of each pixel in a raster.
#' @param mat Input raster matrix.
#' 
#' @examples
#' 
#' \dontrun{
#' 
#' ras <- matrix(eco[["P"]][,1],15,15)
#' image(ras)
#' ras.row <- aue.image2df(ras)
#' ras.row
#' image(matrix(ras.row[,3], 15, 15))
#' 
#' }
#' 
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @export

aue.image2df <- function(mat) {
	xdim <- nrow(mat)
	ydim <- ncol(mat)
	mat.twocol <- data.frame(matrix(0, xdim * ydim, 3))
	colnames(mat.twocol) <- c("x", "y", "z")
	count <- 1
	for(j in 1:ydim) {
		for(i in 1:xdim) {
			mat.twocol[count, ] <- c(i, j, mat[i, j])
			count <- count + 1
		}
	}
	mat.twocol
}

#------------------------------------------------------------------------------#

