#' Generalized Procrustes Analysis 
#'
#' This function aligns a set of landmark configurations using Generalized Procrustes Analysis (GPA).   
#' @param dat a list containing landmark coordinate data of anchor from the specimens of interest
#' @param anchor.index a numeric constant for the anchor of interest; 1 for ventral right; 
#' 2 for ventral left; 3 for dorsal right; 4 for dorsal left
#' @param x a list providing the indices of specimens with anti-clockwise and clockwise orientation of
#' landmarks for the anchor with index \code{anchor.index}
#' @param PrinAxes logical; controls the argument with the same name in \code{gpagen}
#' @param showplot logical; if TRUE, a scatter plot of the GPA coordinates of all specimens of interest is returned
#' @return a list where the components are arrays of GPA coordinates for the specimens of interest
#' @details This function is essentially a wrapper for the \code{gpagen} function 
#' in the \code{geomorph} package (version 3.0.0) to help convert 
#' raw landmark coordinates of monogenean anchors to GPA coordinates for downstream analysis.
#' @seealso \code{\link{matrix2list}}, \code{\link{anglePolygon}}
#' @author Tsung Fei Khang \email{tfkhang@@um.edu.my}
#' @references Khang TF, Soo OYM, Tan WB, Lim LHS. (2016). Monogenean anchor morphometry: systematic value, phylogenetic signal, and evolution. PeerJ 4:e1668.
#'
#' Adams DC, Otarola-Castillo E. (2013). geomorph: an R package for the collection and analysis of geometric
#' morphometric shape data. Methods in Ecology and Evolution 4:393-399.
#' @examples
#' library(geomorph)
#'
#' data(ligophorus_tpsdata)
#'
#' #A data processing step to parse out the orientation of landmarks 
#' #from samples of L.parvicopulatrix
#'
#' O <- matrix(0, length(ligophorus_tpsdata$parvicopulatrix), 4)
#' for(w in 1:length(ligophorus_tpsdata$parvicopulatrix)){
#'	result <- mapply(function(k)
#'	anglePolygon(matrix2list(ligophorus_tpsdata$parvicopulatrix[[w]][(11*(k-1)+1):(11*k),]),
#'	degree=TRUE), k=1:4)
#'
#'	result_angle <- mapply(function(k) list(result[[2*k-1]]), k=1:4)
#'	result_orientation <- mapply(function(k) list(result[[2*k]]), k=1:4)
#'	names(result_angle) <- names(result_orientation) <- c("VR","VL","DR","DL")
#'	O[w,] <- unlist(result_orientation)
#' }
#'
#' mdir <- apply(O, 2, function(k) which(k == "m") )
#' pdir <- apply(O, 2, function(k) which(k == "p") )
#'
#' e <- 1 #Ventral right anchor
#' result <- procrustesFit(ligophorus_tpsdata$parvicopulatrix, e, 
#' list(mdir[[e]], pdir[[e]]), PrinAxes=TRUE, showplot=TRUE)
#'
#' #Standardize the x-coordinate of Landmark 7 by rotating the x-coordinate 
#' #of its mean GPA xy-coordinate to x=0.
#' coordinates <- stdLM(result$coords, reflect=FALSE, swap=TRUE, sgn=c(1,-1))
#'
#' plotLM(coordinates, "VR", pointscale=0.8,axispointscale=0.8,
#' meansize=1.2,polygon.outline=TRUE,c(-.6,.6),c(-.6,.6) )
#'

procrustesFit <- function (dat, anchor.index, x, PrinAxes = FALSE, showplot = FALSE) 
{
    l1 <- length(x[[1]])
    l2 <- length(x[[2]])
    n <- length(dat)
    if (l1 > 0) {
        test <- vector("list", l1)
        for (i in 1:l1) {
            test[[i]] <- dat[[x[[1]][i]]][(11 * (anchor.index - 
                1) + 1):(11 * anchor.index), ]
        }
        V1 <- array(, dim = c(11, 2, l1))
        for (i in 1:l1) {
            V1[, , i] <- test[[i]]
        }
    }
    if (l2 > 0) {
        test <- vector("list", l2)
        for (i in 1:l2) {
            test[[i]] <- dat[[x[[2]][i]]][(11 * (anchor.index - 
                1) + 1):(11 * anchor.index), ]
        }
        V2 <- array(, dim = c(11, 2, l2))
        for (i in 1:l2) {
            V2[, , i] <- test[[i]]
        }
    }
    if ((anchor.index == 2 | anchor.index == 4) & l1 > 0) {
        for (i in 1:l1) {
            V1[, , i][, 1] <- V1[, , i][, 1] * -1
        }
    }
    else if ((anchor.index == 1 | anchor.index == 3) & l2 > 0) {
        for (i in 1:l2) {
            V2[, , i][, 1] <- V2[, , i][, 1] * -1
        }
    }
    if (l1 > 0 & l2 > 0) {
        T <- array(, dim = c(11, 2, n))
        for (i in 1:l1) {
            T[, , i] <- V1[, , i]
        }
        for (i in 1:l2) {
            T[, , (i + l1)] <- V2[, , i]
        }
    }
    else if (l1 == 0) {
        T <- array(, dim = c(11, 2, n))
        for (i in 1:l2) {
            T[, , (i + l1)] <- V2[, , i]
        }
    }
    else if (l2 == 0) {
        T <- array(, dim = c(11, 2, n))
        for (i in 1:l1) {
            T[, , i] <- V1[, , i]
        }
    }
    if (showplot == TRUE) {
        superimpose <- gpagen(T, ProcD = TRUE, 
            PrinAxes = PrinAxes)
        plotAllSpecimens(superimpose$coord)
        xcordmat <- mapply(function(k) superimpose$coord[, , 
            k][, 1], k = 1:n)
        ycordmat <- mapply(function(k) superimpose$coord[, , 
            k][, 2], k = 1:n)
        xcordmat <- t(xcordmat)
        ycordmat <- t(ycordmat)
        bigcircle <- cbind(apply(xcordmat, 2, mean), apply(ycordmat, 
            2, mean))
        polygon(bigcircle[, 1], bigcircle[, 2], lwd = 2, border = "black")
    }
    else if (showplot == FALSE) {
        superimpose <- gpagen(T, ProcD = TRUE, 
            PrinAxes = PrinAxes)
    }
    return(superimpose)
}

