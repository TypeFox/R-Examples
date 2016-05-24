#' Partitioning a matrix by row into objects of a list
#'
#' This function partitions each row of a matrix into objects of a list. 
#' @param x a matrix 
#' @return a list of objects; the number of objects is equal to the number of rows of \code{x}; 
#' each object is a vector of length equal to the number of columns of \code{x}
#' @details This function converts the landmark coordinate data into a format suitable as input for 
#' \code{anglePolygon}.
#' @seealso \code{\link{anglePolygon}}
#' @author Tsung Fei Khang \email{tfkhang@@um.edu.my}
#' @references Khang TF, Soo OYM, Tan WB, Lim LHS. (2016). Monogenean anchor morphometry: systematic value, phylogenetic signal, and evolution. PeerJ 4:e1668.
#' @examples
#' data(ligophorus_tpsdata)
#' 
#' #Check right ventral anchor polygon of the first specimen
#' anglePolygon(matrix2list(ligophorus_tpsdata$bantingensis[[1]][1:11,]), degree=TRUE)
#'
#' #Now check the rest
#' anglePolygon(matrix2list(ligophorus_tpsdata$bantingensis[[1]][12:22,]), degree=TRUE)
#' anglePolygon(matrix2list(ligophorus_tpsdata$bantingensis[[1]][23:33,]), degree=TRUE)
#' anglePolygon(matrix2list(ligophorus_tpsdata$bantingensis[[1]][34:44,]), degree=TRUE)
#'
#' #A more efficient way of doing things
#' result <- mapply(function(k) {
#' anglePolygon(matrix2list(ligophorus_tpsdata$bantingensis[[1]][(11*(k-1)+1):(11*k),]),
#' degree=TRUE)}, k=1:4)
#'
#' result_angle <- mapply(function(k) list(result[[2*k-1]]), k=1:4)
#' result_orientation <- mapply(function(k) list(result[[2*k]]), k=1:4)
#' names(result_angle) <- names(result_orientation) <- c("VR","VL","DR","DL")
#'

matrix2list <- function(x) {
	y <- mapply(function(k) list(x[k,]), k=1:nrow(x))
	return(y)
	}