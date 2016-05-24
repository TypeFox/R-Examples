#' A function to identify neighboring cells.
#' 
#' Returns cell numbers of neighboring pixels to a reference cell
#'
#' @param cell cell number from which to compute neighbor cells
#'
#' @return returns cell numbers of neighboring pixels to the reference one
#'
#' @keywords neighborhood
#'
#' @export
#' 
#' @examples
#' ## Not run neigh_cell(n)

 neigh_cell<-function(cell) {

	nr<-get('nr')
	nc<-get('nc')

	c<-ceiling(cell/nr); r<-(cell+nr-c*nr)
	v.c<-c(c,c+1,c+1,c+1,c,c-1,c-1,c-1)
	v.r<-c(r+1,r+1,r,r-1,r-1,r-1,r,r+1)
	wc.neg<-which(v.c < 1 | v.c > nc)
	if (length(wc.neg) > 0) {
		 v.c<- v.c[-wc.neg] 
		v.r<-v.r[-wc.neg]
	}
	wr.neg<-which(v.r < 1 | v.r > nr)
	if (length(wr.neg) > 0) {
		 v.c<-v.c[-wr.neg]
		 v.r<-v.r[-wr.neg] }
	cell.n<-function(r,c) return(c*nr-nr+r) 
	return(cell.n(r=v.r,c=v.c))
	}
