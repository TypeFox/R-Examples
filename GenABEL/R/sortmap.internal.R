#' Internal function for map-sorting
#' 
#' Internal function for map-sorting, not supposed 
#' to be used directly by user (is open for regression 
#' testing reasons)
#' 
#' @param chrom vector of markers' chromosomes
#' @param map vector of marlers' map ositions
#' @param delta step to do between chroms when building cumulative map
#' 
#' @return list, withe elements 'ix' ('sorted' order), etc.
#' 
#' @author Yurii Aulchenko
#'
"sortmap.internal" <-
		function(chrom,map,delta=1) {
	chnum <- chrom.char2num(chrom)
	
	ix <- order(chnum,map)
	
	map <- map[ix]
	
	off <- c(0,map[1:(length(map)-1)])
	off <- map - off
	off[which(off<=0)] <- delta
	cummap <- cumsum(off)
	
	out <- list()
	out$ix <- ix
	out$cummap <- cummap
	out$chnum <- chnum
	out
}