`is.gap` <-
function(x, gap.char=c("-",".")) {
	if(is.pdbs(x) || class(x)=="fasta") {
		return( colSums( matrix( as.logical(is.na(x$ali)+(x$ali %in% gap.char)), 
			ncol=ncol(x$ali)) ) > 0 )
	} else { 
	  return( as.logical( is.na(x) + (x %in% gap.char) ) )
	}
}
