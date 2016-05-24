
############################################
# BIFIEcdata format for replicate weights
cdata.wgtrep <- function( wgtrep ){
	N <- nrow(wgtrep)
	RR <- ncol(wgtrep)
	longwgtrep <- matrix( wgtrep , nrow=N*RR , ncol=1 )
	unique_wgt <- sort( unique( longwgtrep ) )
	indexwgtrep <- match( longwgtrep , unique_wgt )
	indexwgtrep <- matrix( indexwgtrep , nrow=N , ncol=RR )
	# list with compactly saved replicate weights
	wgtrep_list <- list( "unique_wgt" = unique_wgt , "indexwgtrep" = indexwgtrep )
	return(wgtrep_list)
		}