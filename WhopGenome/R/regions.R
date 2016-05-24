###
###
###		Regions
###
###
###

#
#
#
vcf_getregion <- function( vcffh )	return( .Call("VCF_getRegion", vcffh, PACKAGE="WhopGenome" ) )

#
#	set region to scan - each component is specified individually or as "chr:begin-end" string
#
vcf_setregion <- function( vcffh , tid, from=NA, to=NA )
{
	if( is.na(from) || is.na(to) )
	{
		tid_reg = strsplit( tid, ":" )
		stopifnot( length(tid_reg[[1]]) == 2 )
		tid = tid_reg[[1]][1]
		beg_end = strsplit( tid_reg[[1]][2] , "-" )[[1]]
		from = as.integer( beg_end[1] )
		to = as.integer( beg_end[2] )
	}
	tid = as.character( tid )
	from = as.integer(from)
	to = as.integer(to)
	if( from <= 0 || to <= 0 || to < from ){	stop("begin-end is wrong!") }
	res <- .Call("VCF_setRegion", vcffh, tid, from, to, PACKAGE="WhopGenome" )
	return( res );
}

#
#	scan last set region from start again
#
vcf_restartregion <- function( vcffh )	.Call("VCF_restartRegion",vcffh,PACKAGE="WhopGenome")

