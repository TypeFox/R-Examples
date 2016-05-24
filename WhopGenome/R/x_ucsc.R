###
###
###		UCSC Genome Browser SQL Server query
###
###
###
###



##
##
##	typing-shortcut : the arguments are concatenated into a string and sent as SQL query to the UCSC servers
##
whop.ucsc.query <- function( ... )
{
	return( DBI::dbGetQuery( .WHOPGENenv$ucsc[["connection"]] , .WHOPGENenv$pp(...) ) );
}



##
##
##	e.g.  whop.ucsc.genesForRegion(2,771680000,72000000)
##
whop.ucsc.genesForRegion <- function( chrom, beg, end )
{
	querystring = .WHOPGENenv$pp(	"select * from hg19.refFlat where chrom = 'chr",chrom,
					"' AND txStart > ",as.integer(beg)," AND txEnd < ",as.integer(end)," LIMIT ",.WHOPGENenv$ucsc[["limit"]],";"
				)
	print ( querystring );
	res <- whop.ucsc.query( querystring	)
	return( res )
}


##
##	Get Gene Info by Name
##
##
whop.ucsc.geneInfo <- function(gen,chr=NA) if( is.na(chr) ) whop.ucsc.query( .WHOPGENenv$pp("select * from hg19.refFlat where geneName = '",gen," LIMIT ",.WHOPGENenv$ucsc[["limit"]],";" ) ) else whop.ucsc.query( .WHOPGENenv$pp("select * from hg19.refFlat where chrom = '",chr,"' AND geneName = '",gen,"' LIMIT ",.WHOPGENenv$ucsc[["limit"]],";" ) );



##
##
##
##
whop.ucsc.geneInfoSimilar <- function(gen,chr=NA) if( is.na(chr) ) whop.ucsc.query( .WHOPGENenv$pp("select * from hg19.refFlat where geneName LIKE '%",gen,"%' LIMIT ",.WHOPGENenv$ucsc[["limit"]],";" ) ) else whop.ucsc.query( .WHOPGENenv$pp("select * from hg19.refFlat where chrom = '",chr,"' AND geneName LIKE '%",gen,"%' LIMIT ",.WHOPGENenv$ucsc[["limit"]],";" ) );
