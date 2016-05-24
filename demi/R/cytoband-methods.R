#==============================================================================#
# cytoband-methods.R:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# addCytoband
# findCytoband
# makeUCSCLink
#==============================================================================#

#' Add karyotype information to DEMI differential expression results
#' 
#' The function \code{addCytoband} adds karyotype information to the results of DEMI differential
#' expression 'genome' analysis. It is used internally in DEMI analysis.
#' 
#' @param result A \code{data.frame}. A 'genome' analysis genome section that contains chromosome name,
#' 		  region start and region end coordinates.
#' @param cyto A \code{data.frame}. A \code{data.frame} describing karyotype information of the organism
#' 		  used in the analysis.
#' @return A \code{data.frame} where karyotype information has been added to the input 'result' table.  
#' @seealso \code{\link{findCytoband}} which this function wraps
#' 
#' @author Sten Ilmjarv
#' 
#' @export 
#' @docType methods
#' @rdname addCytoband-methods
"addCytoband" <- function( result, cyto )
{
	cytoband <- apply( result[, c( "chr", "start", "end" ) ], 1, findCytoband, cyto );
	result <- data.frame( result, cytoband );
	return( result );
}# addCytoband

#' Finds cytoband for the specified genome region
#' 
#' The function \code{fincCytoband} finds cytoband for a genome region specified by the chromosome,
#' region start and region end coordinates. It is used internally in DEMI analysis.
#' 
#' @param x A \code{vector}. A vector of "chr", "start" and "end" information about the genome region.
#' @param cytoband A \code{data.frame}. A data.frame containing karyotype information.
#' @return A karyotype \code{character} corresponding to the input genomic region.
#' 
#' @author Sten Ilmjarv
#'  
#' @export 
#' @docType methods
#' @rdname findCytoband-methods
"findCytoband" <- function( x, cytoband = "data.frame" )
{
	start <- as.character( cytoband[ cytoband$chr == x["chr"] & as.numeric( cytoband$start ) <= as.numeric( x["start"] ) & as.numeric( cytoband$end ) > as.numeric( x["start"] ), c( "region" ) ] );
	end <- as.character( cytoband[ cytoband$chr == x["chr"] & as.numeric( cytoband$start ) < as.numeric( x["end"] ) & as.numeric( cytoband$end ) >= as.numeric( x["end"] ), c( "region" ) ] );
	ret <- paste( x["chr"], start, "-", end, sep = "" );
	return( ret );
}#findCytoband

#' Make UCSC link
#' 
#' The function \code{makeUCSCLink} makes a UCSC link of every genomic region in the
#' specified \code{data.frame}. It is used internally in DEMI analysis.
#' 
#' @param result A \code{data.frame}. A \code{data.frame} that consists of chromosome name and start
#' 		  and end coordinates to be used to make the UCSC link.
#' @return The input \code{data.frame} with added UCSC link as the last column.
#' 
#' @author Sten Ilmjarv
#' 
#' @export 
#' @docType methods
#' @rdname makeUCSCLink-methods
"makeUCSCLink" <- function( result )
{
	ucsclink <- paste( result$chr, ":", result$start, "-", result$end, sep = "" );
	result <- data.frame( result, ucsclink );
	result <- result[ ,-c( grep( "chr", colnames( result ) ), grep( "start", colnames( result ) ), grep( "end", colnames( result ) ) )]
	return( result );
}#makeUCSCLink
