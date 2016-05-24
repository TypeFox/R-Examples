#==============================================================================#
# diffexp-methods.R:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# calcHypergeoProb
# calcHypergeoExon
# matchExonTranscript
# diffSpliceScore
# totalMatches_cluster
# totalMatches_all
# adjust4maxprobes
#==============================================================================#


#' Calculates hypergeometric probability in DEMI analysis
#' 
#' Calculates hypergeometric probability in DEMI analysis. It is universal for all
#' DEMI analysis except for 'exon' analysis. It is used internally in DEMI analysis.
#' 
#' @param x A \code{data.frame}.
#' @return A \code{numeric} that represents a the hypergeometric probability p-value.
#' 
#' @author Sten Ilmjarv
#' 
#' @export
#' @docType methods
#' @rdname calcHypergeoProb-methods 
calcHypergeoProb <- function( x ) {
	# dhyper( x, m, n, k, log = FALSE )
	# x - the number of white balls drawn without replacement
	# m - the number of white balls in the urn
	# n - the number of black balls in the urn
	# k - the number of balls drawn from the urn
	#return(sum(dhyper(as.numeric(x[[1]]):as.numeric(x[[2]]),as.numeric(x[[2]]),(as.numeric(x[[3]])-as.numeric(x[[2]])),as.numeric(x[[4]]))))
	return( sum( dhyper( as.numeric( x["countCluster"] ):as.numeric( x["countBlat"] ),
							as.numeric( x["countBlat"] ),
							( as.numeric( x["probesInBlat"] ) - as.numeric( x["countBlat"] ) ),
							as.numeric( x["probesInCluster"] )
					) ) )
}

#' Calculates hypergeometric probability in DEMI analysis
#' 
#' Calculates hypergeometric probability in DEMI analysis. It is only used for 'exon'
#' analysis in DEMI. It is used internally in DEMI analysis.
#' 
#' @param x A \code{data.frame}.
#' @return A \code{numeric} that represents a the hypergeometric probability p-value.
#' 
#' @author Sten Ilmjarv
#' 
#' @export 
#' @docType methods
#' @rdname calcHypergeoExon-methods
calcHypergeoExon <- function( x ) {
	# dhyper( x, m, n, k, log = FALSE )
	# x - the number of white balls drawn without replacement
	# m - the number of white balls in the urn
	# n - the number of black balls in the urn
	# k - the number of balls drawn from the urn
	#return(sum(dhyper(as.numeric(x[[3]]):min(as.numeric(x[[4]]),as.numeric(x[[6]])),as.numeric(x[[6]]),(as.numeric(x[[5]])-as.numeric(x[[6]])),as.numeric(x[[4]]))))
	return( sum( dhyper( as.numeric( x["countCluster"] ):min( as.numeric( x["countBlat"] ), as.numeric( x["geneClust"] ) ),
							as.numeric( x["geneClust"] ),
							( as.numeric( x["geneTotal"] ) - as.numeric( x["geneClust"] ) ),
							as.numeric( x["countBlat"] )
					) ) )
}

#' Matches exons to their corresponding transcripts.
#' 
#' The function \code{matchExonGene} matches exons to their corresponding transcripts. It is
#' used internally in DEMI analysis.
#' 
#' @param cluster A \code{vector}. A \code{vector} of probe ID's in the cluster.
#' @param blatTable A \code{data.frame}. A \code{data.frame} with alignment information.
#' @param annoTable A \code{data.frame}. A \code{data.frame} with annotation information.
#' @return A \code{data.frame} where exons are matched to transcript.
#' 
#' @author Sten Ilmjarv
#' 
#' @export 
#' @docType methods
#' @rdname matchExonGene-methods
matchExonGene <- function( cluster, blatTable, annoTable ) {
	
	matchGene <- NULL;
	
	# Calculate matches for transcript
#	geneTable <- merge( blatTable, annoTable[ , c( "geneID", "targetID" ) ], by = "targetID" );
	geneTable <- plyr::join( blatTable, annoTable[ , c( "geneID", "targetID" ) ], by = "targetID" );
	geneTable <- geneTable[ , c( "probeID", "geneID" ) ];
	geneTable <- unique( geneTable );
	geneTable$geneID <- factor( geneTable$geneID );
	
	# Calculate total matches for transcript
	freqTable <- as.data.frame( table( geneTable$geneID ) );
	colnames( freqTable ) <- c( "geneID", "geneTotal" );
	matchGene <- freqTable;
	
	# Calculate cluster matches for transcript
	geneCluster <- geneTable[ geneTable$probeID %in% cluster, ];
	freqTable <- as.data.frame( table( geneCluster$geneID ) );
	colnames( freqTable ) <- c( "geneID", "geneClust" );
	
	# Merge two datasets together
#	matchGene <- merge( matchGene, freqTable, by = "geneID" );
	matchGene <- plyr::join( matchGene, freqTable, by = "geneID" );
	
	# Add exon back to the dataset
#	matchGene <- merge( matchGene, annoTable[ , c( "geneID", "targetID" ) ], by = "geneID" );
	matchGene <- plyr::join( matchGene, annoTable[ , c( "geneID", "targetID" ) ], by = "geneID" );
	
	# The function returns the dataset
	return( matchGene );
	
}

#' Calculate differential splice scores
#' 
#' The function \code{diffSpliceScore} calculates differential splice scores in DEMI analysis. It
#' is used internally in DEMI analysis. In the current implementation of DEMI it is not used.
#' 
#' @param x A \code{data.frame}.
#' @return Returns the differential splice score as \code{numeric}.
#' 
#' @author Sten Ilmjarv
#' 
#' @docType methods
#' @rdname diffSpliceScore-methods
diffSpliceScore <- function(x) {
	return(sqrt(log(as.numeric(x[[5]])/as.numeric(x[[4]]))) * -2 * log(as.numeric(x[[8]])))
}

#' Calculates the number of matches in the cluster
#' 
#' The function \code{totalMatches_cluster} calculates the number of matches in the cluster in DEMI
#' analysis. It is used internally in DEMI analysis.
#' 
#' @param cluster A \code{vector}. A \code{vector} of probe ID's in the cluster.
#' @param blatTable A \code{data.frame}. A \code{data.frame} with alignment information.
#' @return A \code{data.frame} that represents the number of probes for each target in the cluster.
#' 
#' @author Sten Ilmjarv
#' 
#' @export 
#' @docType methods
#' @rdname totalMatches_cluster-methods
totalMatches_cluster <- function( cluster, blatTable )
{
	
	#cat( "\t\tCalculating matches for probes in cluster ")
	cat( DEMIMessages$diffexp$matchesCluster );
	
	blatTable <- blatTable[ , c( "probeID", "targetID" )]
	
	#blatTable <- merge( blatTable, cluster, by = "probeID" )
	blatTable <- blatTable[ blatTable$probeID %in% cluster, ];
	blatTable <- unique( blatTable )
	
	freqTable <- as.data.frame( table( blatTable$targetID ) )
	colnames( freqTable ) <- c( "targetID", "countCluster" )
	
	#cat( "... Done\n" )
	cat( DEMIMessages$done() );
	
	return ( freqTable )
	
}

#' Calculates the number of matches over all probes
#' 
#' The function \code{totalMatches_all} calculates the number of matches for all probes in DEMI analysis. It
#' is used internally in DEMI analysis.
#' 
#' @param blatTable A \code{data.frame}. A \code{data.frame} with alignment information.
#' @return A \code{data.frame} that represents the number of probes for each target.
#' 
#' @author Sten Ilmjarv
#' 
#' @export
#' @docType methods
#' @rdname totalMatches_all-methods
totalMatches_all <- function( blatTable )
{
	
	#cat( "\t\tCalculating matches over all probes ")
	cat( DEMIMessages$diffexp$matchesAll );
	
	blatTable <- blatTable[ , c( "probeID", "targetID" )]
	blatTable <- unique( blatTable )
	
	freqTable <- as.data.frame( table( blatTable$targetID ) )
	colnames( freqTable ) = c( "targetID", "countBlat" )
	
	#cat( "... Done\n" )
	cat( DEMIMessages$done() );
	
	return ( freqTable )
}

#' Adjust the DEMI analysis by \code{maxprobes} analysis
#' 
#' The function \code{adjust4maxprobes} adjust the number of probes if the \code{maxprobes} has been
#' set in the \code{DEMIExperiment} object. It is used internally in DEMI analysis.
#' 
#' @param targetMatches A \code{data.frame}. The original number of probes per target stored in a \code{data.frame}.
#' @param maxprobes A \code{numeric}. Specifies the maximum number of probes a target a target is adjusted against.
#' @return A \code{data.frame} with the number of probes per target have been adjusted by \code{maxprobes}.
#' 
#' @author Sten Ilmjarv
#' 
#' @export
#' @docType methods
#' @rdname adjust4maxprobes-methods
adjust4maxprobes <- function( targetMatches, maxprobes ) {
	
	targetMatches[ targetMatches$countBlat > maxprobes, c( "countCluster") ] = as.integer( ( targetMatches[ targetMatches$countBlat > maxprobes, c( "countCluster") ] / targetMatches[ targetMatches$countBlat > maxprobes, c( "countBlat") ] ) * maxprobes )
	targetMatches[ targetMatches$countBlat > maxprobes, c( "countBlat" ) ] = maxprobes 
	
	return( targetMatches )
	
}
