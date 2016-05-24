#==============================================================================#
# DEMIDiff-methods.R:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# initialize.DEMIDiff
# getDEMIClust
# getExperiment
# getGroup
# getName
# getResult
# getResultTable
# getProbeLevel
# getTargetProbes
# demisummary
# diffexp
#==============================================================================#

#------------------------------------------------------------------------------#
# DEMIClust initialization:
#------------------------------------------------------------------------------#

#' Initializes the \code{DEMIDiff} object
#' 
#' Initializes the \code{DEMIDiff} object.
#' 
#' @param .Object A DEMIDiff object.
#' @param ... Additional arguments that may never be used.
#' @return Returns a 'DEMDiff' object.
#' @author Sten Ilmjarv
#' @import methods
"initialize.DEMIDiff" <-
		function( .Object, ... ) 
{
	.Object <- callNextMethod( .Object, ... );
	.Object <- diffexp(.Object);
	.Object;
}#initialize.DEMIClust

#' @import methods
setMethod( "initialize", "DEMIDiff", initialize.DEMIDiff );

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# DEMIExperiment get functions:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @rdname getDEMIClust-methods
#' @aliases getDEMIClust,DEMIDiff-method
#' @import methods
setMethod( "getDEMIClust", signature( object = "DEMIDiff" ),
		function( object ) object@cluster
)#getDEMIClust

#' @rdname getExperiment-methods
#' @aliases getExperiment,DEMIDiff-method
#' @import methods
setMethod( "getExperiment", signature( object = "DEMIDiff" ),
		function( object) {
			getExperiment( object@cluster )
		}
)#getExperiment

#' @rdname getGroup-methods
#' @aliases getGroup,DEMIDiff-method
#' @import methods
setMethod( "getGroup", signature( object = "DEMIDiff" ),
		function( object ) object@cluster@group
)#getGroup

#' @rdname getName-methods
#' @aliases getName,DEMIDiff-method
#' @import methods
setMethod( "getName", signature( object = "DEMIDiff" ),
		function( object ) object@name
)#getName

#' @rdname getResult-methods
#' @aliases getResult,DEMIDiff-method
#' @import methods
setMethod( "getResult", signature( object = "DEMIDiff" ),
		function( object ) object@result
)#getName

#' @rdname getResultTable-methods
#' @aliases getResultTable,DEMIDiff-method
#' @import methods
setMethod( "getResultTable", signature( object = "DEMIDiff" ),
		function( object ) makeDEMIResultsTable( list( getResult( object ) ) )	# the function takes in the list of object 'DEMIResult'
)#getResultTable

#' @rdname getProbeLevel-methods
#' @aliases getProbeLevel,DEMIDiff,vector,logical-method
#' @import methods
setMethod( "getProbeLevel", signature( object = "DEMIDiff", probes = "vector", verbose = "logical" ),
		function( object, probes, verbose ) getProbeLevel( getExperiment( object ), probes, TRUE )
)#getProbe

#' @rdname getTargetProbes-methods
#' @aliases getTargetProbes,DEMIDiff,vector-method
#' @import methods
setMethod( "getTargetProbes", signature( object = "DEMIDiff", target = "vector" ),
		function( object, target ) getTargetProbes( getExperiment( object ), target )
)#getGene

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#	DEMIDiff other functions:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @rdname demisummary-methods
#' @aliases demisummary,DEMIDiff-method
#' @import methods
setMethod( "demisummary", signature( object = "DEMIDiff" ),
		function( object, target )
		{
			#	check that the function is correcty run
			if ( missing( target ) == TRUE ) {
				#stop( paste( "parameter ", sQuote( "target" ), " is missing", sep = "" ) );
				stop( DEMIMessages$parameterMissing( "target" ) );
			} else {
				if ( is.vector( target ) == FALSE ) {
					#stop( paste( "parameter ", sQuote( "target" ), " needs to be a vector", sep = "" ) );
					stop( DEMIMessages$parameterNotOfClass( "target", "vector" ) );
				}
			}
			# 	if any targets match the criteria
			if ( isTRUE( check4target( getExperiment( object ), target ) ) ) {
				#	find the targets
				targetProbes <- droplevels( getTargetProbes( object, target ) );
				# create a data.frame from target names
				targetID <- as.vector( unique( droplevels( targetProbes )$targetID ) );
				output <- data.frame( targetID );
#				targetProbes <- as.matrix( targetProbes );
				result <- list( getResult( object ) );
				targetList <- tapply( targetProbes[, grep( "probeID", colnames( targetProbes ) )], targetProbes[, grep( "targetID", colnames( targetProbes ) )], function(x) { return( x ) } );
				#	find if groups are specified
				if ( length( result ) > 0 ) {
					# for each group check if index's exists
					for ( i in 1:length( result ) ) {
						group <- getGroup( result[[i]] );
						if ( length( getIndexA( group ) ) > 0 || length( getIndexB( group ) ) > 0 ) {
							groupA <- getGroupA( group );
							groupB <- getGroupB( group );
							output <- data.frame( output, numeric( nrow( output) ), numeric( nrow( output ) ) );
							colnames( output )[ c( ncol( output ) - 1, ncol( output ) ) ] <- c( groupA, groupB ); 
							for( ii in 1:length( targetList ) ) {
								targetName <- names( targetList[ii]);
								probeLevel <- getProbeLevel( getExperiment( object ), targetList[[ii]], verbose = FALSE );
								output[ output$targetID == targetName, grep( groupA, colnames( output ) ) ] <- mean( probeLevel[ getIndexA( group ) ] );
								output[ output$targetID == targetName, grep( groupB, colnames( output ) ) ] <- mean( probeLevel[ getIndexB( group ) ] );
							}
						}
					}
				}
				# output the whole mean of normalized matrix
				ALL <- do.call( "rbind", lapply( targetList, function( x ) {
									probeLevel <- getProbeLevel( getExperiment( object ), x, verbose = FALSE );
									return( mean( probeLevel ) );
								})
				);
				ALL <- data.frame( ALL, rownames( ALL ) );
				colnames( ALL )[ grep( "rownames.ALL.", colnames( ALL ) ) ] = "targetID";
#				output <- merge( output, ALL, by = "targetID" );
				output <- plyr::join( output, ALL, by = "targetID" );
				return( output );
			} else {
				#stop( "0 specified targets were found in the experiment" );
				stop( DEMIMessages$DEMIDiff$zeroTargetsFound );
			}
		}
);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 	DEMIDiff differential expression functions:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @rdname diffexp-methods
#' @aliases diffexp,DEMIDiff-method
#' @import methods
setMethod( "diffexp", signature( object = "DEMIDiff" ),
		function( object ) {
			
			#cat( "*Calculating differential expression " );
			cat( DEMIMessages$DEMIDiff$calcDiffExp );
			
			demiClust <- getDEMIClust( object );
			
			if ( customObject( demiClust ) == FALSE ) {
				#cat( paste( "for comparison between groups", paste( "'", demiClust@group@groupA, "'", sep = "" ), "and", paste( "'", demiClust@group@groupB, "'", sep = "" ), "\n" ) );
				cat( DEMIMessages$DEMIDiff$betweenGroups( demiClust@group@groupA, demiClust@group@groupB ) );
			} else if ( customObject( demiClust ) == TRUE ) {
				#cat( "on user-defined clusters\n" );
				cat( DEMIMessages$DEMIDiff$onUserDefined );
			}
			
			result <- list();
			
			#	define the analysis parameters
			experiment <- getExperiment( demiClust );
			analysis <- getAnalysis( experiment );
			maxprobes <- getMaxprobes( experiment );
			blatTable <- getAlignment( experiment );
			#	if some probe comes up more then once then we do not need to calculate the frequency of the target more then once
			#	that is only if the analysis is gene for gene can have several transcripts matching to it
			if ( analysis == "gene" ) {
				blatTable = blatTable[, c( "probeID", "targetID" )];
				blatTable = unique( blatTable );
			}
			blatTable$targetID = factor( blatTable$targetID );
			annoTable <- getAnnotation( experiment );
			annoTable$targetID = factor( annoTable$targetID );
			
			# Calculate median number of probes for target if 'maxprobes' is set to 'median'
			if ( is.null( maxprobes ) == FALSE && maxprobes == "median" && length( maxprobes ) != 0 ) {
				maxprobes <- median( as.vector( data.frame( table( blatTable$targetID ) )$Freq ) )
			}
			# Calculate max number of probes for target if 'maxprobes' is set to 'max'
			if ( is.null( maxprobes ) == FALSE && maxprobes == "max" && length( maxprobes ) != 0 ) {
				maxprobes <- max( as.vector( data.frame( table( blatTable$targetID ) )$Freq ) )
			}
			maxprobes <- as.numeric( maxprobes ); # convert maxprobes to numeric
			
			clusters <- getCluster( demiClust );
			clusterID <- character( length( clusters ) );
			# Do the analysis on all clusters in the 'DEMIClust' object
			if ( analysis == "transcript" || analysis == "gene" || analysis == "exon" || analysis == "genome" ) {
				for ( i in 1:length( clusters ) ) {
					
					resultsCluster <- NULL;
					
					#	DEFINE THE DATA FOR CALCULATIONS
					cluster <- clusters[[i]]; # define the clusters
					clusterName <- names( clusters[i] );
					clusterID[i] <- clusterName;
					if ( customObject( demiClust ) == FALSE ) {
						#cat( paste( "\tAnalysing cluster", paste( "'", clusterName, "'", sep = "" ), "in", paste( "'", demiClust@group@groupA, "'", sep = "" ), "compared to", paste( "'", demiClust@group@groupB, "'", sep = "" ), "\n" ) );
						cat( DEMIMessages$DEMIDiff$analyzingClusterNative( clusterName, demiClust@group@groupA, demiClust@group@groupB ) );
					} else if ( customObject( demiClust ) == TRUE ) {
						#cat( paste( "\tAnalysing cluster", paste( "'", clusterName, "'", sep = "" ), "\n" ) );
						cat( DEMIMessages$DEMIDiff$analyzingClusterCustom( clusterName ) );
					}
					
					probesInCluster <- length( cluster ); # the number of unique probes in the cluster
					probesInBlat <- length( unique( blatTable$probeID ) ); # the number of unique probes in the blat table
					totalMatches <- totalMatches_all( blatTable = blatTable ); # calculate the number of matches the target has in the blat table
					clustMatches <- totalMatches_cluster( cluster = cluster, blatTable = blatTable ); # calculate the number of matches the target has in the blat table made up only of probes in the cluster
					
					#	Merge includes all targets because clustMatches can be 0 as well if the target has no matches
					#	in the cluster - therefore everything will be included
					targetMatches <- merge( totalMatches, clustMatches, by = "targetID" );
#					targetMatches <- plyr::join( totalMatches, clustMatches, by = "targetID" ); # this gave an error for some reason
					
					# If maxprobes has been set, adjust targetMatches for maxprobes
					if ( length( ( maxprobes ) ) != 0 && maxprobes != 0 ) {
						targetMatches <- adjust4maxprobes( targetMatches = targetMatches, maxprobes = maxprobes );
					}
					
					#	CALCULATE HYPERGEOMETRIC PROBABILITY
					#cat( "\t\tCalculating hypergeometric probability p-values" );
					cat( DEMIMessages$DEMIDiff$calcHyperGeo );
					
					hypergeoPval <- NULL;
					options( warn = -1 ); # Turn off warnings
					
					if ( nrow( targetMatches ) > 0 ) {
						if ( analysis == "transcript" || analysis == "gene" || analysis == "genome" ) {
							hypergeoPval <- apply( data.frame( targetMatches, probesInBlat, probesInCluster ), 1, calcHypergeoProb );
						} else if ( analysis == "exon" ) {
							# Calculate transcript probe matches
							matchGene <- matchExonGene( cluster = cluster, blatTable = blatTable, annoTable = annoTable );
							targetMatches <- data.frame( merge( targetMatches, matchGene, by = "targetID" ) );
#							targetMatches <- data.frame( plyr::join( targetMatches, matchGene, by = "targetID" ) ); # this gave an error for some reason
							targetMatches <- unique( targetMatches[ , c( "targetID", "countCluster", "countBlat", "geneClust", "geneTotal", "geneID" ) ] );
							hypergeoPval <- apply( targetMatches, 1, calcHypergeoExon );
						}
					}
					options( warn = 1 ) # Turn warnings back on
					#cat( " ... Done\n" );
					cat( DEMIMessages$done() );
					
					#	ADD ANNOTATIONS TO THE TEST RESULTS
					
					FDR <- 1;
					#cat( "\t\tAnnotating results" );
					cat( DEMIMessages$DEMIDiff$annotatingResults );
					if ( nrow( targetMatches ) > 0 ) {
						if ( analysis == "transcript" ) {
							annoTable_temp <- annoTable[ , -c( grep( "start", colnames( annoTable ) ), grep( "length", colnames( annoTable ) ) ) ]; # remove these columns
							resultsCluster <- merge( data.frame( targetMatches[, c( "targetID", "countCluster", "countBlat" )], probesInCluster, probesInBlat, hypergeoPval, FDR ), annoTable_temp, by = "targetID", all.x = TRUE );
#							resultsCluster <- plyr::join( data.frame( targetMatches[, c( "targetID", "countCluster", "countBlat" )], probesInCluster, probesInBlat, hypergeoPval, FDR ), annoTable_temp, by = "targetID", type = "left" );
							#colnames( resultsCluster )[ grep( "targetID", colnames( resultsCluster ) ) ] = "transcriptID";
						} else if ( analysis == "gene" ) {
							annoTable_temp <- annoTable[ , -c( grep( "start", colnames( annoTable ) ), grep( "length", colnames( annoTable ) ), grep( "chr", colnames( annoTable ) ) ) ]; # remove these columns
							annoTableCluster <- annoTable_temp[, -c( grep( "targetID", colnames( annoTable_temp ) ), grep( "peptideID", colnames( annoTable_temp ) ), grep( "biotype", colnames( annoTable_temp ) ) ) ];
							colnames( annoTableCluster )[ grep( "geneID", colnames( annoTableCluster ) ) ] <- "targetID";
							annoTableCluster = unique( annoTableCluster );
							resultsCluster <- merge( data.frame( targetMatches[, c( "targetID", "countCluster", "countBlat" )], probesInCluster, probesInBlat, hypergeoPval, FDR ), annoTableCluster, by = "targetID", all.x = TRUE );
#							resultsCluster <- plyr::join( data.frame( targetMatches[, c( "targetID", "countCluster", "countBlat" )], probesInCluster, probesInBlat, hypergeoPval, FDR ), annoTableCluster, by = "targetID", type = "left" );
							#colnames( resultsCluster )[ grep( "targetID", colnames( resultsCluster ) ) ] = "geneID";
						} else if ( analysis == "exon" ) {
							resultsCluster <- unique( merge( data.frame( targetMatches, hypergeoPval, FDR ), annoTable[ , c( "targetID", "geneSymbol", "description" ) ], by = "targetID", all.x = TRUE ) );
#							resultsCluster <- unique( plyr::join( data.frame( targetMatches, hypergeoPval, FDR ), annoTable[ , c( "targetID", "geneSymbol", "description" ) ], by = "targetID", type = "left" ) );
							#colnames( resultsCluster )[ grep( "targetID", colnames( resultsCluster ) ) ] = "exonID";
						} else if ( analysis == "genome" ) {
							resultsCluster <- merge( data.frame( targetMatches[, c( "targetID", "countCluster", "countBlat" )], probesInCluster, probesInBlat, hypergeoPval, FDR ), annoTable, by = "targetID", all.x = TRUE );
#							resultsCluster <- plyr::join( data.frame( targetMatches[, c( "targetID", "countCluster", "countBlat" )], probesInCluster, probesInBlat, hypergeoPval, FDR ), annoTable, by = "targetID", type = "left" );
							#colnames( resultsCluster )[ grep( "targetID", colnames( resultsCluster ) ) ] = "sectionID";
							# add cytoband if it exists
							cytoband <- getCytoband( experiment );
							if ( nrow( cytoband ) > 0 ) {
								resultsCluster <- addCytoband( resultsCluster, cytoband );
								resultsCluster <- makeUCSCLink( resultsCluster );
							}
						}
					} # nrow( targetMatches ) > 0
					
					# Add proper column names
					colnames( resultsCluster )[ grep( "countBlat", colnames( resultsCluster ) ) ] <- "probesOnTarget";
					colnames( resultsCluster )[ grep( "countCluster", colnames( resultsCluster ) ) ] <- "probesOnTargetInCluster";
					colnames( resultsCluster )[ grep( "probesInBlat", colnames( resultsCluster ) ) ] <- "probesTotal";
					colnames( resultsCluster )[ grep( "hypergeoPval", colnames( resultsCluster ) ) ] <- "P-value";
					colnames( resultsCluster )[ grep( "geneName", colnames( resultsCluster ) ) ] <- "geneSymbol";
					colnames( resultsCluster )[ grep( "transcriptTotal", colnames( resultsCluster ) ) ] <- "probesOnTranscript";
					colnames( resultsCluster )[ grep( "transcriptClust", colnames( resultsCluster ) ) ] <- "probesOnTranscriptInCluster";
					
					result[[clusterID[i]]] <- resultsCluster;
					
					#cat( " ... Done\n" );
					cat( DEMIMessages$done() );
				}
				
				#	ADJUST FOR MULTIPLE CORRECTION
				
				#	Retrieve the table with all the clusters
				resultsAll <- NULL;
				for ( i in 1:length( result ) ) {
					resultsAll <- rbind( resultsAll, data.frame( names( result[i] ), result[[i]] ) );
				}
				colnames( resultsAll )[1] <- "clusterID";
				
				#	This situation should never arise
				if ( is.null( resultsAll$P.value ) == TRUE ) {
					resultsAll$P.value = 1;
				}
				
				FDR <- NULL;
				if ( analysis == "genome" || analysis == "transcript" ) {
					#cat( "\t\tAdjusting p-values for multiple testing by the modified FDR procedure (Benjamini & Yekutieli, 2001)" );
					cat( DEMIMessages$DEMIDiff$multipleCorrectionBY );
					FDR <- p.adjust( as.vector( resultsAll$P.value ), "BY", length( as.vector( resultsAll$P.value ) ) );
				} else {
					#cat( "\t\tAdjusting p-values for multiple testing by the FDR procedure (Benjamini & Hochberg, 1995)" );
					cat( DEMIMessages$DEMIDiff$multipleCorrectionBH );
					FDR <- p.adjust( as.vector( resultsAll$P.value ), "BH", length( as.vector( resultsAll$P.value ) ) );
				}
				resultsAll$FDR <- FDR;
				
				#	Save the results again in the list
				result <- list();
				for ( i in names( table( resultsAll$clusterID ) ) ) {
					result[[i]] <- resultsAll[ resultsAll$clusterID == i, -c( grep( "clusterID", colnames( resultsAll ) ) ) ]
				}
				
				#cat( " ... Done\n" );
				cat( DEMIMessages$done() );
				
			}
			
			#names( result ) <- clusterID;
			
			demiResult <- new( "DEMIResult", group = getGroup( demiClust ), result = result );
			
			#cat( "*Calculating differential expression...Done\n" );
			cat( DEMIMessages$DEMIDiff$diffExpDone );
			
			object@result <- demiResult;
			
			return( object );
		}

)#diffexp
