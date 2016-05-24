#' Add new alignments to the alignment table
#' 
#' The function \code{addCustomTargets} adds additional rows to the alignment table. Before adding
#' new alignments the user needs to verify that the format of the rows corresponds to the format
#' of the alignment table. These new alignments will then be incorporated in the analysis.
#' 
#' @param object A \code{DEMIExperiment} object. Determines the experiment where the additional alignment
#' 		  information will be added to.
#' @param blat A \code{data.frame}. Represents the added alignments of probes to the target sequences.
#' @param anno A \code{data.frame}. Represents the added annotation of the target sequences. If the parameter
#' 		  \code{anno} has not been defined then custom annotation will be populated with NA.
#' @param overwrite A \code{logical}. If FALSE the previous alignment table will be overwritten. By
#' 		  default it is set to FALSE.
#' @return Returns the \code{DEMIExperiment} object where the additional information has been added
#' 		   to the alignment table.
#' @details 
#' 
#' The user needs to make sure that the proper fields in the additional alignment information are not missing. To
#' see which fields are required use the function \code{colnames( getAlignment( x ) )} on the \code{DEMIExperiment}
#' object. The two fields that are always required are 'probeID' and 'targetID', the others are optional. All the
#' other fields will be generated automatically and set to NA. If the user knows the annotation information
#' for the targets in the alignment table then it is recommended to add that information to the annotation table.
#' To see what are the fields of annotation table use the function \code{colnames( getAnnotation( x ) )} on the
#' \code{DEMIExperiment} object. This information will then be seen later in the analysis results. When adding custom
#' annotations the only field that is required is the 'targetID', all other fields are optional however for later use
#' it is better to add some more information about the target like it's description.
#' 
#' @author Sten Ilmjarv
#' @examples 
#' \dontrun{
#' 
#' # To use the example we need to download a subset of CEL files from
#' # http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE9819 published
#' # by Pradervand et al. 2008.
#' 
#' # Set the destination folder where the downloaded files fill be located.
#' # It can be any folder of your choosing.
#' destfolder <- "demitest/testdata/"
#' 
#' # Download packed CEL files and change the names according to the feature
#' # they represent (for example to include UHR or BRAIN in them to denote the
#' # features).
#' # It is good practice to name the files according to their features which
#' # allows easier identification of the files later.
#' 
#' ftpaddress <- "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM247nnn"
#' download.file( paste( ftpaddress, "GSM247694/suppl/GSM247694.CEL.gz", sep = "/" ),
#' 		destfile = paste( destfolder, "UHR01_GSM247694.CEL.gz", sep = "" ) )
#' download.file( paste( ftpaddress, "GSM247695/suppl/GSM247695.CEL.gz", sep = "/" ),
#' 		destfile = paste( destfolder, "UHR02_GSM247695.CEL.gz", sep = "" ) )
#' download.file( paste( ftpaddress, "GSM247698/suppl/GSM247698.CEL.gz", sep = "/" ),
#' 		destfile = paste( destfolder, "UHR03_GSM247698.CEL.gz", sep = "" ) )
#' download.file( paste( ftpaddress, "GSM247699/suppl/GSM247699.CEL.gz", sep = "/" ),
#' 		destfile = paste( destfolder, "UHR04_GSM247699.CEL.gz", sep = "" ) )
#' download.file( paste( ftpaddress, "GSM247696/suppl/GSM247696.CEL.gz", sep = "/" ),
#' 		destfile = paste( destfolder, "BRAIN01_GSM247696.CEL.gz", sep = "" ) )
#' download.file( paste( ftpaddress, "GSM247697/suppl/GSM247697.CEL.gz", sep = "/" ),
#' 		destfile = paste( destfolder, "BRAIN02_GSM247697.CEL.gz", sep = "" ) )
#' download.file( paste( ftpaddress, "GSM247700/suppl/GSM247700.CEL.gz", sep = "/" ),
#' 		destfile = paste( destfolder, "BRAIN03_GSM247700.CEL.gz", sep = "" ) )
#' download.file( paste( ftpaddress, "GSM247701/suppl/GSM247701.CEL.gz", sep = "/" ),
#' 		destfile = paste( destfolder, "BRAIN04_GSM247701.CEL.gz", sep = "" ) )
#' 
#' # We need the gunzip function (located in the R.utils package) to unpack the gz files.
#' # Also we will remove the original unpacked files for we won't need them.
#' library( R.utils )
#' for( i in list.files( destfolder ) ) {
#' 	gunzip( paste( destfolder, i, sep = "" ), remove = TRUE )
#' }
#' 
#' # Now we can continue the example of the function demi
#' 

#' # Set up an experiment
#' demiexp <- DEMIExperiment(analysis = 'gene', celpath = destfolder,
#' 			experiment = 'myexperiment', organism = 'homo_sapiens')
#' 
#' # Create a custom annotation for one target ID
#' anno <- data.frame("RANDOM_ID", "Some kind of description")
#' colnames(anno) <- c("targetID", "description")
#' 
#' # Create a custom alignment for one target ID
#' blat <- data.frame(c(616302,1133177), c("RANDOM_ID","RANDOM_ID"))
#' colnames(blat) <- c("probeID", "targetID")
#' 
#' # Create a custom alignment for one target ID
#' customexp <- addCustomTargets(demiexp, blat, anno, FALSE)
#'  
#' }
#' 
#' 
#' @export 
#' @docType methods
#' @rdname addCustomTargets-methods
"addCustomTargets" <-
function( object = "DEMIExperiment",
		  blat = "data.frame",
		  anno = "data.frame",
		  overwrite = FALSE )
{
		
	analysis = getAnalysis( object );
	includestrand = TRUE;
	#	Check for correct analysis type
	if ( analysis == "gene" || analysis == "transcript" ) {
		
		#	Check if the blat table format is correct
		if ( length( grep( "probeID", colnames( blat ) ) ) == 0 ) {
			#stop( paste( "column", sQuote( "probeID" ), "is missing from the", sQuote( "blat" ), "object." ) );
			stop( DEMIMessages$customTargets$whatMissingFrom( "probeID", "blat" ) );
		}
		if ( length( grep( "targetID", colnames( blat ) ) ) == 0 ) {
			#stop( paste( "column", sQuote( "targetID" ), "is missing from the", sQuote( "blat" ), "object." ) );
			stop( DEMIMessages$customTargets$whatMissingFrom( "targetID", "blat" ) );
		}
		if ( length( grep( "pmsize", colnames( blat ) ) ) == 0 ) {
			#cat( paste( "column", sQuote( "pmsize" ), "is is unspecified in the", sQuote( "blat" ), "object. We will therefore assume that all probe matches to the target are perfect matches with 100% identity.\n" ) );
			cat( DEMIMessages$customTargets$pmsizeMissing );
		}
		if ( length( grep( "start", colnames( blat ) ) ) == 0 ) {
			#cat( paste( "column", sQuote( "start" ), "is is unspecified in the", sQuote( "blat" ), "object. Therefore you can't use these objects for imaging.\n" ) );
			cat( DEMIMessages$customTargets$startMissing );
		}
		if ( length( grep( "strand", colnames( blat ) ) ) == 0 ) {
			includestrand = FALSE;
			#cat( paste( "column", sQuote( "strand" ), "is is unspecified in the", sQuote( "blat" ), "object. We will therefore assume that all probe matches were on the main matching strand.\n" ) );
			cat( DEMIMessages$customTargets$strandMissing );
		}
		
		#	Check if all the probes in the added blat are present on the microarray, exit if not
		notpresent <- which( blat$probeID %in% rownames( getCelMatrix( object ) ) == FALSE );
		if ( length( notpresent ) > 0 ) {
			#stop( "The following probes were not present on the microarray. Remove them from your data and try the function again.\n", paste( notpresent, collapse = "," ) );
			stop( DEMIMessages$customTargets$probesNotPresent( notpresent ) );
		}
		
		#	Add missing columns to the blat so that the two blat tables could be bound
		for ( i in 1:length( colnames( getAlignment( object ) ) ) ) {
			coln <- colnames( getAlignment( object) )[i];
			if ( length( grep( coln, colnames( blat ) ) ) == 0 ) {
				blat <- data.frame( blat, NA );
				colnames( blat )[ length( colnames( blat ) ) ] = colnames( getAlignment( object) )[i]; 
			}
		}
		
		#	Check if annotation has also been added, else generate one
		tempanno <- NULL;
		if( missing( anno ) == FALSE ) {
			#	Check if the annotation format is correct
			if ( length( grep ("targetID", colnames( anno ) ) ) == 0 ) {
				#stop( paste( "Column", sQuote( "targetID" ), "is missing from the", sQuote( "anno" ), "object." ) );
				stop( DEMIMessages$customTargets$whatMissingFrom( "targetID", "anno" ) );
			}
			#	Check if gene analysis
			if ( getAnalysis( object ) == "gene" ) {
				if ( length( grep ("geneID", colnames( anno ) ) ) == 0 ) {
					#cat( paste( "column", sQuote( "geneID" ), "is missing from the", sQuote( "anno" ), "object. Will duplicate", sQuote( "targetID" ), "as", sQuote( "geneID" ), ".\n" ) );
					cat( DEMIMessages$customTargets$geneMissing );
					anno$geneID = anno$targetID;
				}
			}
			tempanno <- unique( anno );
		} else if ( missing( anno ) == TRUE ) {
			#cat( paste( sQuote("anno"), "has not been specified. Creating custom annotation\n" ) );
			cat( DEMIMessages$customTargets$annoMissing );
			tempanno$targetID <- unique( blat$targetID );
			if ( getAnalysis( object ) == "gene" ) {
				tempanno$geneID = tempanno$targetID;
			}
		}
		anno <- as.data.frame(tempanno);
		
		#	Add missing columns to the annotation
		for ( i in 1:length( colnames( getAnnotation( object ) ) ) ) {
			coln <- colnames( getAnnotation( object) )[i];
			if ( length( grep( coln, colnames( anno ) ) ) == 0 ) {
				anno <- data.frame( anno, NA );
				colnames( anno )[ length( colnames( anno ) ) ] = colnames( getAnnotation( object) )[i]; 
			}
		}
		anno <- anno[, colnames( getAnnotation( object ) ) ]; # select only columns in the annotation table
		
		#	Remove those probes that match to a strand with less matches except for 'genome' analysis
		#	since for genome it can match on both strands
		if ( analysis != "genome" ) {
			strand <- NULL;
			strand_table <- NULL;
			if ( overwrite == TRUE & includestrand == FALSE ) {
				strand_table <- NULL; # there is no strand
			} else if ( overwrite == FALSE & includestrand == FALSE ) {
				strand_table <- table( c( as.vector( getAlignment( object )$strand ) ) ); # determines the strand on original table
			} else if ( overwrite == TRUE & includestrand == TRUE ) {
				strand_table <- table( as.vector( blat$strand ) ); # determines the strand on the new table
			} else if ( overwrite == FALSE & includestrand == TRUE ) {
				strand_table <- table( c( as.vector( getAlignment( object )$strand ), as.vector( blat$strand ) ) ); # determines the strand on both tables
			}
			if ( !is.null( strand_table ) ) {
				if ( length( names( strand_table ) ) > 1 ) {
					if ( strand_table[ names( strand_table ) == "+" ] > strand_table[ names( strand_table ) == "-" ] ) {
						strand <- "+";
						#cat( "\tWill ignore all '-' and NA strand matches\n" );
						cat( DEMIMessages$customTargets$ignoreStrand( "-" ) );
					} else {
						strand <- "-";
						#cat( "\tWill ignore all '+' and NA strand matches\n" );
						cat( DEMIMessages$customTargets$ignoreStrand( "+" ) );
					}
				} else if ( length( names( strand_table ) ) == 1 ) {
					strand <- names( strand_table );
				}
			}
			#	Add main strand if strand was not included in the first place, otherwise the whole blat table will be left empty
			if ( !is.null( strand ) ) {
				blat$strand <- strand;
				blat <- blat[ blat$strand == strand, ];
			}
		}
		
		#	Remove probes whose perfect match is less then set by the user
		if ( !is.null( object@pmsize ) && length( grep( "pmsize", colnames( blat ) ) ) != 0 ) {
			blat <- blat[ blat$pmsize >= object@pmsize, ];
			#cat( "\tRemoving probes whose perfect match size is less then initially indicated when creating the DEMIExperiment object\n" );
			cat( DEMIMessages$customTargets$pmsizeSmallerRemove );
		}
		
		blat <- blat[, colnames( getAlignment( object ) ) ];
		
		#	Calculate number of maxtargets for the new probes added if maxtarget has been set
		outProbes <- NULL;
		if ( getMaxtargets( object ) > 0 ) {
			#cat( "\n\tDetermining the number of matches on targets for every probe" );
			cat( DEMIMessages$customTargets$determineMatches );
			
			orgblat <- getAlignment( object );
			tempblat <- rbind( orgblat[ orgblat$targetID %in% blat$targetID, ], blat );
			tempblat <- data.frame( tempblat[, c("probeID", "targetID")] );
			tempblat <- unique( tempblat );
			
			probeTargetCount <- data.frame( table( tempblat$probeID ) );
			outProbes <- as.vector( probeTargetCount[ probeTargetCount$Freq > object@maxtargets, c("Var1") ] );
			
		}

		#	Check if to overwrite or to add the data
		if ( overwrite == TRUE ) {
			object@blatTable <- blat;
			object@annoTable <- anno;
		} else if ( overwrite == FALSE ) {
			object@blatTable <- rbind( object@blatTable, blat );
			object@annoTable <- rbind( object@annoTable, anno );
		}
		
		#	Through out probes from the new blat table that have too many targets
		object@blatTable <- object@blatTable[ !( object@blatTable$probeID %in% outProbes ), ];
		
		#	Needs to be normalized again
		object = celMatrixNormalize( object, object@norm.method );
		
		#	Return the data
		return( object );
		
	} else {
		#stop( "You can't added custom region blat's and annotation to the exon analysis. You can only add them to 'transcript' or 'gene' analysis." );
		stop( DEMIMessages$customTargets$error );
	}
}

#' Functional annotation of DEMI results
#' 
#' The function \code{DEMIPathway} performs functional annotation analysis on DEMI differential expression results
#' stored in the \code{DEMIDiff} object. It takes into account the number of up- and down-regulated targets as well
#' as the total number of targets for each functional category to calculate the statistical significance of the
#' functional annotation. PS! This function can only be used if in the underlying \code{DEMIExperiment} object the \code{analysis}
#' paramater was set as 'gene' or 'transcript' for it will before functional annotation only on genes.
#' 
#' @param object A \code{DEMIDiff} object. The \code{DEMIDiff} object contains the results to differential
#' 		  expression analysis that will be used for functional annotation analysis.
#' @return Returns the results of the functional annotation analysis in a \code{data.frame}.
#' 
#' @author Sten Ilmjarv
#' @examples 
#' \dontrun{
#' 
#' # To use the example we need to download a subset of CEL files from
#' # http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE9819 published
#' # by Pradervand et al. 2008.
#' 
#' # Set the destination folder where the downloaded files fill be located.
#' # It can be any folder of your choosing.
#' destfolder <- "demitest/testdata/"
#' 
#' # Download packed CEL files and change the names according to the feature
#' # they represent (for example to include UHR or BRAIN in them to denote the
#' # features).
#' # It is good practice to name the files according to their features which
#' # allows easier identification of the files later.
#' 
#' ftpaddress <- "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM247nnn"
#' download.file( paste( ftpaddress, "GSM247694/suppl/GSM247694.CEL.gz", sep = "/" ),
#' 		destfile = paste( destfolder, "UHR01_GSM247694.CEL.gz", sep = "" ) )
#' download.file( paste( ftpaddress, "GSM247695/suppl/GSM247695.CEL.gz", sep = "/" ),
#' 		destfile = paste( destfolder, "UHR02_GSM247695.CEL.gz", sep = "" ) )
#' download.file( paste( ftpaddress, "GSM247698/suppl/GSM247698.CEL.gz", sep = "/" ),
#' 		destfile = paste( destfolder, "UHR03_GSM247698.CEL.gz", sep = "" ) )
#' download.file( paste( ftpaddress, "GSM247699/suppl/GSM247699.CEL.gz", sep = "/" ),
#' 		destfile = paste( destfolder, "UHR04_GSM247699.CEL.gz", sep = "" ) )
#' download.file( paste( ftpaddress, "GSM247696/suppl/GSM247696.CEL.gz", sep = "/" ),
#' 		destfile = paste( destfolder, "BRAIN01_GSM247696.CEL.gz", sep = "" ) )
#' download.file( paste( ftpaddress, "GSM247697/suppl/GSM247697.CEL.gz", sep = "/" ),
#' 		destfile = paste( destfolder, "BRAIN02_GSM247697.CEL.gz", sep = "" ) )
#' download.file( paste( ftpaddress, "GSM247700/suppl/GSM247700.CEL.gz", sep = "/" ),
#' 		destfile = paste( destfolder, "BRAIN03_GSM247700.CEL.gz", sep = "" ) )
#' download.file( paste( ftpaddress, "GSM247701/suppl/GSM247701.CEL.gz", sep = "/" ),
#' 		destfile = paste( destfolder, "BRAIN04_GSM247701.CEL.gz", sep = "" ) )
#' 
#' # We need the gunzip function (located in the R.utils package) to unpack the gz files.
#' # Also we will remove the original unpacked files for we won't need them.
#' library( R.utils )
#' for( i in list.files( destfolder ) ) {
#' 	gunzip( paste( destfolder, i, sep = "" ), remove = TRUE )
#' }
#' 
#' # Now we can continue the example of the function demi
#' 
#' # Set up an experiment
#' demiexp <- DEMIExperiment(analysis = 'gene', celpath = destfolder,
#' 			experiment = 'myexperiment', organism = 'homo_sapiens')
#' 
#' # Create clusters with an optimized wilcoxon's rank sum test incorporated within demi that
#' # precalculates the probabilities.
#' demiclust <- DEMIClust( demiexp, group = c( "BRAIN", "UHR" ), clust.method = demi.wilcox.test.fast )
#' 
#' # Calcuate differential expression
#' demidiff <- DEMIDiff( demiclust )
#' 
#' # Perform functiona annotation analysis on the DEMI analysis results
#' demipath <- DEMIPathway( demidiff )
#'  
#' }
#' 
#' @export 
#' @docType methods
#' @rdname DEMIPathway-methods
"DEMIPathway" <-
function( object = "DEMIDiff" )
{
	
	#cat( "*Performing the pathway analysis on the significant genes\n" );
	cat( DEMIMessages$DEMIPathway$main );
	
	if( getAnalysis( getExperiment( object ) ) == "gene" || getAnalysis( getExperiment( object ) ) == "transcript" ) {
		
		experiment <- getExperiment( object );
		#	Remove lines that have no GO:ID
		gos <- getPathway( experiment );
		#	Retrieve unique GO id's
		goIDs <- as.vector( unique( gos$accession ) );
		
		targetName <- NULL;
		if ( getAnalysis( experiment ) == "gene" ) {
			targetName <- "targetID";
		} else if ( getAnalysis ( experiment ) == "transcript" ) {
			targetName <- "geneID"
		}
		
		result <- getResultTable( object );
		clusters <- as.vector( unique( result$clusterID ) );
		output <- NULL;
		
		#	For every cluster calculate the number of significant targets
		for( i in 1:length( clusters ) ) {
			
			sgt <- result[ result$clusterID == clusters[i] & result$FDR <= 0.05, ]; # significantGenesTable
			genesSignificant <- length( unique( as.vector( sgt[ , grep( targetName, colnames( sgt ) ) ] ) ) );
			genesTotal <- length( unique( as.vector( result[ result$clusterID == clusters[i], targetName ] ) ) );
			
			genesInGO <- numeric( length( goIDs ) ); #	create empty vector
			genesSignificantInGO <- numeric( length( goIDs ) ); # create empty vector
			P.value <- numeric( length( goIDs ) ); # create empty vector
			
			go_id_index <- grep( "geneIDs", colnames( gos ) );
			
			for ( ii in 1:length( goIDs ) ) {
				#	get the genes corresponding to the GO id
				genes <- as.vector( strsplit( gos[ gos$accession == goIDs[ii], "geneIDs" ], "," )[[1]] );
				genesInGO[ii] <- length( genes );
				genesSignificantInGO[ii] <- length( which( unique( as.vector( sgt[,targetName] ) ) %in% genes ) );
				#P.value[ii] <- phyper( genesSignificantInGO[ii], genesSignificant, genesTotal - genesSignificant, genesInGO[ii], lower.tail = FALSE );
				P.value[ii] <- sum(dhyper( genesSignificantInGO[ii]:genesInGO[ii], genesSignificant, genesTotal - genesSignificant, genesInGO[ii] ))
			}
			FDR <- 1;
			output <- rbind( output, data.frame( clusters[i], goIDs, genesSignificantInGO, genesInGO, genesSignificant, genesTotal, P.value, FDR ) )
		}
		colnames( output )[1:2] = c( "clusterID", "GO:ID" );
		colnames( gos )[1] = "GO:ID";
#		output <- merge( output, unique( gos[, c(1,3,4,5) ] ), by = "GO:ID" )
		output <- plyr::join( output, unique( gos[, c(1,3,4,5) ] ), by = "GO:ID" )
		output$FDR <- p.adjust( output$P.value, "BY", length( output$P.value ) );
		output <- output[ order( output$FDR ), ];
		colnames( output )[ grep( "GO:ID", colnames( output ) ) ] <- "targetID";
		colnames( output )[ grep( "genesSignificantInGO", colnames( output ) ) ] <- "genesOnTargetInCluster";
		colnames( output )[ grep( "genesInGO", colnames( output ) ) ] <- "genesOnTarget";
		colnames( output )[ grep( "genesSignificant", colnames( output ) ) ] <- "genesInCluster";
		colnames( output )[ grep( "^name$", colnames( output ) ) ] <- "targetName";
		colnames( output )[ grep( "definition", colnames( output ) ) ] <- "description";
		output <- output[, c( "clusterID", "targetID", "genesOnTargetInCluster", "genesOnTarget", "genesInCluster", "genesTotal", "P.value", "FDR", "targetName", "description", "namespace" ) ]
		
		return( output );
	}
	
}

