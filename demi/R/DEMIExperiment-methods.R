#==============================================================================#
# DEMIExperiment-methods.R:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# initialize.DEMIExperiment
# validDemiExperiment
# getAnalysis
# getExperiment
# getCelpath
# getOrganism
# getArraytype
# getMaxtargets
# getMaxprobes
# getAnnotation
# getAlignment
# getCytoband
# getPathway
# getCelMatrix
# getNormMatrix
# getResult
# getResultTable
# getProbeLevel
# getTargetProbes
# getTargetAnnotation
# attachResult
# loadCel
# loadDEMILibrary
# loadAnnotation
# loadBlat
# loadCytoband
# check4probe
# check4target
# demisummary
# checkDEMIExperiment_analysis
# checkDEMIExperiment_experiment
# checkDEMIExperiment_celpath
# checkDEMIExperiment_pmsize
# checkDEMIExperiment_maxtargets
# checkDEMIExperiment_maxprobes
# checkDEMIExperiment_normalization
# DEMIExpression_hist
# probe.hist
# probe.plot
#==============================================================================#

#------------------------------------------------------------------------------#
# DEMIExperiment initialization:
#------------------------------------------------------------------------------#

#' Initializes the \code{DEMIExperiment} object
#' 
#' Initializes the \code{DEMIExperiment} object.
#' 
#' @param .Object A DEMIExperiment object.
#' @param ... Additional arguments that may never be used.
#' @return Returns a 'DEMIExperiment' object.
#' @author Sten Ilmjarv
#' @import methods
"initialize.DEMIExperiment" <-
function(.Object,...) 
{
	.Object <- callNextMethod(.Object,...);
	cat( DEMIMessages$demiStart() );
	.Object <- loadCel(.Object);
	.Object <- loadDEMILibrary(.Object);

	#	Normalize the data matrix
	#	- since CEL files are loaded, normalize them at once as well
	.Object <- celMatrixNormalize( .Object, .Object@norm.method );
	
	.Object;
	
}#initialize.DEMIExperiment

#' @import methods
setMethod( "initialize", "DEMIExperiment", initialize.DEMIExperiment );


#------------------------------------------------------------------------------#
# DEMIExperiment validation:
#------------------------------------------------------------------------------#

#' Validates the \code{DEMIExperiment} object
#' 
#' @param object A \code{DEMIExperiment} object.
#' @return Returns a validated \code{DEMIExperiment} object.
"validDEMIExperiment" <-
function( object )
{
	msg <- NULL;
	
	# check if 'analysis' exists
	analysis <- checkDEMIExperiment_analysis( object@analysis );
	if ( is.null( analysis ) == FALSE ) {
		msg <- paste( msg, "\n\tError:", analysis );
		return( msg );
	}
	# check for correct experiment definition
	experiment <- checkDEMIExperiment_experiment( object@experiment );
	if ( is.null( experiment ) == FALSE ) {
		msg <- paste( msg, "\n\tError:", experiment );
		return( msg );
	}
	# check for correct pmsize
	pmsize <- checkDEMIExperiment_pmsize( object@pmsize );
	if ( is.null( pmsize ) == FALSE ) {
		msg <- paste( msg, "\n\tError:", pmsize );
		return( msg );
	}
	# check for correct maxtargets
	maxtargets <- checkDEMIExperiment_maxtargets( object@maxtargets );
	if ( is.null( maxtargets ) == FALSE) {
		msg <- paste( msg, "\n\tError:", maxtargets );
		return( msg );
	}
	# check for correct maxprobes
	maxprobes <- checkDEMIExperiment_maxprobes( object@maxprobes );
	if ( is.null( maxprobes ) == FALSE) {
		msg <- paste( msg, "\n\tError:", maxprobes );
		return( msg );
	}
	# check for correct normalization method
	normalization <- checkDEMIExperiment_normalization( object@norm.method );
	if ( is.null( normalization ) == FALSE ) {
		msg <- paste( msg, "\n\tError:", normalization );
		return( msg );
	}
	# check for correct celpath
	celpath <- checkDEMIExperiment_celpath( object@celpath );
	if ( is.null( celpath ) == FALSE ) {
		msg <- paste( msg, "\n\tError:", celpath );
		return( celpath );
	}

	# do the verdict
	if ( is.null( msg ) ) TRUE else paste( "\n", msg )
	
}#validDEMIExperiment

#' @import methods
setValidity( "DEMIExperiment", validDEMIExperiment )#setValidity


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# DEMIExperiment get functions:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @rdname getAnalysis-methods
#' @aliases getAnalysis,DEMIExperiment-method
#' @import methods
setMethod( "getAnalysis", signature( object = "DEMIExperiment" ),
		function( object ) object@analysis
)#getAnalysis

#' @rdname getExperiment-methods
#' @aliases getExperiment,DEMIExperiment-method
#' @import methods
setMethod( "getExperiment", signature( object = "DEMIExperiment" ),
		function( object ) object@experiment
)#getExperiment

#' @rdname getCelpath-methods
#' @aliases getCelpath,DEMIExperiment-method
#' @import methods
setMethod( "getCelpath", signature( object = "DEMIExperiment" ),
		function( object ) object@celpath
)#getCelpath

#' @rdname getOrganism-methods
#' @aliases getOrganism,DEMIExperiment-method
#' @import methods
setMethod( "getOrganism", signature( object = "DEMIExperiment" ),
		function( object ) object@organism
)#getOrganism

#' @rdname getArraytype-methods
#' @aliases getArraytype,DEMIExperiment-method
#' @import methods
setMethod( "getArraytype", signature( object = "DEMIExperiment" ),
		function( object ) object@arraytype
)#getArraytype

#' @rdname getMaxtargets-methods
#' @aliases getMaxtargets,DEMIExperiment-method
#' @import methods
setMethod( "getMaxtargets", signature( object = "DEMIExperiment" ),
		function( object ) object@maxtargets
)#getMaxtargets

#' @rdname getMaxprobes-methods
#' @aliases getMaxprobes,DEMIExperiment-method
#' @import methods
setMethod( "getMaxprobes", signature( object = "DEMIExperiment" ),
		function( object ) object@maxprobes
)#getMaxprobes

#' @rdname getAnnotation-methods
#' @aliases getAnnotation,DEMIExperiment-method
#' @import methods
setMethod( "getAnnotation", signature( object = "DEMIExperiment" ),
		function( object ) object@annoTable
)#getAnnotation

#' @rdname getAlignment-methods
#' @aliases getAlignment,DEMIExperiment-method
#' @import methods
setMethod( "getAlignment", signature( object = "DEMIExperiment" ),
		function( object ) object@blatTable
)#getAlignment

#' @rdname getCytoband-methods
#' @aliases getCytoband,DEMIExperiment-method
#' @import methods
setMethod( "getCytoband", signature( object = "DEMIExperiment" ),
		function( object ) object@cytoband
)#getCytoband

#' @rdname getPathway-methods
#' @aliases getPathway,DEMIExperiment-method
#' @import methods
setMethod( "getPathway", signature( object = "DEMIExperiment" ),
		function( object ) object@pathway
)#getPathway

#' @rdname getCelMatrix-methods
#' @aliases getCelMatrix,DEMIExperiment-method
#' @import methods
setMethod( "getCelMatrix", signature( object = "DEMIExperiment" ),
		function( object ) object@exprsData@celMatrix
)#getCelMatrix

#' @rdname getNormMatrix-methods
#' @aliases getNormMatrix,DEMIExperiment-method
#' @import methods
setMethod( "getNormMatrix", signature( object = "DEMIExperiment" ),
		function( object ) object@exprsData@normMatrix
)#getNormMatrix

#' @rdname getResult-methods
#' @aliases getResult,DEMIExperiment-method
#' @import methods
setMethod( "getResult", signature( object = "DEMIExperiment" ),
		function( object ) object@results
)#getResults

#' @rdname getResultTable-methods
#' @aliases getResultTable,DEMIExperiment-method
#' @import methods
setMethod( "getResultTable", signature( object = "DEMIExperiment" ),
		function( object ) makeDEMIResultsTable( getResult( object ) )
)#getResultTable

#' @rdname getProbeLevel-methods
#' @aliases getProbeLevel,DEMIExperiment,vector,logical-method
#' @import methods
setMethod( "getProbeLevel", signature( object = "DEMIExperiment", probes = "vector", verbose = "logical" ),
		function( object, probes, verbose = TRUE ) {
			norm.matrix <- getNormMatrix( object );
			rows <- which( rownames( norm.matrix ) %in% probes );
			excludedProbes <- which( probes %in% rownames( norm.matrix ) == FALSE );
			excludedProbes <- probes[ excludedProbes ];
			if ( length( excludedProbes ) > 0 && verbose == TRUE ) {
				#cat( "Warning: these probes were excluded because they were not present in the analysis:\n", paste( excludedProbes, collapse = ", " ), "\n" )
				cat( DEMIMessages$DEMIExperiment$probesExcluded( excludedProbes ) );
			}
			return( norm.matrix[ rows, ] );
		}
)#getProbeLevel

#' @rdname getTargetProbes-methods
#' @aliases getTargetProbes,DEMIExperiment,vector-method
#' @import methods
setMethod( "getTargetProbes", signature( object = "DEMIExperiment", target = "vector" ),
		function( object, target ) {
			if ( length( target ) == 1 && nchar( target ) == 0 ) {
				#stop( paste( sQuote( "target" ), "can't be an empty string\n" ) );
				stop( DEMIMessages$notEmptyString( "target" ) );
			}
			annotable <- getAnnotation( object );
			blattable <- getAlignment( object );
			rows <- NULL;
			if ( getAnalysis( object ) == "gene" ) {
				rows <- which( tolower( annotable$geneID ) %in% tolower( target ) );
				rows <- c(rows, which( tolower( annotable$targetID ) %in% tolower( target ) ) );
				rows <- c(rows, which( tolower( annotable$geneSymbol ) %in% tolower( target ) ) );
				rows <- c(rows, which( tolower( annotable$peptideID ) %in% tolower( target ) ) );
			} else if ( getAnalysis( object ) == "exon" ) {
				rows <- which( tolower( annotable$geneID ) %in% tolower( target ) );
				rows <- c(rows, which( tolower( annotable$transcriptID ) %in% tolower( target ) ) );
				rows <- c(rows, which( tolower( annotable$targetID ) %in% tolower( target ) ) );
				rows <- c(rows, which( tolower( annotable$geneSymbol ) %in% tolower( target ) ) );
			} else if ( getAnalysis( object ) == "transcript" ) {
				rows <- which( tolower( annotable$geneID ) %in% tolower( target ) );
				rows <- c(rows, which( tolower( annotable$targetID ) %in% tolower( target ) ) );
				rows <- c(rows, which( tolower( annotable$geneSymbol ) %in% tolower( target ) ) );
			} else if ( getAnalysis( object ) == "genome" ) {
				rows <- c(rows, which( tolower( annotable$targetID ) %in% tolower( target ) ) );
			}
			rows <- unique( rows );
			targetIDs <- vector();
			if ( getAnalysis( object ) == "gene" ) {
				targetIDs <- as.vector( annotable[ rows, c( "geneID" )] );
			} else {
				targetIDs <- as.vector( annotable[ rows, c( "targetID" )] );
			}
			if ( length( targetIDs ) == 0 ) {
				#stop( "No targets were found\n" );
				stop( DEMIMessages$DEMIExperiment$noTargetsFound );
			} else {
				probes <- unique( blattable[ ( blattable$targetID %in% targetIDs ) == TRUE, c( "probeID", "targetID" ) ] );
				if ( getAnalysis( object ) == "gene" ) {
					addAnno <- unique( annotable[ which( annotable$geneID %in% probes$targetID ), c( "geneID", "geneSymbol" ) ] );
					colnames( addAnno )[ grep( "geneID", colnames( addAnno ) ) ] <- "targetID";
#					probes <- merge( probes, addAnno, by = "targetID" );
					probes <- plyr::join( probes, addAnno, by = "targetID" );
				} else if ( getAnalysis( object ) == "exon" ) {
					addAnno <- unique( annotable[ which( annotable$targetID %in% probes$targetID ), c( "targetID", "geneSymbol" ) ] );
#					probes <- merge( probes, addAnno, by = "targetID" );
					probes <- plyr::join( probes, addAnno, by = "targetID" );
				} else if ( getAnalysis( object ) == "transcript" ) {
					addAnno <- unique( annotable[ which( annotable$targetID %in% probes$targetID ), c( "targetID", "geneSymbol" ) ] );
#					probes <- merge( probes, addAnno, by = "targetID" );
					probes <- plyr::join( probes, addAnno, by = "targetID" );
				}
				probes <- unique( probes );
				if ( length( probes ) == 0 ) {
					#stop( "No probes measuring the specified targets matched the experiment criteria\n" );
					stop( DEMIMessages$DEMIExperiment$noProbesMeasuring );
				} else {
					return( probes );
				}
			}
		}
)#getTargetProbes

#' @rdname getTargetAnnotation-methods
#' @aliases getTargetAnnotation,DEMIExperiment,vector-method
#' @import methods
setMethod( "getTargetAnnotation", signature( object = "DEMIExperiment", target = "vector" ),
		function( object, target ) {
			if ( length( target ) > 0 ) {
				annotable <- getAnnotation( object );
				rows <- NULL;
				if ( getAnalysis( object ) == "gene" ) {
					rows <- which( tolower( annotable$geneID ) %in% tolower( target ) );
					rows <- c(rows, which( tolower( annotable$targetID ) %in% tolower( target ) ) );
					rows <- c(rows, which( tolower( annotable$geneSymbol ) %in% tolower( target ) ) );
					rows <- c(rows, which( tolower( annotable$peptideID ) %in% tolower( target ) ) );
				} else if ( getAnalysis( object ) == "exon" ) {
					rows <- which( tolower( annotable$geneID ) %in% tolower( target ) );
					rows <- c(rows, which( tolower( annotable$transcriptID ) %in% tolower( target ) ) );
					rows <- c(rows, which( tolower( annotable$targetID ) %in% tolower( target ) ) );
					rows <- c(rows, which( tolower( annotable$geneSymbol ) %in% tolower( target ) ) );
				} else if ( getAnalysis( object ) == "transcript" ) {
					rows <- which( tolower( annotable$geneID ) %in% tolower( target ) );
					rows <- c(rows, which( tolower( annotable$targetID ) %in% tolower( target ) ) );
					rows <- c(rows, which( tolower( annotable$geneSymbol ) %in% tolower( target ) ) );
				} else if ( getAnalysis( object ) == "genome" ) {
					rows <- c(rows, which( tolower( annotable$targetID ) %in% tolower( target ) ) );
				}
				rows <- unique( rows );
				return( annotable[rows, ] );
			}
		}
)#getTargetAnnotation


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# DEMIExperiment set functions:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @rdname attachResult-methods
#' @aliases attachResult,DEMIExperiment,DEMIDiff-method
#' @import methods
setMethod( "attachResult", signature( object = "DEMIExperiment", diffObject = "DEMIDiff" ),
		function( object, diffObject ) {
			result.experiment <- getExperiment( diffObject );
			if ( identical( getNormMatrix( object ), getNormMatrix( result.experiment ) ) == TRUE &&
					identical( getExperiment( object ), getExperiment( result.experiment ) ) == TRUE &&
					identical( getCelMatrix( object ), getCelMatrix( result.experiment ) ) == TRUE ) {
				canAdd <- TRUE;
				for ( i in 1:length( getResult( object ) ) ) {
					if ( identical( getResult( diffObject ), getResult( object )[[ getName( diffObject ) ]] ) == TRUE ) {
						#cat( paste( "Warning: Can't add the variable", sQuote( deparse( substitute( diffObject ) ) ), "of class 'DEMIDiff' to the object of class 'DEMIExperiment' because the results of that object already exists in the 'DEMIExperiment' object\n" ) );
						canAdd <- FALSE;
					}
				}
				if ( canAdd == FALSE ) {
					cat( DEMIMessages$DEMIExperiment$attachErrorResultsExists( diffObject ) );
				}
				else if ( canAdd == TRUE ) {
					object@results[ getName( diffObject ) ] <- getResult( diffObject );
				}
				return( object );
			} else {
				#stop( "Can't add results of class 'DEMIDiff' to this experiment of class 'DEMIExperiment' because they implement different 'DEMIExperiment' objects" );
				stop( DEMIMessages$DEMIExperiment$attachErrorDiffDEMIExp )
			}
		}
)#attachResult


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# DEMIExperiment load functions:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @rdname loadCel-methods
#' @aliases loadCel,DEMIExperiment-method
#' @import methods
setMethod( "loadCel", signature( object = "DEMIExperiment" ),
		function( object ) {
			#	If 'celpath' has been defined then load the cel files
			if ( is.null( checkDEMIExperiment_celpath( object@celpath ) ) == TRUE ) {
				#cat( "*Loading CEL files and creating expression matrix\n" );
				cat( DEMIMessages$DEMIExperiment$loadingCEL );
				#	Read in original data matrix
				orgdata <- NULL;
				options( warn = -1 ) # Turn warnings off
				if ( length( object@celpath ) == 1 ) {
					if ( R.utils::isDirectory( object@celpath ) == TRUE ) {	#	if celpath is a directory
						#cat( paste( "\tThe CEL files included are: ", paste( affy::list.celfiles( object@celpath ), collapse = ", ", sep = "" ), "\n", sep = "" ) );
						cat( DEMIMessages$DEMIExperiment$includedCELFiles( paste( affy::list.celfiles( object@celpath ), collapse = ", ", sep = "" ) ) );
#						orgdata <- ReadAffy( celfile.path = object@celpath );
						orgdata <- oligo::read.celfiles( affy::list.celfiles( object@celpath, full.names = TRUE ) )
                        object@arraytype <- affy::whatcdf( paste( object@celpath, "/", dir( object@celpath )[1], sep = "" ) );
						
					} else {	#	if celpath is a list of files
						#cat( paste( "\tThe CEL files included are: ", paste( object@celpath, collapse = ", ", sep = "" ), "\n", sep = "" ) );
						cat( DEMIMessages$DEMIExperiment$includedCELFiles( paste( object@celpath, collapse = ", ", sep = "" ) ) );
#						orgdata <- ReadAffy( filenames = as.vector( object@celpath ) );
						orgdata <- oligo::read.celfiles( filenames = object@celpath )
						object@arraytype <- affy::whatcdf( ( object@celpath )[1] );
					}
				} else {
					#cat( paste( "\tThe CEL files included are: ", paste( object@celpath, collapse = ", ", sep = "" ), "\n", sep = "" ) );
					cat( DEMIMessages$DEMIExperiment$includedCELFiles( paste( object@celpath, collapse = ", ", sep = "" ) ) );
#					orgdata <- ReadAffy( filenames = as.vector( object@celpath ) );
					orgdata <- oligo::read.celfiles( filenames = object@celpath )
					object@arraytype <- affy::whatcdf( ( object@celpath )[1] );
				}
				options( warn = 1 ) # Turn warnings back on
                orgdata <- affy::exprs( orgdata );
				exprsData <- DEMICel( celMatrix = orgdata );
				object@exprsData <- exprsData;
				#cat( paste( "\tUsing microarray platform ", object@arraytype, "\n", sep = "" ) );
				cat( DEMIMessages$DEMIExperiment$usingMicroarrayPlatform( object@arraytype ) );
				return( object );
			} else {
				#stop( paste( "Error: ", sQuote( "celpath" ), " has not been defined\n", sep = "" ) );
				stop( DEMIMessages$parameterMissing( "celpath" ) );
			}
		}
)#loadCel

#' @rdname loadDEMILibrary-methods
#' @aliases loadDEMILibrary,DEMIExperiment-method
#' @import methods
setMethod( "loadDEMILibrary", signature( object = "DEMIExperiment" ),
		function( object ) {
			pkg <- paste( "demi", cleanorganismname( getOrganism( object ) ), sep = "" );
			#pkg <- paste( "demi", cleanorganismname( "homo_sapiens" ), sep = "" );
			#	check if the library has been loaded - if not load the library
			if ( !( pkg %in% loadedNamespaces() ) ) {
				#	check if library is installed
				if ( length( find.package( pkg, quiet = TRUE ) ) == 0 ) {	# install the library from the repository
					#cat( paste( "Library - package", pkg, "not installed" ) );
					cat( DEMIMessages$DEMIExperiment$libraryNotInstalled( pkg ) );
					#	check if possible to install the package
					#cat( "Now searching the internet repository" );
					cat( DEMIMessages$DEMIExperiment$searchingInternetRepo( pkg ) );
					if ( capabilities()["http/ftp"] ) {
						# check version
						demipkgs <- NULL;
						ver <- paste( version$major, version$minor, sep = ".")
						ver <- sub( "\\.\\d+$", "", rev(ver) )
                        ########################
                        # Now we will do an independent repository, not related to R versions using devtools with specified url
                        # --- Below is old code from version 1.1.1 ---
                        #demipkgs <- available.packages( contriburl = "http://biit.cs.ut.ee/demi/R/src/contrib/" );
                        #if ( !( pkg %in% demipkgs[, "Package"] ) ) {
                        #	#	the package was not found
                        #	#stop( paste( "Package", pkg, "was not found in the http://biit.cs.ut.ee/demi/R/src/contrib/." ) )
                        #	stop( DEMIMessages$DEMIExperiment$packageNotFoundFromRepo( pkg ) );
                        #} else {
                        #	#	install the package
                        #	cat( DEMIMessages$DEMIExperiment$installingPackage( pkg, "http://biit.cs.ut.ee/demi/R" ) );
                        #	install.packages( pkg, repos="http://biit.cs.ut.ee/demi/R", dependencies = FALSE )
                        #}
                        # --- Now will  come new code for version >= 1.1.2 ---
                        ########################
                        # first check if the specified url exist
                        try_install <- try(devtools::install_url(paste("http://biit.cs.ut.ee/demi/R/datapackages/", pkg, "_0.1.tar.gz", sep="")), silent = T)
                        if (class(try_install) == "try-error") {
                            stop( DEMIMessages$DEMIExperiment$packageNotFoundFromRepo( pkg ) );
                        }
					} else {
						#stop( paste( "The current operation could not access the internet. Please",
						#			 "check your internet connection." ) )
						stop( DEMIMessages$DEMIExperiment$checkYourInternet );
					}
				}
				#	load the installed library
				do.call("library", list( pkg, character.only=TRUE ) )
			} # package is already loaded
			
			#	Retrieve the data environment [ might be better solution, might be not ]
			data.env <- pos.to.env( match( paste( "package:", pkg, sep = "" ), search() ) );
			
			#	Load the annotation data from the package
			object@annoTable <- loadAnnotation( object, data.env );
			
			#	Load blat table
			object@blatTable <- loadBlat( object, data.env );
			
			#	Load the cytoband data
			if ( getAnalysis( object ) == "genome" ) {
				object@cytoband <- loadCytoband( object, data.env );
			}
			
			#	Load the pathway data
			if ( getAnalysis( object ) == "gene" || getAnalysis( object ) == "transcript" ) {
				object@pathway <- loadPathway( object, data.env );
			}

			unloadNamespace( pkg );
			gc( verbose=FALSE );

			return( object );
			
		}
)#loadDEMILibrary

#' @rdname loadAnnotation-methods
#' @aliases loadAnnotation,DEMIExperiment,environment-method
#' @import methods
setMethod( "loadAnnotation", signature( object = "DEMIExperiment", pkg = "environment" ),
		function( object, pkg ) {
			annoTable <- NULL;
			anno <- paste( getPackageName( pkg ), "anno", sep = "" );
			if ( anno %in% do.call( "data", list( package = getPackageName( pkg ), envir = pkg ) )$results[,3] ) {
				#	check if the data is already loaded. It can't be however because it is loaded in a temporary environment
				
				#if( !exists( anno, inherits = FALSE ) ) {
				
				#	load the package
				#cat( paste( "*Loading required annotation data from", getPackageName( pkg ), "\n", sep = "" ) );
				cat( DEMIMessages$DEMIExperiment$loadingAnnoSuccess( getPackageName( pkg ) ) );
				# create a new environment for loading, since loading to .GlobalEnv is not allowed
				temp.env <- new.env();
				data( list = anno, envir = temp.env );	# can't load into package environment because for some reason the environment is locked. ask somebody
				annoTable <- get( anno, envir = temp.env );
				#}
				#else {
				#	annoTable <- get( anno );
				#}
				
			} else {
				#stop( paste( "Can't load required annotation table", anno, "from the package", getPackageName( pkg ) ) );
				stop( DEMIMessages$DEMIExperiment$loadingAnnoFail( anno, getPackageName( pkg ) ) );
			}
			#	if it is able to load the annoTable
			analysis <- getAnalysis( object );
			if ( analysis == "exon" && !is.null( exome( annoTable ) ) ) {
				return( exome( annoTable ) );
			} else if ( ( analysis == "transcript" || analysis == "gene" ) && !is.null( transcriptome( annoTable ) ) ) {
				return( transcriptome( annoTable ) );
			} else if ( analysis == "genome" ) {
				if ( !is.null( genome( annoTable, object@sectionsize ) ) ) {
					return( genome( annoTable, object@sectionsize ) );
				} else {
					#stop( paste( "The specified", sQuote( "sectionsize" ), "is not available. Try", paste( sectionsizes( annoTable ), collapse = ", " ) ) );
					stop( DEMIMessages$DEMIExperiment$sectionsizeNotAvail( sectionsizes( annoTable ) ) );
				}
			} else {
				#stop( "Can't load annotation table with the specified parameters" );
				stop( DEMIMessages$DEMIExperiment$cantLoadWhatTable( "anno" ) );
			}
			
		}
)#loadAnnotation

#' @rdname loadBlat-methods
#' @aliases loadBlat,DEMIExperiment,environment-method
#' @import methods
setMethod( "loadBlat", signature( object = "DEMIExperiment", pkg = "environment" ),
		function( object, pkg ) {
			
			analysis <- getAnalysis( object );
			blatdata <- NULL;
            blat <- affy::cleancdfname( getArraytype( object ), addcdf = FALSE );
			if ( blat %in% do.call( "data", list( package = getPackageName( pkg ), envir = pkg ) )$results[,3] ) {
				#	check if the data is already loaded
				
#				if( !exists( blat, inherits = FALSE ) ) {
				#	load the package
				#cat( paste( "*Loading required blat data from ", getPackageName( pkg ), "\n", sep = "" ) );
				cat( DEMIMessages$DEMIExperiment$loadingBlatSuccess( getPackageName( pkg ) ) );
				temp.env <- new.env();
				data( list = blat, envir = temp.env );	# can't load into package environment because for some reason the environment is locked. ask somebody
				blatdata <- get( blat, envir = temp.env );
				
#				}
#				else {
#					blatdata <- get( blat );
#				}
				
			} else {
				#stop( paste( "The required blat table", blat, "for the array", getArraytype( object ),"does not exists in the package", getPackageName( pkg ) ) );
				stop( DEMIMessages$DEMIExperiment$blatDoesNotExists( blat, getArraytype( object ), getPackageName( pkg ) ) );
			}
			
			#	check for correct perfect match sizes
			if ( ! object@pmsize %in% pmsizes( blatdata ) ) {
				#stop( paste( "The specified ", sQuote( "pmsize" ), " is not available. Available pmsize's are ", paste( pmsizes( blatdata ), collapse = ", " ) , sep = "" ) );
				stop( DEMIMessages$DEMIExperiment$pmsizeNotAvail( pmsizes( blatdata ) ) );
			}
			
			#	load the blat data
			blatTable <- NULL;
			maxtargetdata <- NULL;	# PS! maxtargets are only calculated on a specific strand of probes if not 'genome' analysis
			if ( analysis == "exon" ) {
				blatTable <- exome( blatdata )
				maxtargetdata <- maxtarget( blatdata, object@pmsize, "exome" );
			} else if ( analysis == "transcript" || analysis == "gene" ) {
				blatTable <- transcriptome( blatdata );
				maxtargetdata <- maxtarget( blatdata, object@pmsize, "transcriptome" );
			} else if ( analysis == "genome" ) {
				#	check if sectionsize is available
				if ( !is.null( genome( blatdata, object@sectionsize ) ) ) {
					blatTable <- genome( blatdata, object@sectionsize );
					maxtargetdata <- maxtarget( blatdata, object@pmsize, "genome", object@sectionsize );
				}
			} else {
				#stop( "Can't load blat table with the specified parameters" );
				stop( DEMIMessages$DEMIExperiment$cantLoadWhatTable( "blat" ) );
			}
			
			#	Remove those probes that match to a strand with less matches except for 'genome' analysis
			#	since for genome it can match on both strands
			if ( analysis != "genome" ) {
				strand <- NULL;
				strand_table = table( blatTable$strand );
				if ( strand_table[ names( strand_table ) == "+" ] > strand_table[ names( strand_table ) == "-" ] ) {
					strand <- "+";
					#cat( "\tWill ignore the '-' strand matches\n" );
					cat( DEMIMessages$DEMIExperiment$ignoreStrandMatches( "-" ) );
				} else {
					strand <- "-";
					#cat( "\tWill ignore the '+' strand matches\n" );
					cat( DEMIMessages$DEMIExperiment$ignoreStrandMatches( "+" ) );
				}
				blatTable <- blatTable[ blatTable$strand == strand, ];
			}
			
			#	Remove probes whose perfect match is less then set by the user
			if ( !is.null( object@pmsize ) ) {
				blatTable <- blatTable[ blatTable$pmsize >= object@pmsize, ];
			}
			
			#	Remove probes who match to more targets then set by the user
			if ( getMaxtargets( object ) > 0 && !is.null( object@annoTable ) ) {
				#	Read in the maxtarget file
				outProbes <- as.vector( maxtargetdata[ maxtargetdata$targets > getMaxtargets( object ), c( "probeID" ) ] );
				blatTable <- blatTable[ !( blatTable$probeID %in% outProbes ), ];
			}
			
			#	Make a proper blat table for gene analysis
			if ( analysis == "gene" ) {
#				tempBlat <- merge( blatTable[, c( "probeID", "targetID", "start", "strand" )], getAnnotation( object )[, c( "geneID", "targetID" )], by = "targetID" );
				tempBlat <- plyr::join( blatTable[, c( "probeID", "targetID", "start", "strand" )], getAnnotation( object )[, c( "geneID", "targetID" )], by = "targetID" );
				tempBlat <- unique( tempBlat[, c( "probeID", "geneID", "start", "strand" )] );
				colnames( tempBlat )[ grep( "geneID", colnames( tempBlat ) ) ] = "targetID";
				blatTable <- tempBlat;
			}
			
			blatTable <- unique( blatTable ); # recently added - just in case there are duplicacies
			
			return( blatTable );
			
		}
)#loadBlat

#' @rdname loadCytoband-methods
#' @aliases loadCytoband,DEMIExperiment,environment-method
#' @import methods
setMethod( "loadCytoband", signature( object = "DEMIExperiment", pkg = "environment" ),
		function( object, pkg ) {
			
			#cat ( "*Loading cytoband information\n" );
			cat( DEMIMessages$DEMIExperiment$loadingCytoband );
			
			annoTable <- NULL;
			anno <- paste( getPackageName( pkg ), "anno", sep = "" );
			if ( anno %in% do.call( "data", list( package = getPackageName( pkg ), envir = pkg ) )$results[,3] ) {
				#	check if the data is already loaded
#				if( !exists( anno, inherits = FALSE ) ) {

				#	load the package
				#cat( paste( "*Loading required data from ", getPackageName( pkg ), "\n", sep = "" ) );
				temp.env <- new.env()
				data( list = anno, envir = temp.env );	# can't load into package environment because for some reason the environment is locked. ask somebody
				annoTable <- get( anno, envir = temp.env );
				
				return( cytoband( annoTable ) );
				
#				}
#				else {
#					annoTable <- get( anno );
#				}
				
			} else {
				#stop( paste( "Can't load required annotation table", anno, "from the package", getPackageName( pkg ) ) );
				stop( DEMIMessages$DEMIExperiment$loadingAnnoFail( anno, getPackageName( pkg ) ) );
			}
			
		}
)#loadCytoband

#' @rdname loadPathway-methods
#' @aliases loadPathway,DEMIExperiment,environment-method
#' @import methods
setMethod( "loadPathway", signature( object = "DEMIExperiment", pkg = "environment" ),
		function( object, pkg ) {
			
			#cat ( "*Loading pathway information\n" );
			cat( DEMIMessages$DEMIExperiment$loadingPathway );
			
			annoTable <- NULL;
			anno <- paste( getPackageName( pkg ), "anno", sep = "" );
			if ( anno %in% do.call( "data", list( package = getPackageName( pkg ), envir = pkg ) )$results[,3] ) {
				#	check if the data is already loaded
#				if( !exists( anno, inherits = FALSE ) ) {
				#	load the package
				#cat( paste( "*Loading required data from ", getPackageName( pkg ), "\n", sep = "" ) );
				temp.env <- new.env();
				data( list = anno, envir = temp.env );	# can't load into package environment because for some reason the environment is locked. ask somebody
				annoTable <- get( anno, envir = temp.env );
					
#				}
#				else {
#					annoTable <- get( anno );
#				}
				
				return( pathway( annoTable ) );
				
			} else {
				#stop( paste( "Can't load required annotation table", anno, "from the package", getPackageName( pkg ) ) );
				stop( DEMIMessages$DEMIExperiment$loadingAnnoFail( anno, getPackageName( pkg ) ) );
			}
			
		}
)#loadPathway

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#	DEMIExperiment other functions:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @rdname check4probe-methods
#' @aliases check4probe,DEMIExperiment,vector-method
#' @import methods
setMethod( "check4probe", signature( object = "DEMIExperiment", probes = "vector" ),
		function( object, probes ) {
			msg <- NULL;
			blatTable <- getAlignment( object );
			difference <- length( setdiff( probes, blatTable$probeID ) );
			if( difference > 0 ) {
				#msg <- paste( "\tError:\tThere were", difference, "probes that were not located in the blat table.\n" );
				#msg <- paste( msg, "\t\tOnly probes in the blat table can be used to populate custom probe vector cluster\n" );
				#msg <- paste( msg, "\t\tUse the function 'getAlignment' on your 'DEMIExperiment' object to view probes present in blat table\n" );
				stop( DEMIMessages$DEMIExperiment$check4probeError( difference ) );
			}
			return( msg );
		}
)#check4probe

#' @rdname check4target-methods
#' @aliases check4target,DEMIExperiment,vector-method
#' @import methods
setMethod( "check4target", signature( object = "DEMIExperiment", target = "vector" ),
		function( object, target ) {
			targetProbes <- getTargetProbes( object, target );
			if ( nrow( targetProbes ) > 0 ) {
				# check if some target ID's do not exists in the annotation table
				annotable <- getAnnotation( object )
				rows <- NULL
				if ( getAnalysis( object ) == "gene" ) {
					rows <- which( annotable$geneID %in% targetProbes$targetID );
					rows <- c(rows, which( annotable$targetID %in% targetProbes$targetID ) );
					rows <- c(rows, which( annotable$geneSymbol %in% targetProbes$targetID ) );
					rows <- c(rows, which( annotable$peptideID %in% targetProbes$targetID ) );
				} else if ( getAnalysis( object ) == "exon" ) {
					rows <- which( annotable$geneID %in% targetProbes$targetID );
					rows <- c(rows, which( annotable$transcriptID %in% targetProbes$targetID ) );
					rows <- c(rows, which( annotable$targetID %in% targetProbes$targetID ) );
					rows <- c(rows, which( annotable$geneSymbol %in% targetProbes$targetID ) );
				} else if ( getAnalysis( object ) == "transcript" ) {
					rows <- which( annotable$geneID %in% targetProbes$targetID );
					rows <- c(rows, which( annotable$targetID %in% targetProbes$targetID ) );
					rows <- c(rows, which( annotable$geneSymbol %in% targetProbes$targetID ) );
				} else if ( getAnalysis( object ) == "genome" ) {
					rows <- c(rows, which( annotable$targetID %in% targetProbes$targetID ) );
				}
				annotable <- annotable[ unique( rows ), ]
				# check if all the targets exists in the annotable
				if ( getAnalysis( object ) == "gene" ) {
					diff <- intersect( target, annotable$geneID );
					diff <- c( diff, intersect( target, annotable$targetID ) );
					diff <- c( diff, intersect( target, annotable$geneSymbol ) );
					diff <- c( diff, intersect( target, annotable$peptideID ) );
				} else if ( getAnalysis( object ) == "exon" ) {
					diff <- intersect( target, annotable$geneID );
					diff <- c( diff, intersect( target, annotable$transcriptID ) );
					diff <- c( diff, intersect( target, annotable$targetID ) );
					diff <- c( diff, intersect( target, annotable$geneSymbol ) );
				} else if ( getAnalysis( object ) == "transcript" ) {
					diff <- intersect( target, annotable$geneID );
					diff <- c( diff, intersect( target, annotable$targetID ) );
					diff <- c( diff, intersect( target, annotable$geneSymbol ) );
				} else if ( getAnalysis( object ) == "genome" ) {
					diff <- intersect( target, annotable$targetID );
				}
				if ( length( setdiff( target, diff ) ) == 0 ) {
					return( TRUE );
				} else {
					stop( DEMIMessages$DEMIExperiment$check4targetError( setdiff( target, diff ) ) );
				}
			} else {
				#stop( "0 probes were found for the specified targets\n" );
				stop( DEMIMessages$DEMIExperiment$check4targetError( target ) );
			}
		}
)#check4target

#' @rdname demisummary-methods
#' @aliases demisummary,DEMIExperiment-method
#' @import methods
setMethod( "demisummary", signature( object = "DEMIExperiment" ), 
		function( object, target )
		{
			#	check that the function is correcty run
			if ( missing( target ) == TRUE ) {
				#stop( paste( "parameter ", sQuote( "target" ), " is missing", sep = "" ) );
				stop( DEMIMessages$parameterMissing( "target ") );
			} else {
				if ( is.vector( target ) == FALSE ) {
					#stop( paste( "parameter ", sQuote( "target" ), " needs to be a vector", sep = "" ) );
					stop( DEMIMessages$isNotError( "target", "vector" ) );
				}
			}
			# 	if any targets match the criteria
			if ( isTRUE( check4target( object, target ) ) ) {
				#	find the targets
				targetProbes <- droplevels( getTargetProbes( object, target ) );
				# create a data.frame from target names
				targetID <- as.vector( unique( droplevels( targetProbes )$targetID ) );
				output <- data.frame( targetID );
#				targetProbes <- as.matrix( targetProbes );
				result <- getResult( object );
				targetList <- tapply( targetProbes[, grep( "probeID", colnames( targetProbes ) )], targetProbes[, grep( "targetID", colnames( targetProbes ) )], function(x) { return( x ) } );
				#	find if groups are specified
				if ( length( result ) > 0 ) {
					# for each group check if index's exists, otherwise can't produce demisummary over the CEL files
					for ( i in 1:length( result ) ) {
						group <- getGroup( result[[i]] );
						if ( length( getIndexA( group ) ) > 0 || length( getIndexB( group ) ) > 0 ) {
							groupA <- getGroupA( group );
							groupB <- getGroupB( group );
							output <- data.frame( output, numeric( nrow( output) ), numeric( nrow( output ) ) );
							colnames( output )[ c( ncol( output ) - 1, ncol( output ) ) ] <- c( groupA, groupB ); 
							for( ii in 1:length( targetList ) ) {
								targetName <- names( targetList[ii]);
								probeLevel <- getProbeLevel( object, targetList[[ii]], verbose = FALSE );
								#output[ output$targetID == targetName, grep( groupA, colnames( output ) ) ] <- mean( probeLevel[ , getIndexA( group ) ] );
								#output[ output$targetID == targetName, grep( groupB, colnames( output ) ) ] <- mean( probeLevel[ , getIndexB( group ) ] );
								output[ output$targetID == targetName, grep( groupA, colnames( output ) ) ] <- mean( probeLevel[ getIndexA( group ) ] );
								output[ output$targetID == targetName, grep( groupB, colnames( output ) ) ] <- mean( probeLevel[ getIndexB( group ) ] );
							}
						}
					}
				}
				# output the whole mean of normalized matrix
				ALL <- do.call( "rbind", lapply( targetList, function( x ) {
									probeLevel <- getProbeLevel( object, x, verbose = FALSE );
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
				stop( DEMIMessages$DEMIExperiment$zeroTargetsFound );
			}
		}
);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#	DEMIExperiment validation functions:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' Checks if the \code{analysis} is correct
#' 
#' Checks if the \code{analysis} is correct for DEMI analysis. It can be either 'gene', 'transcript',
#' 'exon' or 'genome'. It is used internally in DEMI analysis.
#' 
#' @param analysis A \code{character}.
#' @return Returns NULL if \code{analysis} is ok, else returns an error message.
#' @author Sten Ilmjarv
checkDEMIExperiment_analysis = function( analysis )
{
	msg <- NULL;
	analysisTypes <- c( "transcript", "gene", "exon", "genome" );
	if( is.na( match( analysis, analysisTypes ) ) == TRUE ) {
		#msg <- paste( "invalid", sQuote( "analysis" ), "type. You need to choose", sQuote( "analysis" ), "type from (", paste( analysisTypes, collapse = ", " ),")\n" )
		msg <- DEMIMessages$DEMIExperiment$invalidWhatChooseFrom( "analysis", analysisTypes );
	}
	return( msg );
}#checkDEMIExperiment_analysis

#' Checks if the \code{experiment} is correct
#' 
#' Checks if the \code{experiment} is correct for the DEMI analysis. It is used internally in DEMI
#' analysis.
#' 
#' @param experiment A \code{character}.
#' @return Returns NULL if \code{experiment} is ok, else returns an error message.
#' @author Sten Ilmjarv
checkDEMIExperiment_experiment = function( experiment )
{
	msg <- NULL;
	
	if ( length( experiment ) == 0 ) {
		#msg <- paste( sQuote( "experiment" ), "is unspecified\n" );
		msg <- DEMIMessages$parameterMissing( "experiment" );
	}
	else if ( length( experiment ) > 1 ) {
		#msg <- paste( sQuote( "experiment" ), "can only have name\n" );
		msg <- DEMIMessages$tooManyParameters( "experiment", 1 );
	}
	
	return( msg );
}#checkDEMIExperiment_experiment

#' Checks if \code{celpath} is correct
#' 
#' Checks if \code{celpath} is correct for DEMI analysis. The \code{celpath} stores
#' CEL directory or CEL files. It is used internally in DEMI analysis.
#' 
#' @param celpath A \code{character}.
#' @return Returns NULL if \code{celpath} is ok, else returns an error message.
#' @author Sten Ilmjarv
checkDEMIExperiment_celpath <- function( celpath )
{
	msg <- NULL;
	if ( length( celpath ) > 0 ) {
		if ( length( celpath ) == 1 ) {
            if ( R.utils::isDirectory( celpath ) == TRUE ) {
				if ( file.exists( celpath ) == FALSE ) {
					#msg <- paste( sQuote( "celpath" ), "does not exists\n" );
					msg <- DEMIMessages$DEMIExperiment$notFolder( "celpath" );
				} else if ( file.exists( celpath ) == TRUE ) {
					files <- grep( ".CEL$", dir( celpath ) );
					if ( length( files ) == 0 ) {
						#msg <- paste( sQuote( "celpath" ), "does not contain any CEL files - no CEL suffix found\n" );
						msg <- DEMIMessages$DEMIExperiment$folderEmpty( "celpath" );
					}
				}
			} else {
                if ( affxparser::isCelFile( celpath ) == FALSE ) {
					#msg <- paste( "The ", celpath, " in ", sQuote( "celpath" ), " is not a CEL file!\n", sep = "" );
					msg <- DEMIMessages$DEMIExperiment$notACELFile( celpath, "celpath" );
				}
			}
		} else if ( length( celpath ) > 1 ) {
			for ( i in 1:length( celpath ) ) {
				# For some reason affxparser::isCelFile produces a warning messaage similar to:
				# 3: In readBin(con, what = integer(), size = 4, signed = FALSE,  ... :
  				# 'signed = FALSE' is only valid for integers of sizes 1 and 2
				# 4: closing unused connection 28 (/home/ilmjarv/work/projects/mef/data/CAHun-Ext-02-Control-2-MoEx-1_0-st-v1-a1.CEL)
				if ( affxparser::isCelFile( celpath[i] ) == FALSE ) {
					#msg <- paste( "The ", celpath[i], " in ", sQuote( "celpath" ), " is not a CEL file!\n", sep = "" );
					msg <- paste( msg, DEMIMessages$DEMIExperiment$notACELFile( celpath[i], "celpath" ), sep = "" );
				}
			}
		}
	}
	return( msg );
}#checkDEMIExperiment_celpath

#' Checks if \code{pmsize} is correct
#' 
#' Checks if \code{pmsize} is correct for the DEMI analysis. The 'pmsize' denotes perfect
#' match size meaning the size of continues perfect alignment between the probe and the
#' target sequence. It is used internally in DEMI analysis.
#' 
#' @param pmsize A \code{numeric}.
#' @return Returns NULL if \code{pmsize} is ok, else returns an error message.
#' @author Sten Ilmjarv
checkDEMIExperiment_pmsize <- function( pmsize )
{
	msg <- NULL;
	if ( is.numeric( pmsize ) == FALSE ) {
		#msg <- paste( sQuote( "pmsize" ), "needs to be a positive integer\n" );
		msg <- DEMIMessages$positiveInteger( "pmsize" );
		return( msg );
	}
	else if ( is.numeric( pmsize ) == TRUE ) {
		if ( pmsize < 0 || pmsize > 25 ) {
			#msg <- paste( sQuote( "pmsize" ), "needs to be a positive integer no bigger then 25\n" );
			msg <- DEMIMessages$positiveIntegerUpTo( "pmsize", 25 );
			return( msg );
		}
		if ( pmsize %% 1 != 0 ) {
			#msg <- paste( sQuote( "pmsize" ), "needs to be a positive integer no bigger then 25\n" );
			msg <- DEMIMessages$positiveIntegerUpTo( "pmsize", 25 );
			return( msg );
		}
	}
	
	return( msg );
}#checkDEMIExperiment_pmsize

#' Checks if \code{maxtargets} is correct
#' 
#' Checks if \code{maxtargets} is correct for the DEMI analysis. It is used internally in DEMI
#' analysis.
#' 
#' @param maxtargets A \code{numeric}.
#' @return Returns NULL if \code{maxtargets} is ok, else returns an error message.
#' @author Sten Ilmjarv
checkDEMIExperiment_maxtargets <- function( maxtargets )
{
	msg <- NULL;
	if ( is.numeric( maxtargets ) == FALSE ) {
		#msg <- paste( sQuote( "maxtargets" ), "needs to be a positive integer\n" );
		msg <- DEMIMessages$positiveInteger( "maxtargets" );
		return( msg );
	}
	else if ( is.numeric( maxtargets ) == TRUE ) {
		if ( maxtargets %% 1 != 0 ) {
			#msg <- paste( sQuote( "maxtargets" ), "needs to be a positive integer\n" );
			msg <- DEMIMessages$positiveInteger( "maxtargets" );
			return( msg );
		}
		if ( maxtargets < 0 ) {
			#msg <- paste( sQuote( "maxtargets" ), "needs to be a positive integer\n" );
			msg <- DEMIMessages$positiveInteger( "maxtargets" );
			return( msg );
		}
	}
	return( msg );
}#checkDEMIExperiment_maxtargets

#' Checks if \code{maxprobes} is correct
#' 
#' Checks if \code{maxprobes} is correct for the DEMI analysis. It is used internally in DEMI
#' analysis.
#' 
#' @param maxprobes A \code{numeric}.
#' @return Returns NULL if \code{maxprobes} is ok, else returns an error message.
#' @author Sten Ilmjarv
checkDEMIExperiment_maxprobes <- function( maxprobes )
{
	msg <- NULL;
	allGood = FALSE;

	if ( length( maxprobes ) > 0 ) {
		if ( maxprobes == "median" || maxprobes == "max" || is.null( maxprobes ) == TRUE ) {
			allGood = TRUE
		}
		else if ( as.numeric( maxprobes ) > 0 && as.numeric( maxprobes ) %% 1 == 0 ) {
			allGood = TRUE;
		}
	} else if ( length( maxprobes ) == 0 ) {
		return( msg );
	}
	if ( allGood == FALSE ) {
		#msg <- "The parameter 'maxprobes' has to be a positive integer or set as 'median', 'max' or not defined at all which by default set's it to 'max'.";
		msg <- DEMIMessages$DEMIExperiment$maxprobesError;
	}
	
	return( msg );
}#checkDEMIExperiment_maxprobes

#' Checks if \code{normormalization} is correct
#' 
#' Checks if \code{normalization} is correct for the DEMI analysis. 'normalization' stands for
#' normalization method or 'norm.method' in the \code{DEMIExperiment} object. It is used internally
#' in DEMI analysis.
#' 
#' @param normalization A \code{function}.
#' @return Returns NULL if \code{normalization} is ok, else returns an error message.
#' @seealso \code{DEMIExperiment}
#' @author Sten Ilmjarv
checkDEMIExperiment_normalization <- function( normalization )
{
	msg <- NULL;
	if ( is.function( normalization ) == FALSE ) {
		#msg <- paste( sQuote( "norm.method" ), "is not a proper function\n" );
		msg <- DEMIMessages$hasToBeFunction( "norm.method" );
	}
	return( msg );
}#checkDEMIExperiment_normalization


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#	DEMIExperiment draw functions:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @rdname probe.levels-methods
#' @aliases probe.levels,DEMIExperiment,character-method
#' @import methods
setMethod( "probe.levels", signature( object = "DEMIExperiment", target = "character" ),
	function( object, target )
	{
		if ( length( target ) == 1 ) {
			if( isTRUE( check4target( object, target ) ) ) {
				#library( ggplot2 );
				#library( reshape );
				probes <- getTargetProbes( object, target )$probeID;
				levels <- getProbeLevel( object, probes, FALSE );
				rownames( levels ) <- probes
				colnames( levels ) <- sub( ".CEL", "", colnames( levels ) );
				
				levels <- reshape::melt( levels );
				colnames( levels ) = c( "probeID", "Files", "value" );
				
                pic <- ggplot2::ggplot( data = levels, ggplot2::aes( Files, value ) ) + ggplot2::geom_boxplot( ggplot2::aes( fill = Files, colour = Files ) ) +
						ggplot2::scale_x_discrete("Probe IDs") + ggplot2::facet_wrap(~ probeID ) + ggplot2::scale_y_continuous("Relative rank values") +
						ggplot2::ggtitle( paste( paste( target, collapse = ",", sep = "" ), " gene probe levels" , sep="" ) ) +
						ggplot2::theme( axis.text.x = ggplot2::element_text(, angle = 310, hjust = 0, colour = "grey50") );
				
				return( pic );
			}
		} else {
			#stop( "Error: the target can only be of length 1\n" );
			stop( DEMIMessages$tooManyParameters( "target", 1 ) );
		}
	}
)

#' @rdname probe.plot-methods
#' @aliases probe.plot,DEMIExperiment,character-method
#' @import methods
setMethod( "probe.plot", signature( object = "DEMIExperiment", target = "character" ),
	function( object, target )
	{
		if ( length( target ) == 1 ) {
			if( isTRUE( check4target( object, target ) ) ) {
				#library( ggplot2 );
				#library( reshape );
				probes <- getTargetProbes( object, target )$probeID;
				levels <- getProbeLevel( object, probes, FALSE );
				
				rownames( levels ) <- probes
				colnames( levels ) <- sub( ".CEL", "", colnames( levels ) );
				
				#	Retrieve the annotation for the targets
				anno <- getTargetAnnotation( object, target );
				
				#	Calculate the width of the longest gene region
				end <- max( anno[ , "start" ] + anno[ , "length" ] );
				start <- min( anno[, "start" ] );
				
				levels <- reshape::melt( levels );
				colnames( levels ) = c( "probeID", "Files", "value" );
				
				pic <- ggplot2::ggplot( data = levels, ggplot2::aes( factor( probeID ), value ) )  + ggplot2::geom_point( ggplot2::aes( color = Files ), size = 4 ) +
						ggplot2::scale_x_discrete("Probe IDs") + ggplot2::scale_y_continuous("Relative rank values") +
						ggplot2::ggtitle( paste( paste( target, collapse = ",", sep = "" ), " gene probe plots" , sep="" ) ) +
                        ggplot2::theme( axis.text.x = ggplot2::element_text(, angle = 310, hjust = 0, colour = "grey50") );
				
				return( pic );
			}
		} else {
			#stop( "Error: the target can only be of length 1\n" );
			stop( DEMIMessages$tooManyParameters( "target", 1 ) );
		}
	}
)


## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
##	Old functions - loading data from file
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
#
##
##	Function loadAnnotationInfo
##
##setMethod( "loadAnnotationTable", signature( object = "DEMIExperiment" ),
##		function( object ) {
##			
##			cat ( "*Reading annotation files from '" );
##			
##			annoTable <- data.frame();
##			
##			if ( object@analysis == "transcript" || object@analysis == "gene" ) {
##				cat ( paste( demiConfig$taf, "/", object@organism, ".txt", "'", sep = "" ) );
##				annoTable <- read.table( paste( demiConfig$taf, "/", object@organism, ".txt", sep = "" ), header = TRUE, sep = "\t", comment.char = "#", fill = TRUE, quote = "" );
##				colnames( annoTable )[ grep( "ensembl_gene_id", colnames( annoTable ) ) ] = "geneID";
##				colnames( annoTable )[ grep( "ensembl_transcript_id", colnames( annoTable ) ) ] = "targetID";
##				colnames( annoTable )[ grep( "ensembl_peptide_id", colnames( annoTable ) ) ] = "peptideID";
##				colnames( annoTable )[ grep( "external_gene_id", colnames( annoTable ) ) ] = "geneName";
##				colnames( annoTable )[ grep( "transcript_biotype", colnames( annoTable ) ) ] = "biotype";
##			}
##			else if ( object@analysis == "exon" ) {
##				cat( paste( demiConfig$eaf, "/", object@organism, ".txt", "'", sep = "" ) );
##				annoTable <- read.table( paste( demiConfig$eaf, "/", object@organism, ".txt", sep = "" ), header = TRUE, sep = "\t", comment.char = "#", fill = TRUE, quote = ""  );
##				annoTable[ annoTable$strand == 1, grep( "strand", colnames( annoTable ) ) ] <- "+";
##				annoTable[ annoTable$strand == -1, grep( "strand", colnames( annoTable ) ) ] <- "-";
##				colnames( annoTable )[ grep( "ensembl_gene_id", colnames( annoTable ) ) ] = "geneID";
##				colnames( annoTable )[ grep( "ensembl_transcript_id", colnames( annoTable ) ) ] = "transcriptID";
##				colnames( annoTable )[ grep( "ensembl_exon_id", colnames( annoTable ) ) ] = "targetID";
##				colnames( annoTable )[ grep( "external_gene_id", colnames( annoTable ) ) ] = "geneName";
##			}
##			else if ( object@analysis == "genome" ) {
##				cat ( paste( demiConfig$gf, "/", object@organism, "/sectionAnnotations/", formatC( object@sectionsize, format = "d"), ".txt", "'", sep = "" ) );
##				annoTable <- read.table( paste( demiConfig$gf, "/", object@organism, "/sectionAnnotations/", formatC( object@sectionsize, format = "d"), ".txt", sep = "" ), header = FALSE, sep = "\t", comment.char = "#", fill = TRUE, quote = "" );
##				colnames( annoTable ) <- c( "targetID", "chr", "start", "end", "strand" );
##			}
##			
##			object@annoTable <- annoTable;
##			
##			cat( " ... Done\n" );
##			
##			return( object );
##		}
##)#loadAnnotationTable
##setMethod( "loadAnnotationTable", signature( object = "DEMIExperiment" ),
##		function( object ) {
##			data( DEMIannotation.data );
##		}
##)#loadAnnotationTable
#
##
##	Function readBlat
##
##setMethod( "loadBlatTable", signature( object = "DEMIExperiment" ),
##		function( object ) {
##			
##			cat( "*Loading blat table from '" );
##			
##			#	Define data
##			blatTable <- data.frame();
##			
##			colClasses <- c( rep( "character", 1 ), rep( "integer", 3 ), "character", rep( "integer", 2 ), "character" );
##			
##			#	Read in data
##			if ( object@analysis == "transcript" || object@analysis == "gene" ) {
##				cat ( paste( demiConfig$gpf, "/", object@arraytype, ".txt", "'", sep = "" ) );
##				blatTable <- read.table( paste( demiConfig$tpf, "/", object@arraytype, ".txt", sep = "" ), colClasses = colClasses, header = FALSE, sep = "\t", comment.char = "#" );
##			}
##			else if ( object@analysis == "exon" ) {
##				cat( paste( demiConfig$epf, "/", object@arraytype, ".txt", "'", sep = "" ) );
##				blatTable <- read.table( paste( demiConfig$epf, "/", object@arraytype, ".txt", sep = "" ), colClasses = colClasses, header = FALSE, sep = "\t", comment.char = "#"  );
##			}
##			else if ( object@analysis == "genome" ) {
##				cat ( paste( demiConfig$gf, "/", object@organism, "/sectionBlats/", object@arraytype, "/", formatC( object@sectionsize, format = "d"), ".txt", "'", sep = "" ) );
##				blatTable <- read.table( paste( demiConfig$gf, "/", object@organism, "/sectionBlats/", object@arraytype, "/blat/", formatC( object@sectionsize, format = "d"), ".txt", sep = "" ), colClasses = colClasses, header = FALSE, sep = "\t", comment.char = "#"  );
##			}
##			
##			#	Rename the data column names
##			colnames( blatTable ) = c( "probeID", "pmsize", "probeLength", "random", "targetID", "start", "end", "strand" );
##			
##			#	Remove those probes that match to a strand with less matches except for 'genome' analysis
##			#	since for genome it can match on both strands
##			if ( object@analysis != "genome" ) {
##				strand <- NULL;
##				strand_table = table( blatTable$strand );
##				if ( strand_table[ names( strand_table ) == "+" ] > strand_table[ names( strand_table ) == "-" ] ) {
##					strand <- "+";
##					cat( "\n\tWill ignore the '-' strand matches" );
##				} else {
##					strand <- "-";
##					cat( "\n\tWill ignore the '+' strand matches" );
##				}
##				blatTable <- blatTable[ blatTable$strand == strand, ];
##			}
##			
##			#	Remove probes whose perfect match is less then set by the user
##			if ( is.null( object@pmsize ) == FALSE ) {
##				blatTable <- blatTable[ blatTable$pmsize >= object@pmsize, ];
##			}
##			
##			#	Remove probes who match to more targets then set by the user
##			if ( is.null( object@maxtargets ) == FALSE & is.null( object@annoTable ) == FALSE ) {
##				cat( "\n\tDetermining the number of matches on targets for every probe" );
##				if ( object@analysis == "genome") {
##					#	Read in the maxtarget file
##					probeTargetCount <- read.table( paste( demiConfig$gf, "/", object@organism, "/sectionBlats/", object@arraytype, "/maxtarget/", object@pmsize, "/", formatC( object@sectionsize, format = "d" ), ".txt", sep = "" ), header = FALSE, sep = "\t", comment.char = "#" );
##					colnames( probeTargetCount ) <- c( "probeID", "countTargets" );
##					outProbes <- as.vector( probeTargetCount[ probeTargetCount$countTargets > object@maxtargets, c( "probeID" ) ] );
##					blatTable <- blatTable[ !( blatTable$probeID %in% outProbes ), ];
##				} else {
##					#	If analysis is 'transcript', 'exon' or 'gene'
##					tempBlat <- merge( blatTable[, c( "probeID", "targetID" )], data.frame( object@annoTable[, c( "geneID", "targetID" )] ), by = "targetID" );
##					tempBlat <- data.frame( tempBlat[, c("probeID", "geneID")] );
##					#	This fixes the problem that probes that match to one target several times are not counted
##					#	as having several matches. Instead they are counted only on having one match.
##					tempBlat <- unique( tempBlat );
##					probeTargetCount <- data.frame( table( tempBlat$probeID ) );
##					outProbes <- as.vector( probeTargetCount[ probeTargetCount$Freq > object@maxtargets, c("Var1") ] );
##					if ( object@analysis == "transcript" || object@analysis == "exon" ) {
##						blatTable <- blatTable[ !( blatTable$probeID %in% outProbes ), ];
##					} else if ( object@analysis == "gene" ) {
##						tempBlat <- tempBlat[ !( tempBlat$probeID %in% outProbes ), ];
##						colnames( tempBlat )[ grep( "geneID", colnames( tempBlat ) ) ] = "targetID";
##						blatTable <- tempBlat;
##					}
##				}
##			}
##			
##			object@blatTable <- blatTable;
##			
##			cat( "\n\t... Done\n" );
##			
##			return( object );
##			
##		}
##)#loadBlatTable
#
##
##	Function loadCytoband
##
##setMethod( "loadCytoband", signature( object = "DEMIExperiment" ),
##		function( object ) {
##			
##			cat ( "*Loading cytoband from '" );
##			
##			cytoband <- getCytoband( object );
##			cytofile <- paste( demiConfig$cbf, "/", getOrganism( object ), ".txt", sep = "" );
##			if ( file.exists( cytofile ) == TRUE ) {	
##				cytoband <- read.table( file = cytofile, sep = "\t", comment.char = "#", header = FALSE );
##				colnames( cytoband ) <- c( "chr", "start", "end", "region", "method" );
##				#cytoband$chr <- sub( "chr", "", cytoband$chr ); # [ Unnecessary because the karyotype data from ensembl comes already without the 'chr' prefix ]
##			}
##			object@cytoband <- cytoband;
##			return( object );
##		}
##)#loadCytoband
