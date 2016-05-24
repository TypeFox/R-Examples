#==============================================================================#
# Constructors.R: functions serving as class contructors
# Note: suggested by Martin Morgan to save the user from calling new()
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# DEMIExperiment
# DEMICel
# DEMIGroup
# DEMIClust
# DEMIDiff
# demi
#==============================================================================#

#' Creates a \code{DEMIExperiment} object
#' 
#' This function creates a \code{DEMIExperiment} object. It loads and stores the experiment metadata such
#' as annotation and alignment information and raw expression matrix from CEL files. It then normalizes
#' the raw expression matrix and stores both expression matrices in a \code{DEMICel} object stored under the
#' created \code{DEMIExperiment} object.
#' 
#' @param analysis A \code{character}. Defines the analysis type. It can be either 'transcript',
#' 		  'gene', 'exon' or 'genome'. The default value is 'transcript'. For 'genome' analysis
#' 		  \code{sectionsize} parameter needs to be defined as well.
#' @param celpath A \code{character}. It can point to the directory containing CEL files or is a vector
#' 		  that points directly to the CEL files.
#' @param experiment A \code{character}. A custom name of the experiment defined by the user (e.g. 'myexperiment').
#' @param organism A \code{character}. The name of the species the microarrays are measuring (e.g. 'homo_sapiens'
#' 		  or 'mus_musculus') given in lowercase letters and words are separated by underscore.
#' @param maxtargets A \code{numeric}. The maximum number of allowed targets (e.g. genes or transcripts) one probe
#' 		  can have a match against. If to set it to 1 it means that the probe can match only one gene. If the \code{analysis}
#' 		  is set to 'transcript' the program still calculates the number of matches on genes, not transcripts. Hence
#' 		  a probe matching two transcripts on the same gene would be included but a probe matching two
#' 		  transcripts on different genes would not be included. The value needs to be a positive integer or 0.
#' 		  By default \code{maxtargets} is set to 0.
#' @param maxprobes A \code{character}. Sets the number of unique probes a target is allowed to have a match against. All
#' 		  the targets that yield more alignments to different probes then set by \code{maxprobes} will be scaled
#' 		  down to the number defined by the \code{maxprobes} parameter. It can be either a positive integer or set as
#' 		  'median' or 'max' - 'median' meaning the median number of probes matching to all targets and 'max'
#' 		  meaning the maximum number of probes matching to a target. By default \code{maxprobes} is not set which is
#' 		  the same as setting \code{maxprobes} to 'max'.
#' @param pmsize A \code{numeric}. The minimum number of consecutive nucleotides that need to match perfectly
#' 		  against the target sequence. It can be either 23, 24 or 25. This means that alignments with
#' 		  smaller perfect match size will not be included in the experiment set up. The default value is 25. 
#' @param sectionsize A \code{numeric}. This is only used if the \code{analysis} parameter is set to 'genome'. It defines the
#' 		  length of the genomic target region used in the 'genome' analysis. Currently the only available section
#' 		  sizes are 100000, 500000 and 1000000.
#' @param norm.method A \code{function}. Defines a function used to normalize the raw expression values. The default normalization
#' 		  function is \code{norm.rank}.
#' @param filetag A \code{character}. This is a custom string that can be used to identify the experiment. At the current
#' 		  development stage this parameter is used only when using the function \code{demi}, where the output files will
#' 		  contain the specified filetag.
#' @return A \code{DEMIExperiment} object.
#' @seealso \code{DEMIClust}, \code{DEMIResult}, \code{getResultTable}, \code{getResult}, \code{attachResult}
#' @details
#' 
#' After the analysis has been completed the user can add the results from the analysis to the original \code{DEMIExperiment} object with
#' the function \code{attachResult}. Then the function \code{getResultTable} can be used to retrieve the results from the
#' \code{DEMIExperiment} object. Other useful functions are \code{getNormMatrix} to retrieve normalized expression matrix and
#' \code{getCelMatrix} to retrieve the raw expression matrix. In both cases the probe ID's are present as row names.
#' 
#' Further specification of the parameters:
#' \itemize{
#' 	\item{maxtargets}{
#' 		When \code{analysis} is set to 'gene' then all probes that match to more genes then allowed by \code{maxtargets} parameter will
#' 		not be included in the analysis. For 'transcript' and 'exon' analysis the number is also calculated on a gene
#' 		level. For example if \code{maxtargets} is set to one and a probe matches to two transcripts but on the same gene,
#' 		then this probe will still be used in the analysis. However if the probe matches two transcripts on different
#' 		genes then this probe will not be included in the analysis. For 'genome' analysis the probe in most cases matches
#' 		to two genomic sections because adjacent sections overlap by 50%. However this is considered as one match and the
#' 		probe will still be used in the analysis.
#' 		}
#' 	\item{norm.method}{
#' 		Every user can apply their own normalization method by writing a custom normalization function. The function should
#' 		take in raw expression matrix and return the normalized expression matrix where probe ID's are kept as rownames and
#' 		column names are CEL file names. The normalized expression matrix will then be stored as part of the \code{DEMIExperiment}
#' 		object.
#' 		}
#' 	\item{sectionsize}{
#' 		The \code{sectionsize} parameter defines the length of the genomic target region. Currenlty \code{sectionsize} can be set
#' 		as: 100000, 500000 and 1000000. All adjacent sections, except the ones on chromosome ends, overlap with the
#' 		next adjacent section by 50%. It ensures the all probes matching to genome will be assigned to at least one
#' 		genomic section. This parameter is required when \code{analysis} is set to 'genome'.
#' 		}
#' 	\item{norm.method}{
#' 		The \code{norm.method} defines a function to use for the normalization of raw expression matrix. The user can implement his/her
#' 		own function for the normalization procedure. The function should take in raw expression matrix and return the normalized
#' 		expression matrix where probe ID's are kept as rownames and column names are CEL file names.
#' 		}
#' }
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
#' # Now we can continue the example of the function DEMIExperiment
#' 
#' # Basic experiment set up.
#' demiexp <- DEMIExperiment(analysis = 'gene', celpath = destfolder,
#' 		experiment = 'myexperiment', organism = 'homo_sapiens')
#' 
#' # Run basic experiment set up but this time do 'transcript' analysis.
#' demiexp <- DEMIExperiment(analysis = 'transcript', celpath = destfolder,
#' 		experiment = 'myexperiment', organism = 'homo_sapiens')
#' 
#' # Run basic experiment set up but this time do 'transcript' analysis.
#' demiexp <- DEMIExperiment(analysis = 'exon', celpath = destfolder,
#' 		experiment = 'myexperiment', organism = 'homo_sapiens' )
#' 
#' # For genome analysis do not forget to specify the sectionsize parameter.
#' demiexp <- DEMIExperiment(analysis = 'genome', celpath = destfolder,
#' 		experiment = 'myexperiment', organism = 'homo_sapiens', sectionsize = 500000)
#' 
#' # Specify experiment with specific pmsize; the standard length for Affymetrix microarray
#' # probes is 25 nucleotides.
#' demiexp <- DEMIExperiment(analysis = 'gene', celpath = destfolder,
#' 		experiment = 'myexperiment', organism = 'homo_sapiens', pmsize = 23)
#'
#' # Specify experiment by setting maxtargets to 1.
#' demiexp <- DEMIExperiment(analysis = 'gene', celpath = destfolder,
#' 		experiment = 'myexperiment', organism = 'homo_sapiens', maxtargets = 1)
#'  
#' # Specify experiment by setting maxprobes to 'median'.
#' demiexp <- DEMIExperiment(analysis = 'gene', celpath = destfolder,
#' 		experiment = 'myexperiment', organism = 'homo_sapiens', maxprobes = 'median')
#' 
#' # Retrieve the alignment information from the DEMIExperiment object.
#' head( getAlignment( demiexp ) )
#' 
#' # Retrieve the annotation information from the DEMIExperiment object.
#' head( getAnnotation( demiexp ) )
#' 
#' # Retrieve the raw expression matrix from the DEMIExperiment object.
#' head( getCelMatrix( demiexp ) )
#' 
#' # Retrieve the normalized expression matrix from the DEMIExperiment object.
#' head( getNormMatrix( demiexp ) )
#' 
#' #####################
#' # If the user has done the analysis and wishes to add the results to the original
#' # DEMIExperiment object.
#' #####################
#' 
#' # Create clusters with an optimized wilcoxon's rank sum test incorporated within demi that
#' # precalculates the probabilities.
#' demiclust <- DEMIClust( demiexp, group = c( "BRAIN", "UHR" ), clust.method = demi.wilcox.test.fast )
#' # Calcuate differential expression
#' demidiff <- DEMIDiff( demiclust )
#' 
#' # Attach the results to the original DEMIExperiment object
#' demiexp <- attachResult( demiexp, demidiff )
#' 
#' # Retrieve the results from the DEMIExperiment object
#' head( getResultTable( demiexp ) )
#' 
#' }
#' 
#' @export
#' @docType methods
#' @rdname DEMIExperiment-methods
#' @import methods
"DEMIExperiment" <-
function( analysis 				= "transcript",
		  celpath				= character(),
		  experiment			= character(),
		  organism				= character(),
		  maxtargets			= 0,
		  maxprobes				= character(),
		  pmsize				= 25,
		  sectionsize			= character(),
		  norm.method			= norm.rrank,
		  filetag 				= character() )
{
	
	#	parameter validation
	if ( missing( analysis ) == TRUE ) {
		#stop( paste( "Error:", sQuote( "analysis" ), "is unspecified" ) );
		stop( DEMIMessages$parameterMissing( "analysis" ) );
	} else if ( length( analysis ) > 1 ) {
		#stop( paste( "Error:", sQuote( "analysis" ), "can only be of length 1" ) );
		stop( DEMIMessages$tooManyParameters( "analysis", 1 ) );
	} else if ( analysis == "genome" ) {
		if ( missing( sectionsize ) == TRUE ) {
			#stop( paste( "Error:", sQuote( "sectionsize" ), "has not been specified for 'genome' analysis" ) );
			stop( DEMIMessages$sectionsizeMissing() );
		}
	}
	if ( missing( organism ) == TRUE ) {
		#stop( paste( "Error:", sQuote( "organism" ), "is unspecified" ) );
		stop( DEMIMessages$parameterMissing( "organism" ) );
	} else if ( length( organism ) > 1 ) {
		#stop( paste( "Error:", sQuote( "organism" ), "can only be of length 1" ) );
		stop( DEMIMessages$tooManyParameters( "organism", 1 ) );
	}
	if ( missing( celpath ) == TRUE ) {
		#stop( paste( "Error:", sQuote( "celpath" ), "is unspecified" ) );
		stop( DEMIMessages$parameterMissing( "celpath" ) );
	}
	if ( missing( experiment ) == TRUE ) {
		#stop( paste( "Error:", sQuote( "experiment" ), "is unspecified" ) );
		stop( DEMIMessages$parameterMissing( "experiment" ) );
	} else if ( length( experiment ) > 1 ) {
		#stop( paste( "Error:", sQuote( "experiment" ), "can only be of length 1" ) );
		stop( DEMIMessages$tooManyParameters( "experiment", 1 ) );
	}
	if ( length( maxtargets ) > 0 && is.numeric( maxtargets ) == FALSE ) {
		#stop( paste( "Error:", sQuote( "maxtargets" ), "has to be 0 or a positive integer" ) );
		stop( DEMIMessages$positiveIntegerOrZero( "maxtargets" ) );
	} else if ( length( maxtargets ) > 1 ) {
		#stop( paste( "Error:", sQuote( "maxtargets" ), "can only be of length 1" ) );
		stop( DEMIMessages$tooManyParameters( "maxtargets", 1 ) );
	}
	if ( length( pmsize ) > 0 && is.numeric( pmsize ) == FALSE ) {
		#stop( paste( "Error:", sQuote( "pmsize" ), "has to be a positive integer up to 25" ) );
		stop( DEMIMessages$positiveIntegerUpTo( "pmsize", 25 ) );
	} else if ( length( pmsize ) > 1 ) {
		#stop( paste( "Error:", sQuote( "pmsize" ), "can only be of length 1" ) );
		stop( DEMIMessages$tooManyParameters( "pmsize", 1 ) );
	}
	if ( length( sectionsize ) > 0 && is.numeric( sectionsize ) == FALSE ) {
		#stop( paste( "Error:", sQuote( "sectionsize" ), "has to be a positive integer" ) );
		stop( DEMIMessages$positiveInteger( "sectionsize" ) );
	} else if ( length( sectionsize ) > 1 ) {
		#stop( paste( "Error:", sQuote( "sectionsize" ), "can only be of length 1" ) );
		stop( DEMIMessages$tooManyParameters( "sectionsize", 1 ) );
	}
	if ( is.function( norm.method ) == FALSE ) {
		#stop( paste( "Error:", sQuote( "norm.method" ), "has to be a function" ) );
		stop( DEMIMessages$hasToBeFunction( "norm.method" ) );
	} else if ( length( norm.method ) > 1 ) {
		#stop( paste( "Error:", sQuote( "norm.method" ), "can only be of length 1" ) );
		stop( DEMIMessages$tooManyParameters( "norm.method", 1 ) );
	}
	if ( missing( sectionsize ) == FALSE && analysis != "genome" ) {
		#stop( paste( "Error:", sQuote( "sectionsize" ), "can only be used in 'genome' analysis\n" ) );
		stop( DEMIMessages$sectionsizeUse() )
	}
	
	#	exprsData - will be generated later, during validation
	experiment <- new( "DEMIExperiment",
						analysis			= analysis,
						celpath				= celpath,
						experiment 			= experiment,
						organism 			= organism,
						maxtargets			= as.numeric( maxtargets ),
						maxprobes 			= as.character( maxprobes ),
						pmsize				= as.numeric( pmsize ),
						sectionsize			= as.numeric( sectionsize ),
						norm.method			= as.function( norm.method ),
						filetag 			= filetag
					);
	
	return( experiment );
	
}#DEMIExperiment constructor


#' Creates a \code{DEMICel} object
#' 
#' A \code{DEMICel} holds the raw and normalized expression matrices. It is used
#' internally in DEMI analysis.
#' 
#' @param celMatrix A \code{matrix}. The raw expression matrix.
#' @param normMatrix A \code{matrix}. The normalized expression matrix.
#' @return A \code{DEMICel} object that holds the raw and normalized expression
#' 		   matrices.
#' @details
#' 
#' Both expression matrices store the expression values in their columns
#' and the column name represents the original CEL file name. The row names
#' represent probe ID's which makes it easy to retrieve specific expression
#' data for specific probes.
#' 
#' @author Sten Ilmjarv
#' 
#' @export
#' @docType methods
#' @rdname DEMICel-methods
#' @import methods
"DEMICel" <-
function( celMatrix		= matrix(),
		  normMatrix	= matrix())
{
	celData <- new( "DEMICel",
					celMatrix	= celMatrix,
					normMatrix	= normMatrix
			);
			
	return( celData );
}#DEMICel constructor


#' Creates a \code{DEMIGroup} object
#' 
#' A \code{DEMIGroup} object holds the group annotations such as
#' the column indexes of both groups and names of the groups. It is used
#' internally in DEMI analysis.
#'
#' @param groupA A \code{character}. Holds the name of group A.
#' @param groupB A \code{character}. Holds the name of group B.
#' @param indexA A \code{numeric}. A vector of column indexes belonging to group A.
#' @param indexB A \code{vector}. A vector of column indexes belonging to group B.
#' @param groupNames A \code{character}. Holds the names of custom groups created by
#' 		  the user.
#' @return A \code{DEMIGroup} object that holds the group annotations.
#' @details 
#' 
#' The \code{DEMIGroup} can hold both automatically generated annotations that
#' depend on the group names or custom annotations specified by the user.
#' The automatically generated ones are created by scanning for the specified
#' group names in the column names of the normalized expression matrix. It then
#' retrieves the column indexes where the specified group names occure. The
#' custom group names are just stored in the \code{groupNames} vector and all the other
#' parameters of the \code{DEMIGroup} object will be left empty.
#' 
#' @author Sten Ilmjarv
#' 
#' @export 
#' @docType methods
#' @rdname DEMIGroup-methods
#' @import methods
"DEMIGroup" <-
function( groupA		= character(),
		  groupB		= character(),
		  indexA		= numeric(),
		  indexB		= numeric(),
		  groupNames	= character() )
{
	demiGroups <- new( "DEMIGroup",
					   groupA		= groupA,
					   groupB		= groupB,
					   indexA		= indexA,
					   indexB		= indexB,
					   groupNames	= groupNames
			   );
	return( demiGroups );
}#DEMIGroups constructor


#' Creates a \code{DEMIClust} object
#' 
#' A \code{DEMIClust} object clusters probes by their expression profile. The clustering
#' is done with a function defined by the \code{clust.method} parameter. One could also define custom
#' clusters by defining the \code{cluster} parameter with a list of probes. It then stores the
#' clusters of probes as a \code{DEMIClust} object.
#'
#' @param experiment A \code{DEMIExperiment} object. Holds the \code{DEMIExperiment} object whose metadata
#' 		  (such as normalized expression values) is used to cluster the probes.
#' @param group A \code{character}. Defines the groups that are used for clustering (e.g 'group = c("TEST", "CONTROL")').
#' 		  It uses \code{grep} function to locate the group names from the CEL file names and then builds
#' 		  index vectors determining which files belong to which groups.
#' @param clust.method A \code{function}. Defines the function used for clustering. The user can
#' 		  build a custom clustering function. The input of the custom function needs
#' 		  to be the same \code{DEMIClust} object and the output is a \code{list} of probes, where
#' 		  each list corresponds to a specific cluster. The default function is \code{demi.wilcox.test}
#' 		  that implements the \code{wilcox.test} function. However we recommend to use the function
#' 		  \code{demi.wilcox.test.fast} that uses a custom \code{wilcox.test} and runs a lot faster.
#' @param cluster A \code{list}. Holds the probes of different clusters in a \code{list}.
#' @param cutoff.pvalue A \code{numeric}. Sets the cut-off p-value used for determining statistical
#' 		  significance of the probes when clustering the probes into clusters. Default is 0.05.
#' @return A \code{DEMIClust} object.
#' @seealso \code{DEMIExperiment}, \code{demi.wilcox.test}, \code{demi.wilcox.test.fast}, \code{demi.comp.test},
#' 			\code{wprob}
#' @details
#' 
#' Instead of automatically clustered probes \code{DEMIClust} object can use user defined lists of
#' probes for later calculation of differential expression. This is done by setting the \code{cluster} parameter. It
#' overrides the default behaviour of the \code{DEMIClust} object and no actual clustering occurs. Instead the list of probes defined
#' in the \code{cluster} parameter are considered as already clustered probes. The list needs to contain proper names for
#' probe vectors so that they would be recognizable later. Also instead of using the default clustering method the user can
#' write his/her own function for clustering probes based on the expression values.
#' 
#' Further specification of the parameters:
#' \itemize{
#' 	\item{group}{
#' 		All the CEL files used in the analysis need to contain at least one of the names specified in the \code{group}
#' 		parameter because they determine what groups to compare against each other. It is also a good
#' 		practice to name the CEL files to include their common features. However if a situation arises where the
#' 		group/feature name occurs in all filenames then the user can set group names with specific filenames
#' 		by seperating names in one group with the "|" symbol. For example \code{group = c( "FILENAME1|FILENAME2|FILENAME3",
#' 		"FILENAME4|FILENAME5|FILENAME6" )}. These two groups are then used for clustering the probes expression values.
#' 		}
#' 	\item{clust.method}{
#' 		The user can write his/her own function for clustering probes according to their expression values. The custom
#' 		function should take \code{DEMIClust} object as the only parameter and output a \code{list}. The output list should contain
#' 		the name of the clusters and the corresponding probe ID's. For example \code{return( list( cluster1 = c(1:10), cluster2
#' 		= c(11:20), cluster3 = c(21:30) )}.
#' 		}
#' 	\item{cluster}{
#' 		This parameter allows to calculate differential expression on user defined clusters of probe ID's. It needs to be a
#' 		\code{list} of probe ID's where the \code{list} names correspond to the cluster names. For example \code{list( cluster1 = c(1:10),
#' 		cluster2(1:10) )}. When using this approach you need to make sure that all the probe ID's given in the clusters are available in
#' 		the analysis. Otherwise an error message will be produced and you need to remove those probes that have no alignment
#' 		in the analysis. When setting this parameter the default behaviour will be overridden and no default clustering
#' 		will be applied.
#' 		}
#' }
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
#' # Now we can continue the example of the function DEMIClust
#' 
#' # Set up an experiment.
#' demiexp <- DEMIExperiment(analysis = 'gene', celpath = destfolder,
#' 			experiment = 'myexperiment', organism = 'homo_sapiens')
#' 
#' # Create clusters with default behaviour
#' demiclust <- DEMIClust( demiexp, group = c( "BRAIN", "UHR" ) )
#' 
#' # Create clusters with an optimized wilcoxon's rank sum test incorporated within demi that
#' # precalculates the probabilities.
#' # The user can specify his/her own function for clustering.
#' demiclust <- DEMIClust( demiexp, group = c( "BRAIN", "UHR" ), clust.method = demi.wilcox.test.fast )
#' 
#' # Create a 'DEMIClust' object with custom lists of probeID's
#' demiclust <- DEMIClust( demiexp, cluster = list( customcluster = c(1190, 1998, 2007) ) )
#' 
#' # To retrieve the clusters use
#' getCluster( demiclust )
#' 
#' # To retrieve cluster names use
#' names( getCluster( demiclust ) )
#' 
#' }
#' 
#' @export 
#' @docType methods
#' @rdname DEMIClust-methods
#' @import methods
"DEMIClust" <-
function( experiment	= "DEMIExperiment",
		  group			= character(),
		  clust.method	= function(){},
		  cluster		= list(),
		  cutoff.pvalue	= 0.05 )
{
	
	#	check if experiment has been set
	if ( missing( experiment ) == TRUE ) {
		#stop( paste( "Error:", sQuote( "experiment" ), "of class 'DEMIExperiment' is unspecified\n" ) );
		stop( DEMIMessages$parameterOfClassMissing( "experiment", "DEMIExperiment" ) );
	}
	
	#	check if 'cutoff.pvalue' has been set
	if ( missing( cutoff.pvalue ) == FALSE ) {
		if ( is.numeric( cutoff.pvalue ) == FALSE ) {
			#stop( paste( "Error:", sQuote( "cutoff.pvalue" ), "is not of class 'numeric'\n" ) );
			stop( DEMIMessages$parameterNotOfClass( "cutoff.pvalue", "numeric" ) );
		}
	}
	
	demiGroup <- NULL;
	demiClusters <- list();
	#	if the user overrides the clustering by adding custom probe vector to DEMIClust object
	if ( missing( cluster ) == FALSE ) {
		#	if cluster name has been set
#		if ( missing( group ) == TRUE ) {
#			#stop( paste( "Error:", sQuote( "group" ), "is unspecified. You need to specify the name of you probe cluster\n" ) );
#			stop( DEMIMessages$parameterMissing( "group" ) );
#		} else if ( missing( group ) == FALSE ) {
#			#	if the number of clusters and the names are same
#			if ( length( cluster ) != length( group ) ) {
#				#stop( paste( "Error:", sQuote( "group" ), "and", sQuote( "cluster" ), "are of different length\n" ) );
#				stop( DEMIMessages$paramAndParamNotEqualLength( "group", "cluster" ) );
#			} else {
#				#	check if all the in the cluster are available in the blatTable
#				for ( i in 1:length( cluster ) ) {
#					probesOK = check4probe( experiment, cluster[[i]] );
#					if ( is.null( probesOK ) == FALSE ) {
#						stop( probesOK );
#					}
#				}
#				#	since all the probes were present in the blat table continue
#				for ( i in 1:length( cluster ) ) {
#					demiClusters[i] <- cluster[i];
#				}
#				names( demiClusters ) <- group;
#				demiGroup <- DEMIGroup( groupNames = group );
#			}
#		}
		#	check if all the probes in the clusters are available in the blatTable
		for ( i in 1:length( cluster ) ) {
			probesOK = check4probe( experiment, cluster[[i]] );
			if ( is.null( probesOK ) == FALSE ) {
				stop( probesOK );
			}
		}
		#	since all the probes were present in the blat table continue
#		for ( i in 1:length( cluster ) ) {
#			demiClusters[names( cluster[i] )] <- cluster[i];
#		}
		# no need for the above because cluster already is in the correct format
		demiClusters <- cluster;
		demiGroup <- DEMIGroup( groupNames = names( cluster ) );
	} else if( missing( cluster ) == TRUE ) { #	if the user uses default parameters
		#	check if group's has been set
		if ( missing( group ) == TRUE ) {
			stop( DEMIMessages$parameterOfClassMissing( "group", "character" ) );
		} else if ( length( group ) != 2 ) {
			stop( DEMIMessages$tooManyParameters( "group", 2 ) );
		}
		#	check if custom clustering method has been set
		if ( missing( clust.method ) == TRUE ) {
			#	set demi.wilcox.test as the default clustering method if it is missing
			clust.method = demi.wilcox.test;
		}
		#	check if clust.method is a method
		if ( is.function( clust.method ) == FALSE ) {
			#stop( paste( "Error:", sQuote( "clust.method" ), "has to be a function\n" ) );
			stop( DEMIMessages$hasToBeFunction( "clust.method" ) );
		}
		#	check for correct definitions
		if ( is.vector( group ) == FALSE || length( group ) != 2 ) {
			#stop( paste( "Error:\tThere can only be two character elements in the", sQuote( "group" ), "vector\n",
			#				"\t\tAll group names need to be present in CEL-file names in your CEL-file directory specified by", sQuote( "celpath"), "parameter\n",
			#				"\t\tYou can express one group with regular expression by seperating groups in a string with a '|' symbol, for we use", sQuote( "grep" ), "to find group indexes from CEL file names\n",
			#				"\t\ti.e. '> groups = c( \"TUMOR\", \"NORMAL\" )'\n")
			#);
			stop( DEMIMessages$wrongDefinition() );
		} else {
			demiGroup <- DEMIGroup( groupA = group[1], groupB = group[2] );
		}
	}
	
	demiClust <- new( "DEMIClust",
					  experiment 	= experiment,
					  group			= demiGroup,
					  clust.method	= clust.method,
					  cluster		= demiClusters,
					  cutoff.pvalue	= cutoff.pvalue
			  );
	
	return( demiClust );
	
}#DEMIClust constructor


#' Creates a \code{DEMIDiff} object
#' 
#' The \code{DEMIDiff} object calculates differential expression and holds the analysis results,
#' results of clustering and the original metadata of the experiment. To retrieve the results
#' from the \code{DEMIDiff} object use the the function \code{getResultsTable} that returns the results as
#' a \code{data.frame}.
#' 
#' @param cluster A \code{DEMIClust} object. The \code{DEMIClust} object that holds the clusters used in the
#' 		  analysis.
#' @return A \code{DEMIDiff} object.
#' @seealso \code{DEMIExperiment}, \code{DEMIClust}, \code{DEMIResult}, \code{getResultTable}, \code{getResult}, \code{attachResult}
#' @details
#' 
#' The \code{DEMIDiff} object calculates the differential expression for every cluster in the \code{DEMIClust}
#' object set by the \code{cluster} parameter. The results are then stored in the \code{DEMIDiff} object under
#' the slot \code{result} as a \code{DEMIResult} object. This object can be retrieved with the function \code{getResult}
#' but most of the times it is recommended to use the function \code{getResultTable} which returns the results
#' in a \code{data.frame} sorted by the FDR values.
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
#' # Now we can continue the example of the function DEMIDiff
#' 
#' # Set up an experiment.
#' demiexp <- DEMIExperiment(analysis = 'gene', celpath = destfolder,
#' 		experiment = 'myexperiment', organism = 'homo_sapiens')
#' 
#' # Create clusters with an optimized wilcoxon's rank sum test incorporated within demi that
#' # precalculates the probabilities.
#' demiclust <- DEMIClust( demiexp, group = c( "BRAIN", "UHR" ), clust.method = demi.wilcox.test.fast )
#' 
#' # Calcuate differential expression
#' demidiff <- DEMIDiff( demiclust )
#' 
#' # Retrieve the results in a 'data.frame'
#' head( getResultTable( demidiff ) )
#' 
#' # Attach the results to the original 'DEMIExperiment' object
#' demiexp <- attachResult( demiexp, demidiff )
#' 
#' # Retrieve the results from the 'DEMIExperiment' object
#' head( getResultTable( demiexp ) )
#' 
#' }
#' 
#' @export 
#' @docType methods
#' @rdname DEMIDiff-methods
#' @import methods
"DEMIDiff" <-
function( cluster = "DEMIClust" )
{
	
	# check if clusters are defined
	if ( missing( cluster ) == TRUE ) {
		#stop( paste( "Error:", sQuote( "cluster" ), "of class 'DEMIClust' is unspecified\n" ) );
		stop( DEMIMessages$parameterOfClassMissing( "cluster", "DEMIClust" ) );
	}
	
	diffName <- paste( cluster@group@groupA, "_VS_",cluster@group@groupB, sep = "" );
	
	demiDiff <- new( "DEMIDiff",
					 cluster			= cluster,
					 name				= diffName
			 );
			 
	return( demiDiff );
	
}#DEMIDiff constructor


#' A wrapper for DEMI analysis
#' 
#' Function \code{demi} is a wrapper for the whole DEMI analysis. First it creates a \code{DEMIExperiment} object, then uses
#' it to create a \code{DEMIClust} object that contains the list of clustered probes and then performs differential
#' expression analysis by running the function \code{DEMIDiff} that creates \code{DEMIDiff} object. The latter contains
#' the results of the differential expression analysis. It also prints out the results to the working directory.
#' If parameter \code{pathway} is set to TRUE, it also performs gene ontology analysis on the results in \code{DEMIDiff}
#' object to determine statistically significant gene ontology categories (it also prints out those in the working
#' directory with the file containing the string 'pathway'). It then returns a list containing the \code{DEMIExperiment}
#' object where the results have been attached to and a \code{data.frame} that contains the functional annotation
#' analysis results. NB! The results will be printed out in the working directory.
#' 
#' @param analysis A \code{character}. Defines the analysis type. It can be either 'transcript',
#' 		  'gene', 'exon' or 'genome'. The default value is 'transcript'. For 'genome' analysis
#' 		  \code{sectionsize} parameter needs to be defined as well.
#' @param celpath A \code{character}. It can point to the directory containing CEL files or is a vector
#' 		  that points directly to the CEL files.
#' @param experiment A \code{character}. A custom name of the experiment defined by the user (e.g. 'myexperiment').
#' @param organism A \code{character}. The name of the species the micrroarrays are measuring (e.g. 'homo_sapiens'
#' 		  or 'mus_musculus') given in lowercase and words separated by underscore.
#' @param maxtargets A \code{numeric}. The maximum number of allowed targets (e.g. genes or transcripts) one probe
#' 		  can match against. If to set it to 1 it means that the probe can match only one gene. If the \code{analysis}
#' 		  is set to 'transcript' the program still calculates the number of matches on genes. Hence a probe
#' 		  matching two transcripts on the same gene would be included but a probe matching two transcripts
#' 		  on different genes would not be included. The value needs to be a positive integer or 0. By default
#' 		  \code{maxtargets} is set to 0.
#' @param maxprobes A \code{character}. Sets the number of unique probes a target is allowed to have a match against.
#' 		  All the targets that yield more alignments to different probes then set by \code{maxprobes} will be scaled
#' 		  down to the number defined by the \code{maxprobes} parameter. It can be either a positive integer or set as
#' 		  'median' or 'max' - 'median' meaning the median number of probes matching to all targets and 'max'
#' 		  meaning the maximum number of probes matching to a target. By default \code{maxprobes} is not set which is
#' 		  the same as setting \code{maxprobes} to 'max'.
#' @param pmsize A \code{numeric}. The minimum number of consecutive nucleotides that need to match perfectly
#' 		  against the target sequence. It can be either 23, 24 or 25. This means that alignments with
#' 		  smaller perfect match size will not be included in the experiment set up. The default value is 25.
#' @param sectionsize A \code{numeric}. This is only used if the \code{analysis} parameter is set to 'genome'. It defines the
#' 		  length of the genomic target region used in the 'genome' analysis.
#' @param group A \code{character}. Defines the groups that are used for clustering (e.g 'group = c("test", "control")').
#' 		  It uses \code{grep} function to locate the group names from the CEL file names and then builds
#' 		  index vectors determining which files belong to which groups.
#' @param norm.method A \code{function}. Defines a function used to normalize the raw expression values. The default normalization
#' 		  function is \code{norm.rank}.
#' @param filetag A \code{character}. This is a custom string that can be used to identify the experiment. It incorporates it to
#' 		  the names of the output files.
#' @param cluster A \code{list}. Holds the probes of different clusters in a \code{list}.
#' @param clust.method A \code{function}. Defines the function used for clustering. The user can
#' 		  build a custom clustering function. The input of the custom function needs
#' 		  to be a \code{DEMIClust} object and the output is a \code{list} of probes, where
#' 		  each list corresponds to a specific cluster. The default function is \code{demi.wilcox.test}
#' 		  that implements the \code{wilcox.test} function. However we recommend to use the function
#' 		  \code{demi.wilcox.test.fast} that uses a custom \code{wilcox.test} and runs a lot faster.
#' @param cutoff.pvalue A \code{numeric}. Sets the cut-off p-value used for determining statistical
#' 		  significance of the probes when clustering the probes into clusters.
#' @param pathway A \code{logical}. If set to TRUE the functional annotation analysis is done on top of differential expression
#' 		  analysis.
#' @return A list containing the \code{DEMIExperiment} object where differential expression results have been added to and
#' 		   a \code{data.frame} consisting of the functional annotation analysis results.
#' @seealso \code{DEMIExperiment}, \code{DEMIClust}, \code{DEMIPathway}, \code{DEMIDiff}, \code{demi.wilcox.test.fast}, \code{wilcox.test}
#' @details
#' 
#' Instead of automatically clustered probes \code{DEMIClust} object can use user defined lists of
#' probes for later calculation of differential expression. This is done by setting the \code{cluster} parameter. It
#' overrides the default behaviour and no actual clustering occurs. Instead the list of probes defined
#' in the \code{cluster} parameter are considered as already clustered probes. The list needs to contain proper names for
#' probe vectors so that they would be recognizable later. Also instead of using the default clustering method the user can
#' write his/her own function for clustering probes based on the expression values.
#' 
#' Further specification of the parameters:
#' \itemize{
#' 	\item{maxtargets}{
#' 		When \code{analysis} is set to 'gene' then all probes that match to more genes then allowed by \code{maxtargets} parameter will
#' 		not be included in the analysis. For 'transcript' and 'exon' analysis the number is also calculated on a gene
#' 		level. For example if \code{maxtargets} is set to one and a probe matches to two transcripts but on the same gene,
#' 		then this probe will still be used in the analysis. However if the probe matches two transcripts on different
#' 		genes then this probe will not be included in the analysis. For 'genome' analysis the probe in most cases matches
#' 		to two genomic sections because adjacent sections overlap by 50%. However this is considered as one match and the
#' 		probe will still be used in the analysis.
#' 		}
#' 	\item{norm.method}{
#' 		Every user can apply their own normalization method by writing a custom normalization function. The function should
#' 		take in raw expression matrix and return the normalized expression matrix where probe ID's are kept as rownames and
#' 		column names are CEL file names. The normalized expression matrix will then be stored as part of the \code{DEMIExperiment}
#' 		object.
#' 		}
#' 	\item{sectionsize}{
#' 		The \code{sectionsize} parameter defines the length of the genomic target region. Currenlty \code{sectionsize} can be set
#' 		as: 100000, 500000 and 1000000. All adjacent sections, except the ones on chromosome ends, overlap with the
#' 		next adjacent section by 50%. It ensures the all probes matching to genome will be assigned to at least one
#' 		genomic section. This parameter is required when \code{analysis} is set to 'genome'.
#' 		}
#' 	\item{group}{
#' 		All the CEL files used in the analysis need to contain at least one of the names specified in the \code{group}
#' 		parameter because they determine what groups to compare against each other. It is also a good
#' 		practice to name the CEL files to include their common features. However if a situation arises where the
#' 		group/feature name occurs in all filenames then the user can set group names with specific filenames
#' 		by seperating names in one group with the "|" symbol. For example \code{group = c( "FILENAME1|FILENAME2|FILENAME3",
#' 		"FILENAME4|FILENAME5|FILENAME6" )}. These two groups are then used for clustering the probes expression values.
#' 		}
#' 	\item{norm.method}{
#' 		The \code{norm.method} defines a function to use for the normalization of raw expression matrix. The user can implement his/her
#' 		own function for the normalization procedure. The function should take in raw expression matrix and return the normalized
#' 		expression matrix where probe ID's are kept as rownames and column names are CEL file names.
#' 		}
#' 	\item{clust.method}{
#' 		The user can write his/her own function for clustering probes according to their expression values. The custom
#' 		function should take \code{DEMIClust} object as the only parameter and output a \code{list}. The output list should contain
#' 		the name of the clusters and the corresponding probe ID's. For example \code{return( list( cluster1 = c(1:10), cluster2
#' 		= c(11:20), cluster3 = c(21:30) )}.
#' 		}
#' 	\item{cluster}{
#' 		This parameter allows to calculate differential expression on user defined clusters of probe ID's. It needs to be a
#' 		\code{list} of probe ID's where the \code{list} names correspond to the cluster names. For example \code{list( cluster1 = c(1:10),
#' 		cluster2(1:10) )}. When using this approach you need to make sure that all the probe ID's given in the clusters are available in
#' 		the analysis. Otherwise an error message will be produced and you need to remove those probes that have no alignment
#' 		in the analysis. When setting this parameter the default behaviour will be overridden and no default clustering
#' 		will be applied.
#' 		}
#' }
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
#' # Do DEMI analysis with functional annotation analysis
#' demires <- demi(analysis = 'gene', celpath = destfolder, group = c( "BRAIN", "UHR" ),
#' 		experiment = 'myexperiment', organism = 'homo_sapiens',
#' 		clust.method = demi.wilcox.test.fast, pathway = TRUE)
#'
#' # Do DEMI analysis without functional annotation analysis
#' demires <- demi(analysis = 'gene', celpath = destfolder, group = c( "BRAIN", "UHR" ),
#' 		experiment = 'myexperiment', organism = 'homo_sapiens',
#' 		clust.method = demi.wilcox.test.fast, pathway = FALSE)
#' 
#' # Retrieve results from the created object
#' head( getResultTable( demires$experiment ) )
#' 
#' }
#' 
#' @export 
#' @docType methods
#' @rdname DEMIWrap-methods
#' @import methods
"demi" <-
function( analysis 				= "transcript",
		  celpath				= character(),
		  experiment			= character(),
		  organism				= character(),
		  maxtargets			= 0,
		  maxprobes				= character(),
		  pmsize				= 25,
		  sectionsize			= character(),
		  group 				= character(),
		  norm.method			= norm.rrank,
		  filetag 				= character(),
		  cluster				= list(),
		  clust.method			= function(){},
		  cutoff.pvalue			= 0.05,
		  pathway				= logical() )
{
	
	# DEMIExperiment checking
	if ( missing( analysis ) == TRUE ) {
		#stop( paste( "Error:", sQuote( "analysis" ), "is unspecified" ) );
		stop( DEMIMessages$parameterMissing( "analysis" ) );
	} else if ( length( analysis ) > 1 ) {
		#stop( paste( "Error:", sQuote( "analysis" ), "can only be of length 1" ) );
		stop( DEMIMessages$tooManyParameters( "analysis", 1 ) );
	} else if ( analysis == "genome" ) {
		if ( missing( sectionsize ) == TRUE ) {
			#stop( paste( "Error:", sQuote( "sectionsize" ), "has not been specified for 'genome' analysis" ) );
			stop( DEMIMessages$sectionsizeMissing() );
		}
	}
	if ( missing( organism ) == TRUE ) {
		#stop( paste( "Error:", sQuote( "organism" ), "is unspecified" ) );
		stop( DEMIMessages$parameterMissing( "organism" ) );
	} else if ( length( organism ) > 1 ) {
		#stop( paste( "Error:", sQuote( "organism" ), "can only be of length 1" ) );
		stop( DEMIMessages$tooManyParameters( "organism", 1 ) );
	}
	if ( missing( celpath ) == TRUE ) {
		#stop( paste( "Error:", sQuote( "celpath" ), "is unspecified" ) );
		stop( DEMIMessages$parameterMissing( "celpath" ) );
	}
	if ( missing( experiment ) == TRUE ) {
		#stop( paste( "Error:", sQuote( "experiment" ), "is unspecified" ) );
		stop( DEMIMessages$parameterMissing( "experiment" ) );
	} else if ( length( experiment ) > 1 ) {
		#stop( paste( "Error:", sQuote( "analysis" ), "can only be of length 1" ) );
		stop( DEMIMessages$tooManyParameters( "experiment", 1 ) );
	}
	if ( length( maxtargets ) > 0 && is.numeric( maxtargets ) == FALSE ) {
		#stop( paste( "Error:", sQuote( "maxtargets" ), "has to be 0 or a positive integer" ) );
		stop( DEMIMessages$positiveIntegerOrZero( "maxtargets" ) );
	} else if ( length( maxtargets ) > 1 ) {
		#stop( paste( "Error:", sQuote( "maxtargets" ), "can only be of length 1" ) );
		stop( DEMIMessages$tooManyParameters( "maxtargets", 1 ) );
	}
	if ( length( pmsize ) > 0 && is.numeric( pmsize ) == FALSE ) {
		#stop( paste( "Error:", sQuote( "pmsize" ), "has to be a positive integer up to 25" ) );
		stop( DEMIMessages$positiveIntegerUpTo( "pmsize", 25 ) );
	} else if ( length( pmsize ) > 1 ) {
		#stop( paste( "Error:", sQuote( "pmsize" ), "can only be of length 1" ) );
		stop( DEMIMessages$tooManyParameters( "pmsize", 1 ) );
	}
	if ( length( sectionsize ) > 0 && is.numeric( sectionsize ) == FALSE ) {
		#stop( paste( "Error:", sQuote( "sectionsize" ), "has to be a positive integer" ) );
		stop( DEMIMessages$positiveInteger( "sectionsize" ) );
	} else if ( length( sectionsize ) > 1 ) {
		#stop( paste( "Error:", sQuote( "sectionsize" ), "can only be of length 1" ) );
		stop( DEMIMessages$tooManyParameters( "sectionsize", 1 ) );
	}
	if ( missing( group ) == TRUE ) {
		stop( DEMIMessages$parameterOfClassMissing( "group", "character" ) );
	} else if ( length( group ) != 2 ) {
		stop( DEMIMessages$tooManyParameters( "group", 2 ) );
	}
	if ( is.function( norm.method ) == FALSE ) {
		#stop( paste( "Error:", sQuote( "norm.method" ), "has to be a function" ) );
		stop( DEMIMessages$hasToBeFunction( "norm.method" ) );
	} else if ( length( norm.method ) > 1 ) {
		#stop( paste( "Error:", sQuote( "norm.method" ), "can only be of length 1" ) );
		stop( DEMIMessages$tooManyParameters( "norm.method", 1 ) );
	}
	if ( missing( sectionsize ) == FALSE && analysis != "genome" ) {
		#stop( paste( "Error:", sQuote( "sectionsize" ), "can only be used in 'genome' analysis\n" ) );
		stop( DEMIMessages$sectionsizeUse() );
	}
	if ( missing( cutoff.pvalue ) == FALSE ) {
		if ( is.numeric( cutoff.pvalue ) == FALSE ) {
			#stop( paste( "Error:", sQuote( "cutoff.pvalue" ), "is not of class 'numeric'\n" ) );
			stop( DEMIMessages$parameterNotOfClass( "cutoff.pvalue", "numeric" ) );
		}
	}
	if ( pathway == TRUE ) {
		if ( analysis != "gene" && analysis != "transcript" ) {
			stop( DEMIMessages$DEMIPathway$cantrun );
		}
	}
	
	#	exprsData - will be generated later
	demiexp <- new( "DEMIExperiment",
			analysis			= analysis,
			celpath				= celpath,
			experiment 			= experiment,
			organism 			= organism,
			maxtargets			= as.numeric( maxtargets ),
			maxprobes 			= as.character( maxprobes ),
			pmsize				= as.numeric( pmsize ),
			sectionsize			= as.numeric( sectionsize ),
			norm.method			= as.function( norm.method ),
			filetag 			= filetag
	);
	
	
	# DEMIClust checking
	demiGroup <- NULL;
	demiClusters <- list();
	#	if the user overrides the clustering by adding custom probe vector to DEMIClust object
	if ( missing( cluster ) == FALSE ) {
		#	if cluster name has been set
		if ( missing( group ) == TRUE ) {
			#stop( paste( "Error:", sQuote( "group" ), "is unspecified. You need to specify the name of you probe cluster\n" ) );
			stop( DEMIMessages$parameterMissing( "group" ) );
		} else if ( missing( group ) == FALSE ) {
			#	if the number of clusters and the names are same
			if ( length( cluster ) != length( group ) ) {
				#stop( paste( "Error:", sQuote( "group" ), "and", sQuote( "cluster" ), "are of different length\n" ) );
				stop( DEMIMessages$paramAndParamNotEqualLength( "group", "cluster" ) );
			} else {
				#	check if all the in the cluster are available in the blatTable
				for ( i in 1:length( cluster ) ) {
					probesOK = check4probe( experiment, cluster[[i]] );
					if ( is.null( probesOK ) == FALSE ) {
						stop( probesOK );
					}
				}
				#	since all the probes were present in the blat table continue
				for ( i in 1:length( cluster ) ) {
					demiClusters[i] <- cluster[i];
				}
				names( demiClusters ) <- group;
				demiGroup <- DEMIGroup( groupNames = group );
			}
		}
	} else if( missing( cluster ) == TRUE ) { #	if the user uses default parameters
		#	check if clust.method is a method
		if ( is.function( clust.method ) == FALSE ) {
			#stop( paste( "Error:", sQuote( "clust.method" ), "has to be a function\n" ) );
			stop( DEMIMessages$hasToBeFunction( "clust.method" ) );
		}
		#	set wilcox.test as the default clustering method
		if ( missing( clust.method ) == TRUE ) {
			clust.method = demi.wilcox.test;
		}
		#	check for correct definitions
		if ( is.vector( group ) == FALSE || length( group ) != 2 ) {
			#stop( paste( "Error:\tThere can only be two character elements in the", sQuote( "group" ), "vector\n",
			#				"\t\tAll group names need to be present in CEL-file names in your CEL-file directory specified by", sQuote( "celpath"), "parameter\n",
			#				"\t\tYou can express one group with regular expression by seperating groups in a string with a '|' symbol, for we use", sQuote( "grep" ), "to find group indexes from CEL file names\n",
			#				"\t\ti.e. '> groups = c( \"TUMOR\", \"NORMAL\" )'\n")
			#);
			stop( DEMIMessages$wrongDefinition() );
		} else {
			demiGroup <- DEMIGroup( groupA = group[1], groupB = group[2] );
		}
	}
	
	demiclust <- new( "DEMIClust",
			experiment 		= demiexp,
			group			= demiGroup,
			clust.method	= clust.method,
			cluster			= demiClusters,
			cutoff.pvalue	= cutoff.pvalue
	);

	
	# DEMIDiff checking
	
	diffName <- paste( demiclust@group@groupA, "_VS_", demiclust@group@groupB, sep = "" );
	
	demidiff <- new( "DEMIDiff",
			cluster				= demiclust,
			name				= diffName
	);
	
	finalexp <- attachResult( demiexp, demidiff );
	
	if ( missing( filetag ) == TRUE ) {
		filename <- paste( "demi_", getExperiment( finalexp ), "_", getAnalysis( finalexp ), ".txt", sep = "" );
	} else {
		filename <- paste( "demi_", getExperiment( finalexp ), "_", getAnalysis( finalexp ), "_", finalexp@filetag, ".txt", sep = "" );
	}
	#cat( paste( "*Saving analysis results into file ", filename, " in your working directory\n", sep = "" ) );
	cat( DEMIMessages$savingWhatTo( "analysis results", paste( filename, "in your working directory" ) ) );
	restab <- getResultTable( finalexp );
	write.table( restab, file = filename, sep = "\t", quote = FALSE, row.names = FALSE )
	
	if ( missing( pathway ) == TRUE ) {
		return( list( experiment = finalexp, pathway = NULL ) );
	} else if ( missing( pathway ) == FALSE & pathway == TRUE ) {
		demipath <- DEMIPathway( demidiff );
		if ( missing( filetag ) == TRUE ) {
			#cat( paste( "*Saving pathway results into file ", paste( "demi_", getExperiment( finalexp ), "_pathway.txt", sep = "" ) ," in your working directory\n", sep = "" ) );
			cat( DEMIMessages$savingWhatTo( "pathway results", paste( paste( "demi_", getExperiment( finalexp ), "_pathway.txt", sep = "" ), "in your working directory" ) ) );
			write.table( demipath, file = paste( "demi_", getExperiment( finalexp ), "_pathway.txt", sep = "" ), sep = "\t", quote = FALSE, row.names = FALSE );
		} else {
			#cat( paste( "*Saving pathway results into file ", paste( "demi_", getExperiment( finalexp ), "_pathway_", finalexp@filetag, ".txt", sep = "" ) ," in your working directory\n", sep = "" ) );
			cat( DEMIMessages$savingWhatTo( "pathway results", paste( paste( "demi_", getExperiment( finalexp ), "_pathway_", finalexp@filetag, ".txt", sep = "" ), "in your working directory" ) ) );
			write.table( demipath, file = paste( "demi_", getExperiment( finalexp ), "_pathway_", finalexp@filetag, ".txt", sep = "" ), sep = "\t", quote = FALSE, row.names = FALSE );
		}
		return( list( experiment = finalexp, pathway = demipath ) );
	}
	
	return( list( experiment = finalexp, pathway = NULL ) );
	
}#demi
