#==============================================================================#
# messages.R:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# DEMIMessages - contains all the messages as a list
#==============================================================================#

#' A \code{list} of DEMI messages
#' 
#' This \code{list} contains all the messages shown to the user when running DEMI analysis. It
#' is used internally in DEMI analysis.
#' 
#' @author Sten Ilmjarv, Hendrik Luuk
#' @name DEMIMessages
#' @docType data
#' @keywords data
DEMIMessages <- list(
		#---------------------------------------
		# cluster-methods.R
		#---------------------------------------
		demiequal = list(
				main = "# Clustering probes into 'equal' cluster based on Wilcoxon rank sum test if the groups are bigger than 2\n",
				custom.approach = "\tUsing an optimized implementation to perform Wilcoxon rank sum test\n",
				standard.approach = "\tUsing the native implementation of the Wilcoxon rank sum test (function 'wilcox.test')\n"
				),
		demi.wilcox.test.fast = list(
				main = "# Clustering probes into 'higher' and 'lower' clusters based on Wilcoxon rank sum test\n",
				custom.approach = "\tUsing an optimized implementation to perform Wilcoxon rank sum test\n",
				standard.approach = "\tUsing the native implementation of the Wilcoxon rank sum test (function 'wilcox.test')\n"
				),
		demi.wilcox.test = list(
				main = "# Clustering probes into 'higher' and 'lower' clusters based on Wilcoxon rank sum test\n",
				exact.false = "\tUsing wilcox with parameter 'exact = FALSE' since sample sizes are sufficiently large.\n",
				exact.true = "\tUsing wilcox with parameter 'exact = TRUE'.\n"
				),
		demiExperiment_t.test = list(
				main = paste( "# Clustering probes into 'higher' and 'lower' clusters based on t-test (function 't.test')\n" )
		),
		demi.comp.test = list(
				main = "# Clustering probes into 'higher' and 'lower' clusters\n"
				),
		takestime = function(no.of.probes) { paste( "\tIt can take some time for there are", no.of.probes, "probes to be analyzed.\n" ) },
		probesdone = function(no.of.probes) { paste( "\t\t", no.of.probes, " probes done\n", sep = "" ) },
		#---------------------------------------
		# Constructor.R
		#---------------------------------------
#		sectionsizeMissing = paste( sQuote( "sectionsize" ), "has not been specified for 'genome' analysis" ),
#		sectionsizeMissing = sprintf( "%s has not been specified for 'genome' analysis", sQuote( "sectionsize" ) ),
#		sectionsizeMissing = paste( function(x) {before <- after <- "'"; paste0(before, x, after)}, "has not been specified for 'genome' analysis" ),
		sectionsizeMissing = function() { paste( sQuote( "sectionsize" ), "has not been specified for 'genome' analysis" ) },
		positiveIntegerOrZero = function( parameter ) { paste( sQuote( parameter ), "has to be 0 or a positive integer" ) },
#		sectionsizeUse = paste( sQuote( "sectionsize" ), "can only be used in 'genome' analysis" ),
		sectionsizeUse = function() { paste( sQuote( "sectionsize" ), "can only be used in 'genome' analysis" ) },
		parameterOfClassMissing = function( parameter, cls) { paste( sQuote( parameter ), "of class", sQuote( cls ),"is unspecified" ) },
		paramAndParamNotEqualLength = function( param1, param2 ) { paste( sQuote( param1 ), "and", sQuote( param2 ), "need to be of same length" ) },
		wrongDefinition = function() { paste( paste( "\t\tThere can only be two character elements in the ", sQuote( "group" ), " vector\n", sep = "" ),
								 paste( "\t\tAll group names need to be present in CEL-file names in your CEL-file directory specified by ", sQuote( "celpath"), " parameter\n", sep = "" ),
								 paste( "\t\tYou can express one group with regular expression by separating groups in a string with a '|' symbol, for we use ", sQuote( "grep" ), " to find group indexes from CEL file names\n", sep = "" ),
								 paste( "\t\ti.e. '> groups = c( \"TUMOR\", \"NORMAL\" )'\n", sep = "" ), sep = "" ) },
		savingWhatTo = function( what, to ) { paste( "# Saving", what, "to", to,"\n" ) },
		#---------------------------------------
		# custom-methods.R
		#---------------------------------------
		customTargets = list(
				whatMissingFrom = function( what, from ) { paste( "column", sQuote( what ), "is missing from the", sQuote( from ), "object." ) },
				pmsizeMissing = paste( "# Note: column", sQuote( "pmsize" ), "is unspecified in the", sQuote( "blat" ), "object. We will therefore assume that all probe matches to the target are perfect matches with 100% identity.\n" ),
				startMissing = paste( "# Note: column", sQuote( "start" ), "is unspecified in the", sQuote( "blat" ), "object.\n" ),
				strandMissing = paste( "# Note: column", sQuote( "strand" ), "is unspecified in the", sQuote( "blat" ), "object. We will therefore assume that all probe matches were on the major matching strand.\n" ),
				probesNotPresent = function( notpresent ) { paste( "The following probes were not present on the microarray. Remove them from your data and try again.\n", paste( notpresent, collapse = "," ) ) },
				geneMissing = paste( "# Note: column", sQuote( "geneID" ), "is missing from the", sQuote( "anno" ), "object. Will duplicate", sQuote( "targetID" ), "as", sQuote( "geneID" ), "\n" ),
				annoMissing = paste( "# Note:", sQuote("anno"), "has not been specified. Creating custom annotation\n" ),
				ignoreStrand = function( strand ) { paste( "# Note: will ignore all", strand, "and NA strand matches\n" ) },
				pmsizeSmallerRemove = "# Note: Removing probes with perfect match size less than indicated in the DEMIExperiment object\n",
				determineMatches = "# Determining the number of hits on each target for every probe\n",
				error = "You can't add custom region alignments's and annotation info to the exon analysis. You can only add them to 'transcript' and 'gene' analyses.\n"
				),
		DEMIPathway = list(
				main = "# Performing pathway analysis of differentially expressed genes\n",
				cantrun = "DEMIPathway can only be run if the analysis type in DEMIExperimebt object is a \"gene\" or \"transcript\"\n"
				),
		#---------------------------------------
		# DEMIClust-methods.R
		#---------------------------------------
		DEMIClust = list(
				noCELFilesWithGroupname = function( groupName ) {paste( "\tThere are no CEL files that include the group name", paste( "'", groupName, "'", sep = "" ), "\n",
								"\tThe", sQuote( "groups" ), "is a vector of group-tags that are included in CEL file names to denote the group-membership\n") },
				usingBuiltIn = "# Note: Clustering with demi built-in functions\n",
				usingCreaterOrLower = "# Note: A heuristic-based method will be used to assess probe-level differential expression as one of the group sizes is <= 2.\n",
				usingUserProvided = "# Note: Clustering with user-provided function\n"
				),
		#---------------------------------------
		# DEMIDiff-methods.R
		#---------------------------------------
		DEMIDiff = list(
				zeroTargetsFound = "0 specified targets were found in the experiment",
				calcDiffExp = "# Calculating differential expression ",
				betweenGroups = function( groupA, groupB ) { paste( "for comparison between groups", paste( "'", groupA, "'", sep = "" ), "and", paste( "'", groupB, "'", sep = "" ), "\n" ) },
				onUserDefined = "on user-defined clusters\n",
				analyzingClusterNative = function( clusterName, groupA, groupB ) { paste( "\tAnalysing cluster", paste( "'", clusterName, "'", sep = "" ), "in", paste( "'", groupA, "'", sep = "" ), "compared to", paste( "'", groupB, "'", sep = "" ), "\n" ) },
				analyzingClusterCustom = function( clusterName ) { paste( "\tAnalysing cluster", paste( "'", clusterName, "'", sep = "" ), "\n" ) },
				calcHyperGeo = "\t\tCalculating hypergeometric distribution p-values ... ",
				annotatingResults = "\t\tAnnotating results ...",
				multipleCorrectionBH = "\t\tAdjusting p-values for multiple testing by the FDR method (Benjamini & Hochberg, 1995) ... ",
				multipleCorrectionBY = "\t\tAdjusting p-values for multiple testing by the FDR method under dependence (Benjamini & Yekutieli, 2001) ... ",
				multipleCorrectionBONF = "\t\tAdjusting p-values for multiple testing by the Bonferroni method ... ",
				diffExpDone = "# Calculating differential expression ... DONE\n"
				),
		#---------------------------------------
		# DEMIExperiment-methods.R
		#---------------------------------------
		DEMIExperiment = list (
				#	"Error: " will be added to the front of the variables that start with invalid
				invalidWhatChooseFrom = function( what, chooseFrom ) { paste( sQuote( what ), "is invalid. You need to choose", sQuote( what ), "from (", paste( chooseFrom, collapse = ", " ),")\n" ) },
				notFolder = function( parameter ) { paste( sQuote( parameter ), "is not a folder\n" ) },
				folderEmpty = function( parameter ) { paste( sQuote( parameter ), "does not contain any CEL files - no CEL suffix found\n" ) },
				notACELFile = function( file, parameter ) { paste( "The ", file, " in ", sQuote( parameter ), " is not a CEL file!\n", sep = "" ) },
				maxprobesError = "The parameter 'maxprobes' has to be a positive integer or set as 'median', 'max' or not defined at all which by default set's it to 'max'",
				probesExcluded = function( excludedProbes ) { paste( "Warning: these probes were excluded because they were not present in the analysis:\n", paste( excludedProbes, collapse = ", " ), "\n" ) },
				noTargetsFound = "No targets were found",
				noProbesMeasuring = "No probes measuring the specified targets matched the experiment criteria",
				attachErrorResultsExists = function( variable ) { paste( "Warning: Can't add the variable", sQuote( deparse( substitute( variable ) ) ), "of class 'DEMIDiff' to the object of class 'DEMIExperiment' because the results of that object already exist in the 'DEMIExperiment' object\n" ) },
				attachErrorDiffDEMIExp = "Can't add results of class 'DEMIDiff' to this experiment of class 'DEMIExperiment' because they represent different 'DEMIExperiment' objects",
				loadingCEL = "# Loading CEL files and creating expression matrix\n",
				includedCELFiles = function( celfiles ) { paste( "\tThe CEL files included are: ", celfiles, "\n", sep = "" ) },
				usingMicroarrayPlatform = function( platform ) { paste( "\tUsing microarray platform ", platform, "\n", sep = "" ) },
				libraryNotInstalled = function( packageName ) { paste( "Library - package", packageName, "not installed\n" ) },
				searchingInternetRepo = function( packageName ) { paste( "Searching the internet repository for package ", packageName, "\n", sep = "" ) },
				packageNotFoundFromRepo = function( packageName, package_url ) { paste( "Package", packageName, "was not found from the", url, " url\n" ) },
				installingPackage = function( packageName, repository ) { paste( "Installing package", packageName, "from ", repository, "repository\n") },
				checkYourInternet = "The current operation could not access the internet. Please check your internet connectction",
				loadingAnnoSuccess = function( packageName ) { paste( "# Loading required annotation data from ", packageName, "\n", sep = "" ) },
				loadingAnnoFail = function( annoTableName, packageName ) { paste( "Can't load required annotation table ", annoTableName, " from the package ", packageName, "\n", sep = "" ) },
				sectionsizeNotAvail = function( sectionsizes ) { paste( "The specified", sQuote( "sectionsize" ), "is not available. Try", paste( sectionsizes, collapse = ", " ) ) },
				cantLoadWhatTable = function( whatTable ) { paste( "Can't", whatTable, "table with the specified parameters" ) },
				loadingBlatSuccess = function( packageName ) { paste( "# Loading required alignment data from ", packageName, "\n", sep = "" ) },
				blatDoesNotExists = function( blat, platform, packageName ) { paste( "The required probe alignment table", blat, "for the array", platform,"does not exists in the package", packageName ) },
				pmsizeNotAvail = function( pmsizes ) { paste( "The specified ", sQuote( "pmsize" ), " is not available. Available pmsize's are ", paste( pmsizes, collapse = ", " ) , sep = "" ) },
				ignoreStrandMatches = function( strand ) { paste( "\tWill ignore the '", strand, "' matches\n", sep = "" ) },
				loadingCytoband = "# Loading cytoband information\n",
				loadingPathway = "# Loading pathway information\n",
				check4probeError = function( difference ) { paste( "\tError:\tThere is", difference, "probe(s) that was/were not located in the alignment table.\n",
							"\t\tOnly probes in the alignment table can be used to populate custom probe vector cluster\n",
							"\t\tUse the function 'getAlignment' on your 'DEMIExperiment' object to view probes present in alignment table\n" ) },
				check4targetError = function( notfound ) { paste( "The target(s)", paste( notfound, collapse = ", " ), "either had no probes or was/were not found in annotation information" ) },
				zeroTargetsFound  = "0 specified targets were found in the experiment"
				),
		#---------------------------------------
		# DEMIResults-methods.R
		#---------------------------------------
		DEMIResults = list(
				resultsContain = "the results list can only contain objects of class 'DEMIResult'",
				resultsUndefined = "no results have been defined"
				),
		#---------------------------------------
		# diffexp-methods.R
		#---------------------------------------
		diffexp = list(
				matchesCluster = "\t\tCalculating matches for probes in cluster ... ",
				matchesAll = "\t\tCalculating matches over all probes ... "
				),
		#---------------------------------------
		# normalization-methods.R
		#---------------------------------------
		normalization = list(
				main = "# Normalizing expression values ",
				normrrank = " - using 'relative rank' as the normalization method\n",
				normquantile = "- using 'quantile normalization' as the normalization method\n"
		),
		#---------------------------------------
		# global
		#---------------------------------------
		parameterMissing = function( parameter ) { paste( sQuote( parameter ), "is unspecified" ) },
		parameterNotOfClass = function( parameter, cls) { paste( sQuote( parameter ), "is not of class", sQuote( cls ) ) },
		isNotError = function( parameter, what ) { paste( sQuote( parameter ), "is not a", what ) },
		hasToBeNumericBetween = function( parameter, start, end ) { paste( sQuote( parameter ), "has to be set between", start, "and", end, "\n" ) },
		tooManyParameters = function( parameter, len ) { paste( sQuote( parameter ), "can only be of length", len ) },
		positiveInteger = function( parameter ) { paste( sQuote( parameter ), "has to be a positive integer" ) },
		positiveIntegerUpTo = function( parameter, to ) { paste( sQuote( parameter ), "has to be a positive integer up to", to ) },
		hasToBeFunction = function( parameter ) { paste( sQuote( parameter ), "has to be a function" ) },
		notEmptyString = function( parameter ) { paste( sQuote( parameter ), "can't be an empty string" ) },
		done = function() { "DONE\n" },
		#---------------------------------------
		# random
		#---------------------------------------
		demiStart = function() { paste("\n",
				"\t#-----#        #---------#     #\        #    #",
				"\t|       \\      |               |\\      /|    |",
				"\t|        #     |               | \\    / |    |",
				"\t|        |     |               |  \\  /  |    |",
				"\t|        |     #------         |   \\/   |    |",
				"\t|        |     |               |        |    |",
				"\t|        #     |               |        |    |",
				"\t|       /      |               |        |    |",
				"\t#-----#        #---------#     #        #    |",
				"\tDifferential Expression from Multiple Indicators","",
				"Authors: Sten Ilmjarv <sten.ilmjarv@gmail.com>, Hendrik Luuk <hendrik.luuk@gmail.com>\n\n", sep = "\n" ) }
)
