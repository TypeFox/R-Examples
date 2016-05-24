#==============================================================================#
# AllClasses.R
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# DEMICel
# DEMIGroup
# DEMIExperiment
# DEMIClust
# DEMIResult
# DEMIDiff
#==============================================================================#

#------------------------------------------------------------------------------#
# DEMICel:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' Class \code{DEMICel}
#' 
#' The class \code{DEMICel} holds the raw and normalized expression matrices.
#' 
#' @param celMatrix A \code{matrix}. The raw expression matrix.
#' @param normMatrix A \code{matrix}. The normalized expression matrix.
#' 
#' @author Sten Ilmjarv
#' 
#' @name DEMICel-class
#' @rdname DEMICel-class
#' @exportClass DEMICel
#' @import methods
setClass( "DEMICel",
		representation( celMatrix	= "matrix",
						normMatrix	= "matrix"
		),
		prototype( celMatrix	= matrix( nr = 0 , nc = 0 ),
				   normMatrix	= matrix( nr = 0 , nc = 0 )
		)
)


#------------------------------------------------------------------------------#
# DEMIGroup:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' Class \code{DEMIGroup}
#' 
#' The class \code{DEMIGroup} holds the information about the groups. Such as
#' the column indexes of both groups and names of the groups.
#' 
#' @param groupA A \code{character}. Holds the name of group A.
#' @param groupB A \code{character}. Holds the name of group B.
#' @param indexA A \code{numeric}. A vector of column indexes belonging to group A.
#' @param indexB A \code{numeric}. A vector of column indexes belonging to group B.
#' @param groupNames A \code{character}. Holds the names of custom groups created by
#' 		  the user.
#' 
#' @author Sten Ilmjarv
#' 
#' @name DEMIGroup-class
#' @rdname DEMIGroup-class
#' @exportClass DEMIGroup
#' @import methods
setClass( "DEMIGroup",
		representation( groupA	= "character",
						groupB	= "character",
						indexA	= "numeric",
						indexB	= "numeric",
						groupNames = "character"
		),
		prototype( groupA = character(),
				   groupB = character(),
				   indexA = numeric(),
				   indexB = numeric(),
				   groupNames = character()
		)
)#DEMIGroup


#------------------------------------------------------------------------------#
# DEMIExperiment:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' Class \code{DEMIExperiment}
#' 
#' The class \code{DEMIExperiment} defines an experiment. It holds the raw and normalized expression
#' data as well as annotation information for selected analysis (either 'gene', 'transcript',
#' 'exon' or 'genome'). It can be used to hold all the analysis results (\code{DEMIDiff} objects) done
#' on the same \code{DEMIExperiment} object.
#' 
#' @param analysis A \code{character}. Defines the analysis type. It can be either 'transcript',
#' 		  'gene', 'exon' or 'genome'. The default value is 'transcript'. For 'genome' analysis
#' 		  \code{sectionsize} parameter needs to be defined as well.
#' @param celpath A \code{character}. It can point to the directory containing the CEL files or is a vector
#' 		  that points directly to the CEL files.
#' @param experiment A \code{character}. A custom name of the experiment defined by the user (e.g. 'myexperiment').
#' @param organism A character. The name of the species the micrroarrays are measuring (e.g. 'homo_sapiens'
#' 		  or 'mus_musculus') given in lowercase and words separated by underscore.
#' @param arraytype A \code{character}. Holds the platform name of the microarrays used in the analysis.
#' @param maxtargets A \code{numeric}. The maximum number of allowed targets (e.g. genes or transcripts) one probe
#' 		  can have a match against. If to set it to 1 it means that the probe can match only one gene. If the \code{analysis}
#' 		  is set to 'transcript' the program still calculates the number of matches on genes. Hence a probe
#' 		  matching two transcripts on the same gene would be included but a probe matching two transcripts
#' 		  on different genes would not be included. The value needs to be a positive integer or 0.
#' @param maxprobes A \code{character}. Sets the number of unique probes a target is allowed to have a match against.
#' 		  All the targets that yield more alignments to different probes then set by \code{maxprobes} will be scaled
#' 		  down to the number defined by the \code{maxprobes} parameter. It can be either a positive integer or set as
#' 		  'median' or 'max' - 'median' meaning the median number of probes matching to all targets and 'max'
#' 		  meaning the maximum number of probes matching to a target.
#' @param pmsize A \code{numeric}. The minimum number of consecutive nucleotides that need to match perfectly
#' 		  against the target sequence. It can be either 23, 24 or 25. This means that alignments with
#' 		  smaller perfect match size will not be included in the experiment. 
#' @param sectionsize A \code{numeric}. This is only used if the \code{analysis} parameter is set to 'genome'. It defines the
#' 		  length of the genomic target region used in the 'genome' analysis.
#' @param norm.method A \code{function}. Defines a function used to normalize the raw expression values. The function should take
#' 		  in raw expression matrix and return the normalized expression matrix where probe ID's are kept as rownames and
#' 		  column names are CEL file names.
#' @param filetag A \code{character}. This is a custom string that can be used to identify the experiment. At the current
#' 		  development stage this parameter is used only when using the function \code{demi}, where the output files will
#' 		  contain the specified filetag.
#' @param annoTable A \code{data.frame}. Holds the annotation information used in the experiment.
#' @param blatTable A \code{data.frame}. Holds the alignment information of probes and their corresponding targets.
#' @param cytoband A \code{data.frame}. Only used in the 'genome' analysis. Holds the karyotype information for every
#' 		  chromosome of the species specified by the \code{organism} parameter.
#' @param pathway A \code{data.frame}. Only used in the 'gene' and 'transcript' analysis. Holds the genes for every gene
#' 		  ontology category.
#' @param exprsdata A \code{DEMICel} object. Holds the raw and normalized expression matrices in a \code{DEMICel} object.
#' @param results A \code{list}. Can be used to store all the results as \code{DEMIDiff} objects done on the same \code{DEMIExperiment}
#' 		  object.
#' 
#' @author Sten Ilmjarv
#' 
#' @name DEMIExperiment-class
#' @rdname DEMIExperiment-class
#' @exportClass DEMIExperiment
#' @import methods
setClass( "DEMIExperiment",
		representation( analysis			= "character",
						celpath				= "character",
						experiment			= "character",
						organism			= "character",
						arraytype			= "character",
						maxtargets			= "numeric",
						maxprobes			= "character",
						pmsize				= "numeric",
						sectionsize			= "numeric",
						norm.method			= "function",
						filetag 			= "character",
						annoTable			= "data.frame",
						blatTable			= "data.frame",
						cytoband			= "data.frame",
						pathway				= "data.frame",
						exprsData			= "DEMICel",
						results				= "list"
		),
		prototype( analysis				= "transcript",
				   celpath 				= character(),
				   experiment 			= character(),
				   organism 			= character(),
				   arraytype 			= character(),
				   maxtargets 			= as.numeric( 0 ),
				   maxprobes 			= character(),
				   pmsize 				= as.numeric( 25 ),
				   sectionsize			= numeric(),
				   norm.method			= function(){},
				   filetag 				= character(),
				   annoTable			= data.frame( matrix( nr = 0 , nc = 0 ) ),
				   blatTable			= data.frame( matrix( nr = 0 , nc = 0 ) ),
				   cytoband				= data.frame( matrix( nr = 0 , nc = 0 ) ),
				   pathway				= data.frame( matrix( nr = 0 , nc = 0 ) ),
				   exprsData			= new( "DEMICel" ),
				   results				= list()
		)
)#DEMIExperiment


#------------------------------------------------------------------------------#
# DEMIClust:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' Class \code{DEMIClust}
#' 
#' The class \code{DEMIClust} stores the probe clusters in a \code{DEMIClust} object.
#'
#' @param experiment A \code{DEMIExperiment} object. Holds the \code{DEMIExperiment} object whose metadata
#' 		  (such as normalized expression values) is used to cluster the probes.
#' @param group A \code{DEMIGroup} object. Defines the groups that are used for clustering. The \code{DEMIGroup}
#' 		  object uses CEL file names to determine which files belong to which group. It uses
#' 		  \code{grep} function to locate the group names from the CEL file names.
#' @param clust.method A \code{function}. Defines the function used for clustering. The user can
#' 		  build a custom clustering function. The input of the custom function needs
#' 		  to be the same \code{DEMIClust} object and the output is a list of probes, where
#' 		  each list corresponds to a specific cluster. The default function is \code{demi.wilcox.test} that
#' 		  utilizes \code{wilcox.test}.
#' @param cluster A \code{list}. Holds the probes of different clusters in a \code{list}.
#' @param cutoff.pvalue A \code{numeric}. Sets the cut-off p-value used for determining statistical
#' 		  significance of the probes when clustering the probes into clusters. Default is 0.05.
#' 
#' @author Sten Ilmjarv
#' 
#' @name DEMIClust-class
#' @rdname DEMIClust-class
#' @exportClass DEMIClust
#' @import methods
setClass( "DEMIClust",
		representation( experiment		= "DEMIExperiment",
						group			= "DEMIGroup",
						clust.method	= "function",
						cluster			= "list",
						cutoff.pvalue	= "numeric"
		),
		prototype( experiment		= NULL,
				   group			= NULL,
				   clust.method		= function(){},
				   cluster			= list(),
				   cutoff.pvalue	= 0.05
		)
)#DEMIClust


#------------------------------------------------------------------------------#
# DEMIResult:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' Class \code{DEMIResult}
#' 
#' The class \code{DEMIResult} holds the results of the specified clusters.
#' 
#' @param group A \code{DEMIGroup} object. Holds the information about the groups that were used
#' 		  for clustering.
#' @param result A \code{list}. Holds the analysis results for each cluster.
#' 
#' @author Sten Ilmjarv
#' 
#' @name DEMIResult-class
#' @rdname DEMIResult-class
#' @exportClass DEMIResult
#' @import methods
setClass( "DEMIResult",
		representation( group	= "DEMIGroup",
						result	= "list"
		),
		prototype( group	= NULL,
				   result	= list()
		)
)


#------------------------------------------------------------------------------#
# DEMIDiff:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' Class \code{DEMIDiff}
#' 
#' The class \code{DEMIDiff} holds the analysis results, results of clustering and the
#' original metadata of the experiment.
#' 
#' @param cluster A \code{DEMIClust} object. Holds information about the clusters in a \code{DEMIClust}
#' 		  object that it itself contains the \code{DEMIExperiment} object were experiment
#' 		  metadata such as alignment information and annotation is held.
#' @param name A \code{character}. A specific name of the differential expression
#' 		  analysis (e.g. 'BRAINvsUHR') that is generated according to the group names.
#' @param result A \code{DEMIResult} object. Holds the \code{DEMIResult} object of the analysis.
#' 
#' @author Sten Ilmjarv
#' 
#' @name DEMIDiff-class
#' @rdname DEMIDiff-class
#' @exportClass DEMIDiff
#' @import methods
setClass( "DEMIDiff",
		representation( cluster				= "DEMIClust",
						name				= "character",
						result				= "DEMIResult"
						
		),
		prototype( cluster				= NULL,
				   name					= character(),
				   result				= new( "DEMIResult" )
				   
		)
)#DEMIDiff

