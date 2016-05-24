# generic methods



utils::globalVariables(c("exome","transcriptome","genome","sectionsizes","pmsizes","maxtarget","cytoband","pathway","Files","value","probeID"))

#
#	DEMIGroup
#

#------------------------------------------------------------------------------#
# DEMIGroup get functions
#------------------------------------------------------------------------------#

#' Returns the \code{indexA} parameter
#' 
#' Returns the \code{indexA} parameter of the \code{DEMIGroup} object. It is a \code{numeric} that
#' represents the column indexes of group A.
#' 
#' @param object A \code{DEMIGroup} object.
#' @return Returns the \code{indexA} parameter of the \code{DEMIGroup} object that is a \code{numeric}.
#' @seealso \code{DEMIGroup}
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
#' # Now we can continue the example of the function getIndexA
#' 
#' # Set up an experiment
#' demiexp <- DEMIExperiment( analysis = 'gene', celpath = destfolder,
#' 			experiment = 'myexperiment', organism = 'homo_sapiens' )
#' 
#' # Create clusters with an optimized wilcoxon's rank sum test incorporated within demi that
#' # precalculates the probabilities
#' demiclust <- DEMIClust( demiexp, group = c( "BRAIN", "UHR" ), clust.method = demi.wilcox.test.fast )
#' 
#' # Calcuate differential expression
#' demidiff <- DEMIDiff( demiclust )
#' 
#' # Get group A indexes
#' getIndexA( getGroup( demiclust ) )
#' getIndexA( getGroup( demidiff ) )
#' 
#' }
#' 
#' @export
#' @docType methods
#' @rdname getIndexA-methods
#' @import methods
setGeneric("getIndexA",					function(object) standardGeneric("getIndexA"));

#' Returns the \code{indexB} parameter
#' 
#' Returns the \code{indexB} parameter of the \code{DEMIGroup} object. It is a \code{numeric} that
#' represents the column indexes of group B.
#' 
#' @param object A \code{DEMIGroup} object.
#' @return Returns the \code{indexB} parameter of the \code{DEMIGroup} object that is a \code{numeric}.
#' @seealso \code{DEMIGroup}
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
#' # Now we can continue the example of the function getIndexB
#' 
#' # Set up an experiment
#' demiexp <- DEMIExperiment( analysis = 'gene', celpath = destfolder,
#' 			experiment = 'myexperiment', organism = 'homo_sapiens' )
#' 
#' # Create clusters with an optimized wilcoxon's rank sum test incorporated within demi that
#' # precalculates the probabilities
#' demiclust <- DEMIClust( demiexp, group = c( "BRAIN", "UHR" ), clust.method = demi.wilcox.test.fast )
#' 
#' # Calcuate differential expression
#' demidiff <- DEMIDiff( demiclust )
#' 
#' # Get group B indexes
#' getIndexB( getGroup( demiclust ) )
#' getIndexB( getGroup( demidiff ) )
#' 
#' }
#' 
#' @export
#' @docType methods
#' @rdname getIndexB-methods
#' @import methods
setGeneric("getIndexB",					function(object) standardGeneric("getIndexB"));

#' Returns the \code{groupA} parameter
#' 
#' Returns the \code{groupA} parameter of the \code{DEMIGroup} object. It is a \code{character} that
#' represents the name of group A.
#'  
#' @param object A \code{DEMIGroup} object.
#' @return Returns the \code{groupA} parameter of the \code{DEMIGroup} object that is a \code{character}.
#' @seealso \code{DEMIGroup}
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
#' # Now we can continue the example of the function getGroupA
#' 
#' # Set up an experiment
#' demiexp <- DEMIExperiment( analysis = 'gene', celpath = destfolder,
#' 			experiment = 'myexperiment', organism = 'homo_sapiens' )
#' 
#' # Create clusters with an optimized wilcoxon's rank sum test incorporated within demi that
#' # precalculates the probabilities
#' demiclust <- DEMIClust( demiexp, group = c( "BRAIN", "UHR" ), clust.method = demi.wilcox.test.fast )
#' 
#' # Calcuate differential expression
#' demidiff <- DEMIDiff( demiclust )
#' 
#' # Get group A name
#' getGroupA( getGroup( demiclust ) )
#' getGroupA( getGroup( demidiff ) )
#' 
#' }
#' 
#' @export
#' @docType methods
#' @rdname getGroupA-methods
#' @import methods
setGeneric("getGroupA",					function(object) standardGeneric("getGroupA"));

#' Returns the \code{groupB} parameter
#' 
#' Returns the \code{groupB} parameter of the \code{DEMIGroup} object. It is a \code{character} that
#' represents the name of group B.
#'  
#' @param object A \code{DEMIGroup} object.
#' @return Returns the \code{groupB} parameter of the \code{DEMIGroup} object that is a \code{character}.
#' @seealso \code{DEMIGroup}
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
#' # Now we can continue the example of the function getGroupB
#' 
#' # Set up an experiment
#' demiexp <- DEMIExperiment( analysis = 'gene', celpath = destfolder,
#' 			experiment = 'myexperiment', organism = 'homo_sapiens' )
#' 
#' # Create clusters with an optimized wilcoxon's rank sum test incorporated within demi that
#' # precalculates the probabilities
#' demiclust <- DEMIClust( demiexp, group = c( "BRAIN", "UHR" ), clust.method = demi.wilcox.test.fast )
#' 
#' # Calcuate differential expression
#' demidiff <- DEMIDiff( demiclust )
#' 
#' # Get group B name
#' getGroupB( getGroup( demiclust ) )
#' getGroupB( getGroup( demidiff ) )
#' 
#' }
#' 
#' @export
#' @docType methods
#' @rdname getGroupB-methods
#' @import methods
setGeneric("getGroupB",					function(object) standardGeneric("getGroupB"));

#' Returns the \code{groupNames} parameter
#' 
#' Returns the \code{groupNames} parameter of the \code{DEMIGroup} object. It is a \code{charater} that
#' represent the custom group names. The \code{groupNames} parameter is only stored when the user defines
#' his/her own clusters.
#'  
#' @param object A DEMIGroup object.
#' @return Returns the \code{groupNames} parameter of the \code{DEMIGroup} object that is a \code{character}.
#' @seealso \code{DEMIGroup}
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
#' # Now we can continue the example of the function getGroupNames
#' 
#' # Set up an experiment
#' demiexp <- DEMIExperiment( analysis = 'gene', celpath = destfolder,
#' 			experiment = 'myexperiment', organism = 'homo_sapiens' )
#' 
#' # Define a custom 'DEMIClust' object with user defined clusters
#' demiclust_custom <- DEMIClust( demiexp, cluster = list( customcluster = c(1190, 1998, 2007) ) )
#' 
#' # Calcuate differential expression
#' demidiff_custom <- DEMIDiff( demiclust_custom )
#' 
#' # Get group B name
#' getGroupNames( getGroup( demiclust_custom ) )
#' getGroupNames( getGroup( demidiff_custom ) )
#' 
#' }
#' 
#' @export
#' @docType methods
#' @rdname getGroupNames-methods
#' @import methods
setGeneric("getGroupNames",				function(object) standardGeneric("getGroupNames"));


#
#	DEMIExperiment
#

#------------------------------------------------------------------------------#
# DEMIExperiment get functions
#------------------------------------------------------------------------------#

#' Returns the \code{analysis} parameter
#' 
#' Returns the \code{analysis} parameter of the \code{DEMIExperiment} object. It is a \code{character}
#' that represents the type of DEMI analysis. It can be either 'gene', 'transcript', 'exon' or 'genome'.
#' 
#' @param object A \code{DEMIExperiment} object.
#' @return Returns the \code{analysis} parameter of the \code{DEMIExperiment} object that is a \code{character}.
#' @seealso \code{DEMIExperiment}
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
#' # Now we can continue the example of the function getAnalysis
#' 
#' # Set up an experiment
#' demiexp <- DEMIExperiment( analysis = 'gene', celpath = destfolder,
#' 			experiment = 'myexperiment', organism = 'homo_sapiens' )
#' 
#' # Retrieve the 'analysis' parameter
#' getAnalysis( demiexp )
#'
#' }
#' 
#' @export
#' @docType methods
#' @rdname getAnalysis-methods
#' @import methods
setGeneric("getAnalysis",				function(object) standardGeneric("getAnalysis"));

#' Returns the \code{celpath} parameter
#' 
#' Returns the \code{celpath} parameter of the \code{DEMIExperiment} object. It is a \code{character} that
#' represents the directory where CEL files are located or is a vector of individual CEL files used in
#' the DEMI analysis.
#' 
#' @param object A \code{DEMIExperiment} object.
#' @return Returns the \code{celpath} parameter of the \code{DEMIExperiment} object that is a \code{character}.
#' @seealso \code{DEMIExperiment}
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
#' # Now we can continue the example of the function getCelPath
#' 
#' # Set up an experiment
#' demiexp <- DEMIExperiment( analysis = 'gene', celpath = destfolder,
#' 			experiment = 'myexperiment', organism = 'homo_sapiens' )
#' 
#' # Retrieve the 'celpath' parameter
#' getCelpath( demiexp )
#'
#' }
#' 
#' @export
#' @docType methods
#' @rdname getCelpath-methods
#' @import methods
setGeneric("getCelpath",				function(object) standardGeneric("getCelpath"));

#' Returns the \code{organism} parameter
#' 
#' Returns the \code{organism} parameter of the \code{DEMIExperiment} object. It is a \code{character} that
#' represents the species used in DEMI analysis.
#' 
#' @param object A \code{DEMIExperiment} object.
#' @return Returns the \code{organism} parameter of the \code{DEMIExperiment} object that is a \code{character}.
#' @seealso \code{DEMIExperiment}
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
#' # Now we can continue the example of the function getOrganism
#' 
#' # Set up an experiment
#' demiexp <- DEMIExperiment( analysis = 'gene', celpath = destfolder,
#' 			experiment = 'myexperiment', organism = 'homo_sapiens' )
#' 
#' # Retrieve the 'organism' parameter
#' getOrganism( demiexp )
#'
#' }
#' 
#' @export
#' @docType methods
#' @rdname getOrganism-methods
#' @import methods
setGeneric("getOrganism",				function(object) standardGeneric("getOrganism"));

#' Returns the \code{arraytype} parameter
#' 
#' Returns the \code{arraytype} parameter of the \code{DEMIExperiment} object. It is a \code{character} that
#' represents the microarray platform used in DEMI analysis.
#' 
#' @param object A \code{DEMIExperiment} object.
#' @return Returns the \code{arraytype} parameter of the \code{DEMIExperiment} object that is a \code{character}.
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
#' # Now we can continue the example of the function getArraytype
#' 
#' # Set up an experiment
#' demiexp <- DEMIExperiment( analysis = 'gene', celpath = destfolder,
#' 			experiment = 'myexperiment', organism = 'homo_sapiens' )
#' 
#' # Retrieve the 'arraytype' parameter
#' getArraytype( demiexp )
#'
#' }
#' 
#' @export
#' @docType methods
#' @rdname getArraytype-methods
#' @import methods
setGeneric("getArraytype",				function(object) standardGeneric("getArraytype"));

#' Returns the \code{maxtargets} parameter
#' 
#' Returns the \code{maxtargets} parameter of the \code{DEMIExperiment} object. It is a \code{numeric} that represents
#' the maximum number of allowed targets (e.g. genes or transcripts) a probe can have a match against.
#' 
#' @param object A \code{DEMIExperiment} object.
#' @return Returns the \code{maxtargets} parameter of the \code{DEMIExperiment} object that is a \code{numeric}.
#' @details 
#' 
#' If the \code{analysis} in \code{DEMIExperiment} object is set to 'transcript' the program still calculates the number
#' of matches on genes. Hence a probe matching two transcripts on the same gene would be included but a probe matching
#' two transcripts on different genes would not be included if \code{maxtargets} would be set to 1.
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
#' # Now we can continue the example of the function getMaxtargets
#' 
#' # Set up an experiment
#' demiexp <- DEMIExperiment( analysis = 'gene', celpath = destfolder,
#' 			experiment = 'myexperiment', organism = 'homo_sapiens' )
#' 
#' # Retrieve the 'maxtargets' parameter
#' getMaxtargets( demiexp )
#'
#' }
#' 
#' @export
#' @docType methods
#' @rdname getMaxtargets-methods
#' @import methods
setGeneric("getMaxtargets",				function(object) standardGeneric("getMaxtargets"));

#' Returns the \code{maxprobes} parameter
#' 
#' Returns the \code{maxprobes} of the \code{DEMIExperiment} object. It is a \code{character} that represents
#' the maximum number of probes a target is allowed to have a match against in DEMI analysis.
#' 
#' @param object A \code{DEMIExperiment} object.
#' @return Returns the \code{maxprobes} parameter of the \code{DEMIExperiment} object.
#' @seealso \code{DEMIExperiment}
#' @details 
#' 
#' If the \code{maxprobes} in \code{DEMIExperiment} object is set to 'median' or some integer larger then 0, then all
#' targets that yield more alignments to different probes then defined by \code{maxprobes} will be scaled down to the
#' number set in the \code{maxprobes} parameter.
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
#' # Now we can continue the example of the function maxprobes
#' 
#' # Set up an experiment
#' demiexp <- DEMIExperiment( analysis = 'gene', celpath = destfolder,
#' 			experiment = 'myexperiment', organism = 'homo_sapiens' )
#' 
#' # Retrieve the 'maxprobes' parameter
#' getMaxprobes( demiexp )
#'
#' }
#' 
#' @export
#' @docType methods
#' @rdname getMaxprobes-methods
#' @import methods
setGeneric("getMaxprobes",				function(object) standardGeneric("getMaxprobes"));

#' Returns the \code{annoTable} parameter representing annotation information
#' 
#' Returns the \code{annoTable} of the \code{DEMIExperiment} object. It is a \code{data.frame} that stores the information
#' about target annotations (such as it's ID's and description) used in DEMI analysis.
#' 
#' @param object A \code{DEMIExperiment} object.
#' @return Returns the \code{annoTable} parameter of the \code{DEMIExperiment} object that is a \code{data.frame}.
#' @seealso \code{DEMIExperiment}
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
#' # Now we can continue the example of the function getAnnotation
#' 
#' # Set up an experiment
#' demiexp <- DEMIExperiment( analysis = 'gene', celpath = destfolder,
#' 			experiment = 'myexperiment', organism = 'homo_sapiens' )
#' 
#' # Retrieve the 'annotation' parameter representing annotation information
#' head( getAnnotation( demiexp ) )
#'
#' }
#' 
#' @export
#' @docType methods
#' @rdname getAnnotation-methods
#' @import methods
setGeneric("getAnnotation",				function(object) standardGeneric("getAnnotation"));

#' Returns the \code{blatTable} parameter representing alignment information
#' 
#' Returns the \code{blatTable} of the \code{DEMIExperiment} object. It is a \code{data.frame} that stores the alignment
#' information (such as probe matches on targets) used in DEMI analysis.
#' 
#' @param object A \code{DEMIExperiment} object.
#' @return Returns the \code{blatTable} parameter of the \code{DEMIExperiment} object that is a \code{data.frame}.
#' @seealso \code{DEMIExperiment}
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
#' # Now we can continue the example of the function getAlignment
#' 
#' # Set up an experiment
#' demiexp <- DEMIExperiment( analysis = 'gene', celpath = destfolder,
#' 			experiment = 'myexperiment', organism = 'homo_sapiens' )
#' 
#' # Retrieve the 'blatTable' parameter representing alignment information
#' head( getAlignment( demiexp ) )
#'
#' }
#' 
#' @export
#' @docType methods
#' @rdname getAlignment-methods
#' @import methods
setGeneric("getAlignment", 				function(object) standardGeneric("getAlignment"));

#' Returns the \code{cytoband} parameter representing karyotype information
#' 
#' Returns the \code{cytoband} parameter of the \code{DEMIExperiment} object. It is a \code{data.frame} that stores
#' the karyotype information of the chromosomes.
#' 
#' @param object A \code{DEMIExperiment} object.
#' @return Returns the \code{cytoband} paramter of the \code{DEMIExperiment} object that is a \code{data.frame}.
#' @seealso \code{DEMIExperiment}
#' @details 
#' 
#' If the \code{analysis} parameter in \code{DEMIExperiment} object is set to 'genome' then genome sections are
#' being annotated by their karyotypes. The annotation information is stored in the \code{cytoband} parameter of
#' the \code{DEMIExperiment} object.
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
#' # Now we can continue the example of the function getCytoband
#' 
#' # Set up an experiment. Note that the cytoband can only retrieved when the analysis
#' # has been set to genome.
#' demiexp_genome <- DEMIExperiment( analysis = 'genome', celpath = destfolder,
#' 		experiment = 'myexperiment', organism = 'homo_sapiens', sectionsize = 500000 )
#' 
#' # Retrieve the 'cytoband' parameter representing karyotype information
#' head( getCytoband( demiexp_genome ) )
#'
#' }
#' 
#' @export
#' @docType methods
#' @rdname getCytoband-methods
#' @import methods
setGeneric("getCytoband", 				function(object) standardGeneric("getCytoband"));

#' Returns the \code{pathway} parameter representing functional annotation information
#' 
#' Returns the \code{pathway} parameter of the \code{DEMIExperiment} object. It is a \code{data.frame} that stores
#' information about gene ontology categories.
#' 
#' @param object A \code{DEMIExperiment} object.
#' @return Returns the \code{pathway} parameter of the \code{DEMIExperiment} object that is a \code{data.frame}.
#' @seealso \code{DEMIExperiment}, \code{DEMIPathway}
#' @details 
#' 
#' The information about gene ontology categories is used when the user runs pathway analysis on DEMI differential
#' expression results with the function \code{DEMIPathway}.
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
#' # Now we can continue the example of the function getPathway. Note that pathway can only
#' # be retrieved if the analysis is set to gene or transcript.
#' 
#' # Set up an experiment
#' demiexp <- DEMIExperiment( analysis = 'gene', celpath = destfolder,
#' 			experiment = 'myexperiment', organism = 'homo_sapiens' )
#' 
#' # Retrieve the 'pathway' parameter representing functional annotation information
#' head( getPathway( demiexp ) )
#'
#' }
#' 
#' @export
#' @docType methods
#' @rdname getPathway-methods
#' @import methods
setGeneric("getPathway", 				function(object) standardGeneric("getPathway"));

#' Returns the raw expression matrix
#' 
#' Returns the raw expression matrix of the \code{DEMIExperiment} object. It is a \code{matrix} where column
#' names indicate different file names and row names indicate probe ID's.
#' 
#' @param object A \code{DEMIExperiment} object.
#' @return Returns the raw expression matrix of the \code{DEMIExperiment} object that is a \code{matrix}.
#' @seealso \code{DEMIExperiment}, \code{DEMICel}
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
#' # Now we can continue the example of the function getCelMatrix
#' 
#' # Set up an experiment
#' demiexp <- DEMIExperiment( analysis = 'gene', celpath = destfolder,
#' 			experiment = 'myexperiment', organism = 'homo_sapiens' )
#' 
#' # Retrieve the raw expression matrix
#' head( getCelMatrix( demiexp ) )
#'
#' }
#' 
#' @export
#' @docType methods
#' @rdname getCelMatrix-methods
#' @import methods
setGeneric("getCelMatrix",				function(object) standardGeneric("getCelMatrix"));

#' Returns the normalized expression matrix
#' 
#' Returns the normalized expression matrix of the \code{DEMIExperiment} object. It is a \code{matrix} where
#' column names indicate different file names and row names indicate probe ID's.
#' 
#' @param object A \code{DEMIExperiment} object.
#' @return Returns the normalized expression matrix of the \code{DEMIExperiment} object that is a \code{matrix}.
#' @seealso \code{DEMIExperiment}, \code{DEMICel}
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
#' # Now we can continue the example of the function getNormMatrix
#' 
#' # Set up an experiment
#' demiexp <- DEMIExperiment( analysis = 'gene', celpath = destfolder,
#' 			experiment = 'myexperiment', organism = 'homo_sapiens' )
#' 
#' # Retrieve the normalized expression matrix
#' head( getNormMatrix( demiexp ) )
#'
#' }
#' 
#' @export
#' @docType methods
#' @rdname getNormMatrix-methods
#' @import methods
setGeneric("getNormMatrix",				function(object) standardGeneric("getNormMatrix"));

#' Returns annotation information for the specified targets
#' 
#' Returns annotation information for the specified targets from a \code{DEMIExperiment} object. Depending on
#' the \code{analysis} parameter in the \code{DEMIExperiment} object the \code{target} parameter can
#' be an ensembl gene ID or gene symbol (e.g. 'MAOB'), ensembl transcript ID, ensembl peptide ID or
#' genomic region ID.
#' 
#' @param object A \code{DEMIExperiment} object.
#' @param target A \code{vector}. A vector of targets whose annotation information should be
#' 		  returned. Depending on the analysis the \code{target} can be ensembl gene ID or gene
#' 		  symbol (e.g. 'MAOB'), ensembl transcript ID, ensembl peptide ID or genomic region ID.
#' @return Returns annotation information from \code{DEMIExperiment} object specified by the targets.
#' @seealso \code{DEMIExperiment}
#' @details 
#' 
#' To see available targets used in the analysis you can try \code{head(getAnnotation(x))} where x is an
#' object of class \code{DEMIExperiment}.
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
#' # Now we can continue the example of the function getTargetAnnotation
#' 
#' # Set up an experiment
#' demiexp <- DEMIExperiment( analysis = 'gene', celpath = destfolder,
#' 			experiment = 'myexperiment', organism = 'homo_sapiens' )
#' 
#' # Retrieve the target annotation for gene symbol 'MAOB' and 'NGB'
#' getTargetAnnotation( demiexp, c( "MAOB", "NGB" ) )
#'
#' }
#' 
#' @export
#' @docType methods
#' @rdname getTargetAnnotation-methods
#' @import methods
setGeneric("getTargetAnnotation",		function(object, target) standardGeneric("getTargetAnnotation"));

#------------------------------------------------------------------------------#
# DEMIExperiment set functions
#------------------------------------------------------------------------------#

#' Attach results from \code{DEMIDiff} object to \code{DEMIExperiment} object
#' 
#' The function \code{attachResult} attaches results stored in a \code{DEMIDiff} object to the underlying
#' \code{DEMIExperiment} object. This function is useful because \code{DEMIDiff} can store results only for
#' one differential expression analysis run whereas \code{DEMIExperiment} object can store all the results
#' done on the same metadata stored in the \code{DEMIExperiment} object. So the user is allowed to keep
#' several DEMI differential expression analysis results in one \code{DEMIExperiment} object for ease of use.
#' 
#' @param object A \code{DEMIExperiment} object. The user needs to make sure that the \code{DEMIExperiment}
#' 		  object where the results will be added is identical to the \code{DEMIExperiment} object whose metadata
#' 		  was used to calculate differential expression.
#' @param diffObject A \code{DEMIDiff} object. The results from the \code{diffObject} parameter will be added to
#' 		  the results of the \code{DEMIExperiment} object in the \code{object} parameter.
#' @return Returns a \code{DEMIExperiment} updated with the results from \code{DEMIDiff} object.
#' @seealso \code{DEMIExperiment},\code{DEMIDiff},\code{getExperiment},\code{identical}
#' @details 
#' 
#' When adding results to \code{DEMIExperiment} object from a \code{DEMIDiff} object the user needs to make sure that
#' the \code{DEMIExperiment} object that is stored under \code{DEMIDiff} object is identical to the \code{DEMIExperiment}
#' object where the results will be added to. You can access the \code{DEMIExperiment} object from the \code{DEMIDiff}
#' object with the function \code{getExperiment(x)} where x is a \code{DEMIDiff} object. With the function \code{identical}
#' you can check if the \code{DEMIExperiment} objects are indeed identical.
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
#' # Now we can continue the example of the function attachResult.
#' 
#' # Set up an experiment
#' demiexp <- DEMIExperiment( analysis = 'gene', celpath = destfolder,
#' 			experiment = 'myexperiment', organism = 'homo_sapiens' )
#' 
#' # Create clusters with an optimized wilcoxon's rank sum test incorporated within demi that
#' # precalculates the probabilities
#' demiclust <- DEMIClust( demiexp, group = c( "BRAIN", "UHR" ), clust.method = demi.wilcox.test.fast )
#' 
#' # Calcuate differential expression
#' demidiff <- DEMIDiff( demiclust )
#' 
#' # Attach the differential expression analysis results to the original 'DEMIExperiment' object
#' demiexp <- attachResult( demiexp, demidiff )
#'
#' }
#' 
#' @export
#' @docType methods
#' @rdname attachResult-methods
#' @import methods
setGeneric("attachResult",				function(object, diffObject) standardGeneric("attachResult"));

#------------------------------------------------------------------------------#
# DEMIExperiment other functions:
#------------------------------------------------------------------------------#

#' Checks if the probes are available
#' 
#' The function \code{check4probe} checks if the probe ID's specified in the probes vector are present in the
#' alignment data of the specified \code{DEMIExperiment} object.
#' 
#' @param object A \code{DEMIExperiment} object.
#' @param probes A \code{vector}. A vector of probe ID's.
#' @return Returns NULL if all the probes are exist in the alignment data, else returns an error message.
#' @seealso \code{DEMIExperiment}
#' @details 
#' 
#' To see which probes are available in the alignment data use the function getAlignment(x)$probeID where
#' x is an object of class \code{DEMIExperiment}.
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
#' # Now we can continue the example of the function check4probe
#' 
#' # Set up an experiment
#' demiexp <- DEMIExperiment( analysis = 'gene', celpath = destfolder,
#' 			experiment = 'myexperiment', organism = 'homo_sapiens' )
#' 
#' # Create clusters with an optimized wilcoxon's rank sum test incorporated within demi that
#' # precalculates the probabilities
#' demiclust <- DEMIClust( demiexp, group = c( "BRAIN", "UHR" ), clust.method = demi.wilcox.test.fast )
#' 
#' # Check for probe ID's
#' check4probe( demiexp, c( 1155955, 100210 ) )
#' 
#' # To see what probes gave the error
#' setdiff( c( 1155955, 100210 ), getAlignment( demiexp )$probeID )
#' 
#' }
#' 
#' @export
#' @docType methods
#' @rdname check4probe-methods
#' @import methods
setGeneric("check4probe",				function(object, probes) standardGeneric("check4probe"));

#' Checks if the targets are available
#' 
#' The function \code{check4target} checks if the targets specified in the target vector are present in the
#' alignment data of the specified \code{DEMIExperiment} object.
#' 
#' @param object A \code{DEMIExperiment} object.
#' @param target A \code{vector}. Depending on the analysis the \code{target} can be an ensembl gene
#' 		  ID or gene symbol (e.g. 'MAOB'), ensembl transcript ID, ensembl peptide ID or genomic region ID.
#' @return Returns TRUE if all the targets exists, else stops with an error message.
#' @seealso \code{DEMIExperiment}
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
#' # Now we can continue the example of the function check4target
#' 
#' # Set up an experiment
#' demiexp <- DEMIExperiment( analysis = 'gene', celpath = destfolder,
#' 			experiment = 'myexperiment', organism = 'homo_sapiens' )
#' 
#' # Create clusters with an optimized wilcoxon's rank sum test incorporated within demi that
#' # precalculates the probabilities
#' demiclust <- DEMIClust( demiexp, group = c( "BRAIN", "UHR" ), clust.method = demi.wilcox.test.fast )
#' 
#' # Check the target for gene symbols 'MAOB', 'EMPTY' and 'NGB'
#' check4target( demiexp, c( "MAOB", "EMPTY", "NGB" ) )
#' 
#' # Since 'EMPTY' gave an error try without it
#' check4target( demiexp, c( "MAOB", "NGB" ) )
#' 
#' }
#' 
#' @export
#' @docType methods
#' @rdname check4target-methods
#' @import methods
setGeneric("check4target",				function(object, target) standardGeneric("check4target"));

#' Initializes the normalization of the raw expression matrix
#' 
#' The function \code{celMatrixNormalize} initializes the normalization of the raw expression matrix
#' in the \code{DEMIExperiment} object. It is used internally in DEMI analysis.
#' 
#' @param object A \code{DEMIExperiment} object. It stores the raw expression matrix.
#' @param fun A \code{function}. The function used to normalize the raw expression matrix.
#' @return Returns a \code{DEMIExperiment} object updated with normalized expression matrix.
#' 
#' @author Sten Ilmjarv
#' 
#' @export
#' @docType methods
#' @rdname celMatrixNormalize-methods
#' @import methods
setGeneric("celMatrixNormalize",		function(object, fun) standardGeneric("celMatrixNormalize"));

#' Draws a histogram of the normalized expression levels of the specified targets
#' 
#' The function \code{probe.levels} draws a histogram of the normalized expression levels for
#' the specified targets. Depending on the analysis the \code{target} can be ensembl gene ID
#' or gene symbol (e.g. 'MAOB'), ensembl transcript ID, ensembl peptide ID or genomic region
#' ID.
#' 
#' @param object A \code{DEMIExperiment} object.
#' @param target A \code{vector}. Depending on the analysis the \code{target} can be ensembl gene
#' 		  ID or gene symbol (e.g. 'MAOB'), ensembl transcript ID, ensembl peptide ID or genomic
#' 		  region ID.
#' @return A \code{ggplot} object.
#' @seealso \code{DEMIExperiment}
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
#' # Now we can continue the example of the function probe.level
#' 
#' # Set up an experiment
#' demiexp <- DEMIExperiment( analysis = 'gene', celpath = destfolder,
#' 			experiment = 'myexperiment', organism = 'homo_sapiens' )
#' 
#' # Draw probes levels measuring the gene 'MAOB'
#' pdf( "MAOB_probe_levels.pdf", width=8, height=8 )
#' probe.levels( demiexp, "MAOB" )
#' dev.off() 
#' 
#' }
#' 
#' @export
#' @docType methods
#' @rdname probe.levels-methods
#' @import methods
setGeneric("probe.levels",				function(object, target) standardGeneric("probe.levels"));

#' Draws a plot of the normalized expression levels of the specified targets
#' 
#' The function \code{probe.plot} draws a plot of the normalized expression levels for the
#' specified targets. Depending on the analysis the \code{target} can be ensembl gene ID
#' or gene symbol (e.g. 'MAOB'), ensembl transcript ID, ensembl peptide ID or genomic region
#' ID.
#' 
#' @param object A \code{DEMIExperiment} object.
#' @param target A \code{vector}. Depending on the analysis the \code{target} can be ensembl gene
#' 		  ID or gene symbol (e.g. 'MAOB'), ensembl transcript ID, ensembl peptide ID or genomic region
#' 		  ID.
#' @return A \code{ggplot} object.
#' @seealso \code{DEMIExperiment}
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
#' # Now we can continue the example of the function probe.plot
#' 
#' # Set up an experiment
#' demiexp <- DEMIExperiment( analysis = 'gene', celpath = destfolder,
#' 			experiment = 'myexperiment', organism = 'homo_sapiens' )
#' 
#' # Draw probes levels measuring the gene 'MAOB'
#' pdf( "MAOB_probe_plot.pdf", width=8, height=8 )
#' probe.plot( demiexp, "MAOB" )
#' dev.off() 
#' 
#' }
#' 
#' @export
#' @docType methods
#' @rdname probe.plot-methods
#' @import methods
setGeneric("probe.plot",				function(object, target) standardGeneric("probe.plot"));

#------------------------------------------------------------------------------#
# DEMIExperiment load functions:
#------------------------------------------------------------------------------#

#' Loads the DEMI annotation package specified by the \code{DEMIExperiment} object
#' 
#' The function \code{loadDEMILibrary} loads the DEMI annotation packages specified by the \code{organism} parameter
#' in the \code{DEMIExperiment} object. It is used internally in DEMI analysis.
#' 
#' @param object A \code{DEMIExperiment} object.
#' @return Returns a \code{DEMIExperiment} object updated with annotation and alignment information for the
#' 		   specified microarray platform and species. If the \code{analysis} parameter of the \code{DEMIExperiment}
#' 		   object is set to 'genome' it also attaches the cytoband information and if the \code{analysis} parameter
#' 		   of the \code{DEMIExperiment} object is set to 'gene' or 'transcript' it additionally loads the pathway
#' 		   information.
#' @seealso \code{DEMIExperiment}
#' 
#' @author Sten Ilmjarv
#' 
#' @export
#' @docType methods
#' @rdname loadDEMILibrary-methods
#' @import methods
setGeneric("loadDEMILibrary",			function(object) standardGeneric("loadDEMILibrary"));

#' Loads the annotation information specified by the \code{DEMIExperiment} object
#' 
#' The function \code{loadAnnotation} loads the annotation information for the specified \code{DEMIExperiment}
#' object. It is used internally in DEMI analysis.
#' 
#' @param object A \code{DEMIExperiment} object.
#' @param pkg An \code{environment}. Specifies the environment where to load the data from.
#' @return Returns a \code{data.frame} with annotation information.
#' @seealso \code{DEMIExperiment}, \code{environment}
#' 
#' @author Sten Ilmjarv
#' 
#' @export
#' @docType methods
#' @rdname loadAnnotation-methods
#' @import methods
setGeneric("loadAnnotation",			function(object, pkg) standardGeneric("loadAnnotation"));

#' Loads the alignment information specified by the \code{DEMIExperiment} object
#' 
#' The function \code{loadBlat} loads the alignment information for the specified \code{DEMIExperiment} object. It
#' is used internally in DEMI analysis.
#' 
#' @param object A \code{DEMIExperiment} object.
#' @param pkg An \code{environment}. Specifies the environment where to load the data from.
#' @return Returns a \code{data.frame} with alignment information.
#' @seealso \code{DEMIExperiment}, \code{environment}
#' 
#' @author Sten Ilmjarv
#' 
#' @export
#' @docType methods
#' @rdname loadBlat-methods
#' @import methods
setGeneric("loadBlat",					function(object, pkg) standardGeneric("loadBlat"));

#' Loads the karyotype information specified by the \code{DEMIExperiment} object
#' 
#' The function \code{loadCytoband} loads the karyotype information for the specified \code{DEMIExperiment} object. It
#' is used internally in DEMI analysis.
#' 
#' @param object A \code{DEMIExperiment} object.
#' @param pkg An \code{environment}. Specifies the environment where to load the data from.
#' @return Returns a \code{data.frame} with karyotype information.
#' @seealso \code{DEMIExperiment}, \code{environment}
#' 
#' @author Sten Ilmjarv
#' 
#' @export
#' @docType methods
#' @rdname loadCytoband-methods
#' @import methods
setGeneric("loadCytoband",				function(object, pkg) standardGeneric("loadCytoband"));

#' Loads the pathway information specified by the \code{DEMIExperiment} object
#' 
#' the function \code{loadPathway} loads the pathway information for the specified \code{DEMIExperiment} object. It
#' is used internally in DEMI analysis.
#' 
#' @param object A \code{DEMIExperiment} object.
#' @param pkg An \code{environment}. Specifies the environment where to load the data from.
#' @return Returns a \code{data.frame} with pathway information.
#' @seealso \code{DEMIExperiment}, \code{environment}
#' 
#' @author Sten Ilmjarv
#' 
#' @export
#' @docType methods
#' @rdname loadPathway-methods
#' @import methods
setGeneric("loadPathway",				function(object, pkg) standardGeneric("loadPathway"));

#' Loads the raw expression matrix into a \code{DEMIExperiment} object
#' 
#' The function \code{loadCel} loads the raw expression matrix from CEL files into a \code{DEMIExperiment}
#' object. It is used internally in DEMI analysis.
#' 
#' @param object A \code{DEMIExperiment} object.
#' @return Returns a \code{DEMIExperiment} object updated with a \code{DEMICel} object attached to the slot
#' 		   \code{exprsData} that contains the raw expression matrix loaded from the CEL files.
#' @seealso \code{DEMIExperiment}, \code{DEMICel}
#' 
#' @author Sten Ilmjarv
#' 
#' @export
#' @docType methods
#' @rdname loadCel-methods
#' @import methods
setGeneric("loadCel",					function(object) standardGeneric("loadCel"));

#------------------------------------------------------------------------------#
# DEMIExperiment normalization functions:
#------------------------------------------------------------------------------#

#' Relative rank normalization function
#' 
#' The function \code{norm.rank} normalizes the raw expression matrix by relative ranking. It
#' is used internally in DEMI analysis.
#'
#' @param object A \code{matrix} or \code{numeric}. The raw expression matrix or a single expression vector.
#' @return A \code{data.frame} representing the normalized expression matrix.
#' 
#' @author Sten Ilmjarv
#' @examples 
#' \dontrun{
#' 
#' # Create a matrix with 1000 values that represents raw expression values
#' rawmatrix <- matrix(rexp(1000, rate=1), ncol=8)
#' 
#' # Normalize the raw expression matrix
#' normmatrix <- norm.rrank( rawmatrix )
#' 
#' }
#' 
#' @export
#' @docType methods
#' @rdname norm.rrank-methods
#' @import methods
setGeneric("norm.rrank",				function(object) standardGeneric("norm.rrank"));

#' Quantile normalization function
#' 
#' A function for normalizing the expression matrix with quantiles. In the current It tries to
#' mimic rma quantile normalization. In the current state it is not used in DEMI analysis.
#'
#' @param object A \code{matrix}. The raw expression matrix.
#' @return A \code{data.frame} representing the normalized expression matrix.
#' 
#' @author Sten Ilmjarv
#' @examples 
#' \dontrun{
#' 
#' # Create a matrix with 1000 values that represents raw expression values
#' rawmatrix <- matrix(rexp(1000, rate=1), ncol=8)
#' 
#' # Normalize the raw expression matrix
#' normmatrix <- norm.quantile( rawmatrix )
#' 
#' }
#' 
#' @docType methods
#' @rdname norm.quantile-methods
#' @import methods
setGeneric("norm.quantile",				function(object) standardGeneric("norm.quantile"));


#
#	DEMIClust
#

#------------------------------------------------------------------------------#
# DEMIClust blat grouping functions:
#------------------------------------------------------------------------------#

#' Creates a \code{DEMIGroup} object
#' 
#' The function \code{createGroup} creates a \code{DEMIGroup} object for \code{DEMIClust}
#' object. This function is used internally by DEMI analysis when clusters are created
#' automatically from normalized expression matrix.
#' 
#' @param object A \code{DEMIClust} object.
#' @return Returns a \code{DEMIClust} object that includes a \code{DEMIGroup} objec as the
#' 		   \code{group} parameter of the object.
#' @seealso \code{DEMIClust}, \code{DEMIGroup}
#' 
#' @author Sten Ilmjarv
#' 
#' @export
#' @docType methods
#' @rdname createGroup-methods
#' @import methods
setGeneric("createGroup",				function(object) standardGeneric("createGroup"));

#------------------------------------------------------------------------------#
# DEMIClust get functions:
#------------------------------------------------------------------------------#

#' Returns the \code{clust.method} parameter
#' 
#' Returns the \code{clust.method} parameter of the \code{DEMIClust} object. It is a
#' function that is used for clustering the probes into clusters.
#' 
#' @param object A \code{DEMIClust} object.
#' @return Returns the \code{clust.method} parameter of the \code{DEMIClust} object that is a \code{function}.
#' @seealso \code{DEMIClust}
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
#' # Now we can continue the example of the function getClustMethod
#' 
#' # Set up an experiment
#' demiexp <- DEMIExperiment( analysis = 'gene', celpath = destfolder,
#' 			experiment = 'myexperiment', organism = 'homo_sapiens' )
#' 
#' # Create clusters with an optimized wilcoxon's rank sum test incorporated within demi that
#' # precalculates the probabilities
#' demiclust <- DEMIClust( demiexp, group = c( "BRAIN", "UHR" ), clust.method = demi.wilcox.test.fast )
#' 
#' # Retrieve the 'clust.method' parameter that is a 'function'
#' getClustMethod( demiclust )
#'
#' }
#' 
#' @export
#' @docType methods
#' @rdname getClustMethod-methods
#' @import methods
setGeneric("getClustMethod",			function(object) standardGeneric("getClustMethod"));

#' Returns the \code{cutoff.pvalue} parameter
#' 
#' Returns the \code{cutoff.pvalue} parameter of the 'DEMIClust' object. It is used to
#' determine the cutoff significance level of probe signalling when applying the
#' clustering function.
#' 
#' @param object A \code{DEMIClust} object.
#' @return Returns the \code{cutoff.pvalue} parameter of the 'DEMIClust' object.
#' @seealso \code{DEMIClust}
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
#' # Now we can continue the example of the function getCutoffPvalue
#' 
#' # Set up an experiment
#' demiexp <- DEMIExperiment( analysis = 'gene', celpath = destfolder,
#' 			experiment = 'myexperiment', organism = 'homo_sapiens' )
#' 
#' # Create clusters with an optimized wilcoxon's rank sum test incorporated within demi that
#' # precalculates the probabilities
#' demiclust <- DEMIClust( demiexp, group = c( "BRAIN", "UHR" ), clust.method = demi.wilcox.test.fast )
#' 
#' # Retrieve the 'cutoff.pvalue' parameter
#' getCutoffPvalue( demiclust )
#'
#' }
#' 
#' @export
#' @docType methods
#' @rdname getCutoffPvalue-methods
#' @import methods
setGeneric("getCutoffPvalue",			function(object) standardGeneric("getCutoffPvalue"));

#' Returns the \code{cluster} parameter
#' 
#' Returns the \code{cluster} parameter of the \code{DEMIClust} object. It is a list that represents
#' clusters containing probes. 
#' 
#' @param object A \code{DEMIClust} object.
#' @return Returns \code{cutoff.pvalue} parameter of the \code{DEMIClust} object.
#' @seealso \code{DEMIClust}
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
#' # Now we can continue the example of the function getCluster
#' 
#' # Set up an experiment
#' demiexp <- DEMIExperiment( analysis = 'gene', celpath = destfolder,
#' 			experiment = 'myexperiment', organism = 'homo_sapiens' )
#' 
#' # Create clusters with an optimized wilcoxon's rank sum test incorporated within demi that
#' # precalculates the probabilities
#' demiclust <- DEMIClust( demiexp, group = c( "BRAIN", "UHR" ), clust.method = demi.wilcox.test.fast )
#' 
#' # Retrieve the 'cluster' parameter
#' getCluster( demiclust )
#'
#' }
#' 
#' @export
#' @docType methods
#' @rdname getCluster-methods
#' @import methods
setGeneric("getCluster",				function(object) standardGeneric("getCluster"));

#------------------------------------------------------------------------------#
# DEMIClust cluster functions:
#------------------------------------------------------------------------------#

#' Initializes the clustering of probes into clusters
#' 
#' The function \code{cluster} clusters probes with the function specified in the
#' in the \code{clust.method} parameter of the \code{DEMIClust} object. This function is
#' used internally by DEMI analysis when clusters are created automatically from
#' normalized expression matrix.
#' 
#' @param object A \code{DEMIClust} object.
#' @return Returns a \code{DEMIClust} object that is updated with a list of clustered
#' 		   probes.
#' @seealso \code{DEMIClust}
#' 
#' @author Sten Ilmjarv
#' 
#' @export
#' @docType methods
#' @rdname cluster-methods
#' @import methods
setGeneric("cluster",					function(object) standardGeneric("cluster"));

#------------------------------------------------------------------------------#
# DEMIClust other functions:
#------------------------------------------------------------------------------#

#' Checks if the \code{DEMIClust} object is user defined or automatically generated
#' 
#' The function \code{customObject} determines if the \code{DEMIClust} is a custom object
#' defined by the users clusters or built automatically by a clustering method. It is used
#' internally in DEMI analysis.
#' 
#' @param object A \code{DEMIClust} object.
#' @return Returns FALSE if the \code{DEMIClust} object was built automatically. Returns
#' 		   TRUE if the \code{DEMIClust} is user defined.
#' @seealso \code{DEMIClust}
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
#' # Now we can continue the example of the function customObject
#' 
#' # Set up an experiment
#' demiexp <- DEMIExperiment( analysis = 'gene', celpath = destfolder,
#' 			experiment = 'myexperiment', organism = 'homo_sapiens' )
#' 
#' # Create clusters with an optimized wilcoxon's rank sum test incorporated within demi that
#' # precalculates the probabilities
#' demiclust <- DEMIClust( demiexp, group = c( "BRAIN", "UHR" ), clust.method = demi.wilcox.test.fast )
#' 
#' # Check if user defined 'DEMIClust' object
#' customObject( demiclust )
#' 
#' # Define a custom 'DEMIClust' object with user defined clusters
#' demiclust_custom <- DEMIClust( demiexp, cluster = list( customcluster = c(1190, 1998, 2007) ) )
#' 
#' # Check if user defined 'DEMIClust' object
#' customObject( demiclust_custom )
#'
#' }
#' 
#' @export
#' @docType methods
#' @rdname customObject-methods
#' @import methods
setGeneric("customObject",				function(object) standardGeneric("customObject"));


#
#	DEMIDiff
#

#------------------------------------------------------------------------------#
# DEMIDiff get functions:
#------------------------------------------------------------------------------#

#' Returns the \code{cluster} parameter
#' 
#' Returns the \code{cluster} of the \code{DEMIDiff} object. It is a \code{DEMIClust} object.
#' 
#' @param object A \code{DEMIDiff} object.
#' @return Returns the \code{cluster} parameter of the \code{DEMIDiff} object which is
#' 		   a \code{DEMIClust} object
#' @seealso \code{DEMIClust}, \code{DEMIDiff}
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
#' # Now we can continue the example of the function getDEMIClust
#' 
#' # Set up an experiment
#' demiexp <- DEMIExperiment( analysis = 'gene', celpath = destfolder,
#' 			experiment = 'myexperiment', organism = 'homo_sapiens' )
#' 
#' # Create clusters with an optimized wilcoxon's rank sum test incorporated within demi that
#' # precalculates the probabilities
#' demiclust <- DEMIClust( demiexp, group = c( "BRAIN", "UHR" ), clust.method = demi.wilcox.test.fast )
#' 
#' # Calcuate differential expression
#' demidiff <- DEMIDiff( demiclust )
#' 
#' # Retrieve the 'cluster' parameter that is 'DEMIClust' object
#' getDEMIClust( demidiff )
#'
#' }
#' 
#' @export
#' @docType methods
#' @rdname getDEMIClust-methods
#' @import methods
setGeneric("getDEMIClust",				function(object) standardGeneric("getDEMIClust"));

#' Returns the \code{name} parameter
#' 
#' Returns the \code{name} parameter of the \code{DEMIDiff} object. It is a \code{character}
#' that represents the name of the differential expression analysis stored in the \code{DEMIDiff}
#' object.
#' 
#' @param object A \code{DEMIDiff} object.
#' @return Returns the \code{name} parameter of the \code{DEMIDiff} object which is a
#' 		   \code{character}.
#' @seealso \code{DEMIDiff}
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
#' # Now we can continue the example of the function getName
#' 
#' # Set up an experiment
#' demiexp <- DEMIExperiment( analysis = 'gene', celpath = destfolder,
#' 			experiment = 'myexperiment', organism = 'homo_sapiens' )
#' 
#' # Create clusters with an optimized wilcoxon's rank sum test incorporated within demi that
#' # precalculates the probabilities
#' demiclust <- DEMIClust( demiexp, group = c( "BRAIN", "UHR" ), clust.method = demi.wilcox.test.fast )
#' 
#' # Calcuate differential expression
#' demidiff <- DEMIDiff( demiclust )
#' 
#' # Retrieve the 'cluster' parameter that is 'DEMIClust' object
#' getName( demidiff )
#'
#' }
#' 
#' @export
#' @docType methods
#' @rdname getName-methods
#' @import methods
setGeneric("getName",					function(object) standardGeneric("getName"));


#------------------------------------------------------------------------------#
# DEMIDiff diffexp functions:
#------------------------------------------------------------------------------#

#' Initializes the differential expression analysis
#' 
#' The function \code{diffexp} performs differential expression analysis. This function
#' is used internally by DEMI analysis.
#' 
#' @param object A \code{DEMIDiff} object.
#' @return Returns the \code{DEMIDiff} object that is updated with the results from differential
#' 		   expression analysis.
#' @seealso \code{DEMIDiff}
#' 
#' @author Sten Ilmjarv
#' 
#' @export
#' @docType methods
#' @rdname diffexp-methods
#' @import methods
setGeneric("diffexp",					function(object) standardGeneric("diffexp"));


#
#	Multiple Classes
#

#' Returns the \code{group} parameter
#' 
#' Returns the \code{group} parameter which is a \code{DEMIGroup} object.
#' 
#' @param object A \code{DEMIClust}, \code{DEMIDiff} or \code{DEMIResult} object.
#' @return Returns the \code{group} parameter that is a \code{DEMIGroup} object.
#' @seealso \code{DEMIClust}, \code{DEMIDiff}, \code{DEMIResult}
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
#' # Now we can continue the example of the function getGroup
#' 
#' # Set up an experiment
#' demiexp <- DEMIExperiment( analysis = 'gene', celpath = destfolder,
#' 			experiment = 'myexperiment', organism = 'homo_sapiens' )
#' 
#' # Create clusters with an optimized wilcoxon's rank sum test incorporated within demi that
#' # precalculates the probabilities
#' demiclust <- DEMIClust( demiexp, group = c( "BRAIN", "UHR" ), clust.method = demi.wilcox.test.fast )
#' 
#' # Calcuate differential expression
#' demidiff <- DEMIDiff( demiclust )
#' 
#' # Retrieve results from 'DEMIDiff' object that returns 'DEMIResult' object
#' demiresult <- getResult( demidiff )
#' 
#' # Retrieve the 'group' parameter that is 'DEMIGroup' object
#' getGroup( demiclust )
#' getGroup( demidiff )
#' getGroup( demiresult )
#' 
#' }
#' 
#' @export
#' @docType methods
#' @rdname getGroup-methods
#' @import methods
setGeneric("getGroup",	 				function(object) standardGeneric("getGroup"));

#' Returns the \code{experiment} parameter
#' 
#' Returns the \code{experiment} parameter of the specified object. For object of class \code{DEMIExperiment}
#' it returns the name given to the experiment. For objects of class \code{DEMIClust} or \code{DEMIDiff}
#' it returns the initial \code{DEMIExperiment} object. This function can be useful if the user wants
#' to access metadata, such as annotations and alignments from other DEMI analysis objects. As well
#' as accessing the name of the analysis.
#' 
#' @param object A \code{DEMIExperiment}, \code{DEMIClust} or \code{DEMIDiff} object.
#' @return Returns the \code{experiment} parameter. If the input object is \code{DEMIExperiment} it returns a
#' 		   character, if the input object is either \code{DEMIClust} or \code{DEMIDiff} it returns a
#' 		   \code{DEMIExperiment} object.
#' @seealso \code{DEMIExperiment}, \code{DEMIClust}, \code{DEMIDiff}
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
#' # Now we can continue the example of the function getExperiment
#' 
#' # Set up an experiment
#' demiexp <- DEMIExperiment( analysis = 'gene', celpath = destfolder,
#' 			experiment = 'myexperiment', organism = 'homo_sapiens' )
#' 
#' # Create clusters with an optimized wilcoxon's rank sum test incorporated within demi that
#' # precalculates the probabilities
#' demiclust <- DEMIClust( demiexp, group = c( "BRAIN", "UHR" ), clust.method = demi.wilcox.test.fast )
#' 
#' # Calcuate differential expression
#' demidiff <- DEMIDiff( demiclust )
#' 
#' # Retrieve the 'experiment' parameter
#' getExperiment( demiexp )
#' getExperiment( demiclust )
#' getExperiment( demidiff )
#' 
#' }
#' 
#' @export
#' @docType methods
#' @rdname getExperiment-methods
#' @import methods
setGeneric("getExperiment",				function(object) standardGeneric("getExperiment"));

#' Returns the \code{result} parameter
#' 
#' Returns the \code{result} parameter stored in the specified object. If the object is of class
#' \code{DEMIExperiment} then it returns a \code{list} of \code{DEMIResult} objects. If the
#' object is of class \code{DEMIDiff} then it returns only one \code{DEMIResult} object. But if
#' the object is of class \code{DEMIResult} then the function returns a \code{list} that
#' contains the results for every cluster in a \code{data.frame}. However instead of using this
#' function it maybe easier to use the function \code{getResultTable} that returns the \code{result}
#' parameter as a \code{data.frame}.
#' 
#' @param object A \code{DEMIExperiment}, \code{DEMIDiff} or \code{DEMIResult} object.
#' @return Returns the \code{result} parameter. For objects of class \code{DEMIExperiment} it returns
#' 		   a \code{list} of \code{DEMIResult} objects. For objects of class \code{DEMIDiff} it returns
#' 		   a single \code{DEMIResult} object and for objects of class \code{DEMIResult} it returns
#' 		   a \code{list}.
#' @seealso \code{DEMIExperiment}, \code{DEMIDiff}, \code{DEMIResult}, \code{getResultTable}
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
#' # Now we can continue the example of the function getResult
#' 
#' # Set up an experiment
#' demiexp <- DEMIExperiment( analysis = 'gene', celpath = destfolder,
#' 			experiment = 'myexperiment', organism = 'homo_sapiens' )
#' 
#' # Create clusters with an optimized wilcoxon's rank sum test incorporated within demi that
#' # precalculates the probabilities
#' demiclust <- DEMIClust( demiexp, group = c( "BRAIN", "UHR" ), clust.method = demi.wilcox.test.fast )
#' 
#' # Calcuate differential expression
#' demidiff <- DEMIDiff( demiclust )
#' 
#' # Attach results from 'DEMIDiff' object to 'DEMIExperiment' object
#' demiexp_attached <- attachResult( demiexp, demidiff )
#' 
#' # Retrieve the 'result' parameter
#' getResult( demiexp_attached )
#' getResult( demidiff )
#' 
#' }
#' 
#' @export
#' @docType methods
#' @rdname getResult-methods
#' @import methods
setGeneric("getResult",					function(object) standardGeneric("getResult"));

#' Retruns the DEMI analysis results as a \code{data.frame}
#' 
#' The function \code{getResultTable} returns the DEMI analysis results as a \code{data.frame}. It retrieves
#' the \code{result} parameter of the specified object with the function \code{getResult} and converts it into
#' a \code{data.frame} for convenient viewing.
#' 
#' @param object A \code{DEMIExperiment} or \code{DEMIDiff} object.
#' @return Returns the \code{result} parameter of the specified object as a \code{data.frame}.
#' @seealso \code{DEMIExperiment}, \code{DEMIDiff}, \code{getResult}
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
#' # Now we can continue the example of the function getResultTable
#' 
#' # Set up an experiment
#' demiexp <- DEMIExperiment( analysis = 'gene', celpath = destfolder,
#' 			experiment = 'myexperiment', organism = 'homo_sapiens' )
#' 
#' # Create clusters with an optimized wilcoxon's rank sum test incorporated within demi that
#' # precalculates the probabilities
#' demiclust <- DEMIClust( demiexp, group = c( "BRAIN", "UHR" ), clust.method = demi.wilcox.test.fast )
#' 
#' # Calcuate differential expression
#' demidiff <- DEMIDiff( demiclust )
#' 
#' # Attach results from 'DEMIDiff' object to 'DEMIExperiment' object
#' demiexp_attached <- attachResult( demiexp, demidiff )
#' 
#' # Retrieve the 'result' parameter as a 'data.frame'
#' head( getResultTable( demiexp_attached ) )
#' head( getResultTable( demidiff ) )
#' 
#' }
#' 
#' @export
#' @docType methods
#' @rdname getResultTable-methods
#' @import methods
setGeneric("getResultTable",			function(object) standardGeneric("getResultTable"));

#' Returns the probe levels from the normalized expression matrix for the specified probes
#' 
#' The function \code{getProbeLevel} returns the probe levels in the normalized expression matrix specified
#' by the probe ID's.
#' 
#' @param object A \code{DEMIExperiment} or \code{DEMIDiff} object.
#' @param probes A \code{vector}. A vector of probe ID's whose expression levels should be returned.
#' @param verbose A \code{logical}. If TRUE it will print out the probe ID's that were not found in
#' 		  normalized expression matrix.
#' @return Returns the probe levels in the normalized expression matrix for the specified probes.
#' @seealso \code{DEMIExperiment}, \code{DEMIDiff}
#' @details 
#' 
#' To see what are the available probes in the normalized expression matrix you can try \code{row.names(getNormMatrix(x))}
#' where x is an object of class \code{DEMIExperiment}. 
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
#' # Now we can continue the example of the function getProbeLevel
#' 
#' # Set up an experiment
#' demiexp <- DEMIExperiment( analysis = 'gene', celpath = destfolder,
#' 			experiment = 'myexperiment', organism = 'homo_sapiens' )
#' 
#' # Create clusters with an optimized wilcoxon's rank sum test incorporated within demi that
#' # precalculates the probabilities
#' demiclust <- DEMIClust( demiexp, group = c( "BRAIN", "UHR" ), clust.method = demi.wilcox.test.fast )
#' 
#' # Calcuate differential expression
#' demidiff <- DEMIDiff( demiclust )
#' 
#' # Retrieve the probe levels specified by probe ID's of the normalized expression matrix
#' getProbeLevel( demiexp, c( 1171,1182 ), TRUE )
#' getProbeLevel( demidiff, c( 1171,1182 ), TRUE )
#' 
#' }
#'
#' @export
#' @docType methods
#' @rdname getProbeLevel-methods
#' @import methods
setGeneric("getProbeLevel",				function(object, probes, verbose) standardGeneric("getProbeLevel"));

#' Returns the probe ID's of the specified targets
#' 
#' The function \code{getTargetProbes} returns the probe ID's of the specified targets. Depending on
#' the \code{analysis} parameter in the underlying \code{DEMIExperiment} object the \code{target} parameter can
#' be an ensembl gene ID or gene symbol (e.g. 'MAOB'), ensembl transcript ID, ensembl peptide ID or
#' genomic region ID.
#' 
#' @param object A \code{DEMIExperiment} or \code{DEMIDiff} object.
#' @param target A \code{vector}. A vector of targets whose probe ID's should be returned. Depending on
#' 		  the analysis the \code{target} can be ensembl gene ID or gene symbol (e.g. 'MAOB'), ensembl
#' 		  transcript ID, ensembl peptide ID or genomic region ID.
#' @return Returns the probes ID's specified by the targets.
#' @seealso \code{DEMIExperiment}, \code{DEMIDiff}
#' @details 
#' 
#' To see available targets used in the analysis you can try \code{head(getAnnotation(x))} where x is an
#' object of class \code{DEMIExperiment}. Alternatively you could use \code{head(getAnnotation(getExperiment(y)))}
#' where y is of class \code{DEMIDiff}.
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
#' # Now we can continue the example of the function getTargetProbes
#' 
#' # Set up an experiment
#' demiexp <- DEMIExperiment( analysis = 'gene', celpath = destfolder,
#' 			experiment = 'myexperiment', organism = 'homo_sapiens' )
#' 
#' # Create clusters with an optimized wilcoxon's rank sum test incorporated within demi that
#' # precalculates the probabilities
#' demiclust <- DEMIClust( demiexp, group = c( "BRAIN", "UHR" ), clust.method = demi.wilcox.test.fast )
#' 
#' # Calcuate differential expression
#' demidiff <- DEMIDiff( demiclust )
#' 
#' # Retrieve the probe ID's of the specified targets
#' getTargetProbes( demiexp, "MAOB" )
#' getTargetProbes( demidiff, "MAOB" )
#' 
#' }
#' 
#' @export
#' @docType methods
#' @rdname getTargetProbes-methods
#' @import methods
setGeneric("getTargetProbes",			function(object, target) standardGeneric("getTargetProbes"));

#' Returns the mean normalized expression levels for the specified targets
#' 
#' The function \code{demisummary} returns the mean normalized expression levels for the
#' specified targets. It returns the mean expression values for the whole dataset as well
#' as for individual groups. Depending on the \code{analysis} parameter of the underlying
#' \code{DEMIExperiment} object the \code{target} can be ensembl gene ID or gene symbol
#' (e.g. 'MAOB'), ensembl transcript ID, ensembl peptide ID or genomic region ID.
#' 
#' @param object A \code{DEMIExperiment}, \code{DEMIDiff} object.
#' @param target A \code{vector}. Depending on the analysis the \code{target} can be ensembl gene ID or
#' 		  gene symbol (e.g. 'MAOB'), ensembl transcript ID, ensembl peptide ID or genomic region ID.
#' @return Returns the mean normalized expression levels of the specified targets.
#' @seealso \code{DEMIExperiment},\code{DEMIDiff},\code{attachResult}
#' @details 
#' 
#' To see available targets used in the analysis you can try \code{head(getAnnotation(x))} where x is an
#' object of class \code{DEMIExperiment}. Alternatively you could use \code{head(getAnnotation(getExperiment(y)))}
#' where y is of class \code{DEMIDiff}.
#' 
#' If no results have been attached to the \code{DEMIExperiment} object then it only returns the mean normalized
#' expression values for the whole dataset not for individual groups. To attach results to \code{DEMIExperiment}
#' object use the function \code{attachResult(x,y)} where x is an object of class \code{DEMIExperiment} and y is
#' an object of class \code{DEMIDiff} that stores the results.  
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
#' # Now we can continue the example of the function demisummary
#' 
#' # Set up an experiment
#' demiexp <- DEMIExperiment( analysis = 'gene', celpath = destfolder,
#' 			experiment = 'myexperiment', organism = 'homo_sapiens' )
#' 
#' # Create clusters with an optimized wilcoxon's rank sum test incorporated within demi that
#' # precalculates the probabilities
#' demiclust <- DEMIClust( demiexp, group = c( "BRAIN", "UHR" ), clust.method = demi.wilcox.test.fast )
#' 
#' # Calcuate differential expression
#' demidiff <- DEMIDiff( demiclust )
#' 
#' # Retrieve the mean normalized expression values for the specified targets
#' demisummary( demiexp, c( "MAOB" ) )
#' demisummary( demidiff, "MAOB" )
#' 
#' # Attach results from 'DEMIDiff' object to 'DEMIExperiment' object
#' demiexp_attached <- attachResult( demiexp, demidiff )
#' 
#' # Retrieve mean normalized expression values again and note these are also retrieved for specific
#' # groups
#' demisummary( demiexp_attached, "MAOB" )
#' 
#' }
#' 
#' @export
#' @docType methods
#' @rdname demisummary-methods
#' @import methods
setGeneric("demisummary",				function(object, target) standardGeneric("demisummary"));
