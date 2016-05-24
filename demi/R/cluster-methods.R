#==============================================================================#
# cluster-methods.R:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# wprob
# demiequal
# demi.wilcox.test.fast
# demi.wilcox.test
# demi.t.test
# demi.comp.test
#==============================================================================#


#' Calculates wilcoxon's upper and lower probabilities
#' 
#' Calculates the wilcoxon's upper and lower probabilities for each
#' possible rank sum defined by the size of the test and reference samples.
#' 
#' @param m An \code{integer}. Defines the test sample size.
#' @param n An \code{integer}. Defines the reference sample size.
#' @return A \code{list}. Returns a \code{list} of all possible lower and upper tail
#' 		   p-values defined by the sum of the possible rank combinations.
#' 
#' @author Sten Ilmjarv
#' @examples 
#' 
#' # For test sample 4 and reference sample 6 calculate wilcoxon's upper and lower probabilities
#' wprob( 4, 6 )
#' 
#' @export
#' @docType methods
#' @rdname wprob-methods
"wprob" <-
function( m, n ) {
	# m - TEST sample size
	# n - REFERENCE sample size
	
	# unique m+n by m combinations
	cmb <- combn(m+n,m)
	
	# calculate rank sums
	rsum <- apply(cmb,2,sum)
	
	# get the frequency of each rank sum
	frsum <- table(rsum) / sum(table(rsum))
	
	# calculate lower-tail cumulative probabilities
	lowert <- cumsum(frsum)
	
	# calculate lower-tail cumulative probabilities
	uppert <- cumsum(rev(frsum))
	
	list(lower=lowert,upper=uppert)
}


#' Cluster probes that have no statistically significant differential signalling
#' 
#' Performs \code{wilcox.test} on normalized expression value matrix defined in \code{DEMIClust}
#' object and selects only these probes that have no differential signalling.
#' 
#' @param x A \code{DEMIClust} object. The DEMIClust object containing normalized expression
#' 		  values used for statistical significance test on differential signalling
#' 		  of probes. The object contains the column indexes of groups (e.g. 'test'
#' 		  and 'control') used in the analysis.
#' @return A \code{list}. Returns a \code{list} containing probes that did not have statistically
#' 		   significant differential signalling.
#' @seealso \code{\link{wilcox.test}} which this function wraps.
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
#' # Now we can continue the example of the function demiequal
#' 
#' # Basic experiment set up
#' demiexp <- DEMIExperiment(analysis = 'gene', celpath = destfolder,
#' 			experiment = 'myexperiment', organism = 'homo_sapiens')
#' 
#' # Create clusters with default behaviour
#' demiclust <- DEMIClust( demiexp, group = c( "BRAIN", "UHR" ) )
#' 
#' # Retrieve probes whose differential signalling was not statistically significant
#' nosigprobes <- demiequal( demiclust )
#' 
#' # However it makes more sense to incorporate the method straight into \code{DEMIClust} object
#' demiclust <- DEMIClust( demiexp, group = c( "BRAIN", "UHR" ), clust.method = demiequal )
#' 
#' # Retrieve the probes whose differential signalling was not statistically significant
#' nosigprobes <- getCluster( demiclust )
#' 
#' # Since the function only produces one cluster with a sign '[E]' derming equal
#' head( nosigprobes[[grep("\\[E\\]", names( nosigprobes ))]] )
#' 
#' }
#' 
#' @export 
#' @docType methods
#' @rdname demiequal-methods
"demiequal" <-
function( x = "DEMIClust" )
{
	
	#cat( "*Clustering probes into 'equal' cluster based on differential probe signal using Wilcox rank sum test if the groups are bigger then 2\n" );
	cat( DEMIMessages$demiequal$main );
	
	ndata <- getNormMatrix( getExperiment( x ) );
	rows <- nrow( ndata );
	
	group1 <- getGroup( x )@indexA;
	group2 <- getGroup( x )@indexB;
	groupA <- getGroup( x )@groupA;
	groupB <- getGroup( x )@groupB;
	
	H <- numeric( rows );
	L <- numeric( rows );
	
	options( warn = -1 ) # Turn warnings off
	
	#	calculate the probabilities of the wilcoxon's rank sum test if sample size is below 15
	wprobs <- NULL;
	if ( length( group1 ) < 15 && length( group2 ) < 15 ) {
		#cat( "\tUsing custom approach to do Wilcoxon rank sum test\n" );
		cat( DEMIMessages$demiequal$custom.approach );
		wprobs <- wprob( length( group1 ), length( group2 ) );
	} else if ( length( group1 ) >= 14 || length( group2 ) >= 14 ) {
		#cat( "\tUsing the standard Wilcoxon rank sum test function 'wilcox.test'\n" );
		cat( DEMIMessages$demiequal$standard.approach );
	}
	
	cat( DEMIMessages$takestime( rows ) );
	
	if ( length( group1 ) > 2 && length( group2 ) > 2 ) {
		# VERSION 1
		#	this function runs a lot quicker on matrix
		#	it uses a for cycle
		for ( i in 1:rows ) {
			if ( i %% 100000 == 0) {
				#cat( paste( "\t\t", i, " probes done\n", sep = "" ) );
				cat( DEMIMessages$probesdone( i ) );
			}
			
			#	if we use our custom wilcoxon's rank sum test
			if ( is.null( wprobs ) == FALSE ) {
				#	first check if the values are not duplicated
				if ( TRUE %in% duplicated( ndata[i,] ) ) {
					#	if the values are duplicated use the standard wilcox test
					L[i] <- wilcox.test( ndata[i, group1], ndata[i, group2], alternative = "less" )$p.value;
					H[i] <- wilcox.test( ndata[i, group1], ndata[i, group2], alternative = "greater" )$p.value;
				} else {
					#	if the values are not duplicated use the custom wilcox test
					#	rank the probe values over all samples
					pranks <- rank( ndata[i, c( group1, group2 )] );
					#	calculate the rank sum of the group A
					rsum <- sum( pranks[ seq( 1, length( group1 ), 1 ) ] );
					#	calculate the lower tail probability
					L[i] <- wprobs$lower[ which( names( wprobs$lower ) == rsum ) ]
					#	calculate the upper tail probability
					H[i] <- wprobs$upper[ which( names( wprobs$upper ) == rsum ) ]
				}
			} else {
				L[i] <- wilcox.test( ndata[i, group1], ndata[i, group2], alternative = "less", exact = FALSE )$p.value;
				H[i] <- wilcox.test( ndata[i, group1], ndata[i, group2], alternative = "greater", exact = FALSE )$p.value;
			}
		}
	} else {
		for ( i in 1:nrow( ndata ) ) {
			if ( i %% 100000 == 0) {
				cat( paste ( "\t", i, " probes done\n", sep = "" ) );
			}
			if ( max( ndata[i, group1] ) < min( ndata[i, group2] ) ) {
				L[i] <- 0;
				H[i] <- 1;
			} else if ( min( ndata[i, group1] ) > max( ndata[i, group2] ) ) {
				H[i] <- 0;
				L[i] <- 1;
			} else {
				L[i] <- 1;
				H[i] <- 1;
			}
		}
	}
	
	R <- list( equal = as.vector( rownames( ndata )[ intersect( which( L > getCutoffPvalue( x ) ), which( H > getCutoffPvalue( x ) ) ) ] ) );
	names( R ) <- paste( groupA,"[E]_", groupB, sep = "" );

	options( warn = 1 ); # Turn warnings back on
	
	return( R );
	
}#demiequal


#' Cluster probes into higher and lower clusters based on their differential signalling
#' 
#' Performs a modified \code{wilcox.test} on normalized expression value matrix defined in \code{DEMIClust}
#' object. It precalculates the probabilities of the rank sums and makes the algorithm run a lot quicker.
#' 
#' @param x A \code{DEMIClust} object. The \code{DEMIClust} object containing normalized expression
#' 		  values used for statistical significance test on differential signalling
#' 		  of probes. The object contains the column indexes of groups (e.g. 'test'
#' 		  and 'control') used in the analysis.
#' @return A \code{list}. Returns a \code{list} containing different sets of probes that behave
#' 		   similarly under current statistical test (e.g. up- or down-regulated probes).
#' @seealso \code{\link{wilcox.test}} which this function mimics and \code{\link{wprob}} which this function
#' 		  	implements.
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
#' # Now we can continue the example of the function demi.wilcox.test.fast
#' 
#' # Basic experiment set up.
#' demiexp <- DEMIExperiment(analysis = 'gene', celpath = destfolder,
#' 			experiment = 'myexperiment', organism = 'homo_sapiens')
#' 
#' # Create clusters with default behaviour
#' demiclust <- DEMIClust( demiexp, group = c( "BRAIN", "UHR" ) )
#' 
#' # Retrieve probes whose differential signalling was statistically significant
#' sigprobes <- demi.wilcox.test.fast( demiclust )
#' 
#' # However it makes more sense to incorporate the method straight into \code{DEMIClust} object
#' demiclust <- DEMIClust( demiexp, group = c( "BRAIN", "UHR" ), clust.method = demi.wilcox.test.fast )
#' 
#' # Retrieve the probes whose differential signalling was statistically significant
#' sigprobes <- getCluster( demiclust )
#' 
#' # Retrieve the cluster names since we have both up-regulated and down-regulated probe clusters
#' names( sigprobes )
#' 
#' # Retrieve the up-regulated probes whose cluster names contain the sign '[H]'
#' head( sigprobes[[grep("\\[H\\]", names( sigprobes ))]] )
#' 
#' # Retrieve the down-regulated probes whose cluster names contain the sign '[L]'
#' head( sigprobes[[grep("\\[L\\]", names( sigprobes ))]] )
#' 
#' }
#' 
#' @export 
#' @docType methods
#' @rdname demi.wilcox.test.fast-methods
"demi.wilcox.test.fast" <-
function( x = "DEMIClust" )
{
	
	#cat( "*Clustering probes into 'higher' and 'lower' clusters based on differential probe signal using Wilcox rank sum test\n" );
	cat( DEMIMessages$demi.wilcox.test.fast$main );
	
	ndata <- getNormMatrix( getExperiment( x ) );
	rows <- nrow( ndata );
	
	group1 <- getGroup( x )@indexA;
	group2 <- getGroup( x )@indexB;
	groupA <- getGroup( x )@groupA;
	groupB <- getGroup( x )@groupB;
	
	H <- numeric( rows );
	L <- numeric( rows );
	
	options( warn = -1 ) # Turn warnings off
	
	#	calculate the probabilities of the wilcoxon's rank sum test if sample size is below 15
	wprobs <- NULL;
	if ( length( group1 ) < 15 && length( group2 ) < 15 ) {
		#cat( "\tUsing custom approach to do Wilcoxon rank sum test\n" );
		cat( DEMIMessages$demi.wilcox.test.fast$custom.aproach );
		wprobs <- wprob( length( group1 ), length( group2 ) );
	} else if ( length( group1 ) >= 14 || length( group2 ) >= 14 ) {
		#cat( "\tUsing the standard Wilcoxon rank sum test function 'wilcox.test'\n" );
		cat( DEMIMessages$demi.wilcox.test.fast$standard.approach );
	}
	
	#cat( paste( "\tIt can take some time for there are a total of", rows, "probes to be analyzed.\n" ) );
	cat( DEMIMessages$takestime( rows ) );
	
	# VERSION 1
	#	this function runs a lot quicker on matrix
	#	it uses a for cycle
	for ( i in 1:rows ) {
		if ( i %% 100000 == 0) {
			#cat( paste( "\t\t", i, " probes done\n", sep = "" ) );
			cat( DEMIMessages$probesdone( i ) );
		}
		
		#	if we use our custom wilcoxon's rank sum test
		if ( is.null( wprobs ) == FALSE ) {
			#	first check if the values are not duplicated
			if ( TRUE %in% duplicated( ndata[i,] ) ) {
				#	if the values are duplicated use the standard wilcox test
				L[i] <- wilcox.test( ndata[i, group1], ndata[i, group2], alternative = "less" )$p.value;
				H[i] <- wilcox.test( ndata[i, group1], ndata[i, group2], alternative = "greater" )$p.value;
			} else {
				#	if the values are not duplicated use the custom wilcox test
				#	rank the probe values over all samples
				pranks <- rank( ndata[i, c( group1, group2 )] );
				#	calculate the rank sum of the group A
				rsum <- sum( pranks[ seq( 1, length( group1 ), 1 ) ] );
				#	calculate the lower tail probability
				L[i] <- wprobs$lower[ which( names( wprobs$lower ) == rsum ) ]
				#	calculate the upper tail probability
				H[i] <- wprobs$upper[ which( names( wprobs$upper ) == rsum ) ]
			}
		} else {
			L[i] <- wilcox.test( ndata[i, group1], ndata[i, group2], alternative = "less", exact = FALSE )$p.value;
			H[i] <- wilcox.test( ndata[i, group1], ndata[i, group2], alternative = "greater", exact = FALSE )$p.value;
		}
	}
	
	# VERSION 2 - uncompleted.
	#	use apply and pre-calculated rank table
#	pranks <- apply( ndata, 1, rank )
#	
#	if ( is.null( wprobs ) == FALSE ) {
#		pvals <- apply( pranks, 2, function(x) {
#					if ( TRUE %in% duplicated( x ) ) {
#						lower <- wilcox.test( x[group1], x[group2], alternative = "less" )$p.value;
#						upper <- wilcox.test( x[group1], x[group2], alternative = "greater" )$p.value;
#						c( lower, upper )
#					} else {
#						rsum <- sum( x[group1] )
#						lower <- 1
#						upper <- 1
#						lower <- wprobs$lower[ which( names( wprobs$lower ) == rsum ) ]
#						upper <- wprobs$upper[ which( names( wprobs$upper ) == rsum ) ]
#						c( lower, upper )
#					}
#				})
#	} else {
#		pvals <- apply( pranks, 2, function(x) {
#					lower <- wilcox.test( x[group1], x[group2], alternative = "less" )$p.value;
#					upper <- wilcox.test( x[group1], x[group2], alternative = "greater" )$p.value;
#					c( lower, upper )
#				})
#	}
#	
#	L <- as.numeric( pvals[1,] );
#	H <- as.numeric( pvals[2,] );
	
	# CONTINUES FOR BOTH CONDITIONS
	#	results object - contains both H and L
	R <- list( up_in_A_VS_B = H, dn_in_A_VS_B = L );
	names(R$up_in_A_VS_B) <- rownames( ndata );
	names(R$dn_in_A_VS_B) <- rownames( ndata );
	
	R$up_in_A_VS_B <- as.vector( names( which( R$up_in_A_VS_B <= getCutoffPvalue( x ) ) ) );
	R$dn_in_A_VS_B <- as.vector( names( which( R$dn_in_A_VS_B <= getCutoffPvalue( x ) ) ) );
	
	names( R )[ grep( "up_in_A_VS_B", names( R ) )] <- paste( groupA,"[H]_", groupB, sep = "" );
	names( R )[ grep( "dn_in_A_VS_B", names( R ) )] <- paste( groupA,"[L]_", groupB, sep = "" );
	
	options( warn = 1 ); # Turn warnings back on
	
	return( R );
	
}#demi.wilcox.test.fast


#' Cluster probes into higher and lower clusters based on their differential signalling
#' 
#' Performs \code{wilcox.test} on normalized expression value matrix defined in \code{DEMIClust}
#' object.
#' 
#' @param x A \code{DEMIClust} object. The \code{DEMIClust} object containing normalized expression
#' 		  values used for statistical significance test on differential signalling
#' 		  of probes. The object contains the column indexes of groups (e.g. 'test'
#' 		  and 'control') used in the analysis.
#' @return A \code{list}. Returns a \code{list} containing different sets of probes that behave
#' 		   similarly under current statistical test (e.g. up- or down-regulated probes).
#' @seealso \code{\link{wilcox.test}} which this function wraps.
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
#' # Now we can continue the example of the function demi.wilcox.test
#' 
#' # Basic experiment set up
#' demiexp <- DEMIExperiment(analysis = 'gene', celpath = destfolder,
#' 			experiment = 'myexperiment', organism = 'homo_sapiens')
#' 
#' # Create clusters with default behaviour
#' demiclust <- DEMIClust( demiexp, group = c( "BRAIN", "UHR" ) )
#' 
#' # Retrieve probes whose differential signalling was statistically significant
#' sigprobes <- demi.wilcox.test( demiclust )
#' 
#' # However it makes more sense to incorporate the method straight into \code{DEMIClust} object
#' demiclust <- DEMIClust( demiexp, group = c( "BRAIN", "UHR" ), clust.method = demi.wilcox.test )
#' 
#' # Retrieve the probes whose differential signalling was statistically significant
#' sigprobes <- getCluster( demiclust )
#' 
#' # Retrieve the cluster names since we have both up-regulated and down-regulated probe clusters
#' names( sigprobes )
#' 
#' # Retrieve the up-regulated probes whose cluster names contain the sign '[H]'
#' head( sigprobes[[grep("\\[H\\]", names( sigprobes ))]] )
#' 
#' # Retrieve the down-regulated probes whose cluster names contain the sign '[L]'
#' head( sigprobes[[grep("\\[L\\]", names( sigprobes ))]] )
#' 
#' }
#' 
#' @export 
#' @docType methods
#' @rdname demi.wilcox.test-methods
"demi.wilcox.test" <-
function( x = "DEMIClust" )
{
	
	#cat( "*Clustering probes into 'higher' and 'lower' clusters based on differential probe signal using Wilcox rank sum test\n" );
	cat( DEMIMessages$demi.wilcox.test$main );
	
	ndata <- getNormMatrix( getExperiment( x ) );
	rows <- nrow( ndata );
	
	group1 <- getGroup( x )@indexA;
	group2 <- getGroup( x )@indexB;
	groupA <- getGroup( x )@groupA;
	groupB <- getGroup( x )@groupB;
	
	H <- numeric( rows );
	L <- numeric( rows );
	
	options( warn = -1 ) # Turn warnings off
	
	#cat( paste( "\tIt can take some time for there are a total of", rows, "probes to be analyzed.\n" ) );
	cat( DEMIMessages$takestime( rows ) );
	
	#	this function runs a lot quicker on matrix
	#	it uses a for cycle
	if ( length( group1 ) > 15 && length( group2 ) > 15 ) {
		#cat( paste( "\tUsing wilcox with parameter 'exact = FALSE' for there are many samples in both groups.\n" ) );
		cat( DEMIMessages$demi.wilcox.test$exact.false );
		for ( i in 1:rows ) {
			if ( i %% 100000 == 0) {
				#cat( paste( "\t\t", i, " probes done\n", sep = "" ) );
				cat( DEMIMessages$probesdone( i ) );
			}
			L[i] <- wilcox.test( ndata[i, group1], ndata[i, group2], alternative = "less", exact = FALSE )$p.value;
			H[i] <- wilcox.test( ndata[i, group1], ndata[i, group2], alternative = "greater", exact = FALSE )$p.value;
		}
	} else {
		#cat( paste( "\tUsing wilcox with parameter 'exact = TRUE'.\n" ) );
		cat( DEMIMessages$demi.wilcox.test$exact.true );
		for ( i in 1:rows ) {
			if ( i %% 100000 == 0) {
				#cat( paste( "\t\t", i, " probes done\n", sep = "" ) );
				cat( DEMIMessages$probesdone( i ) );
			}
			L[i] <- wilcox.test( ndata[i, group1], ndata[i, group2], alternative = "less" )$p.value;
			H[i] <- wilcox.test( ndata[i, group1], ndata[i, group2], alternative = "greater" )$p.value;
		}
	}
	
	# CONTINUES FOR BOTH CONDITIONS
	#	results object - contains both H and L
	R <- list( up_in_A_VS_B = H, dn_in_A_VS_B = L );
	names(R$up_in_A_VS_B) <- rownames( ndata );
	names(R$dn_in_A_VS_B) <- rownames( ndata );
	
	R$up_in_A_VS_B <- as.vector( names( which( R$up_in_A_VS_B <= getCutoffPvalue( x ) ) ) );
	R$dn_in_A_VS_B <- as.vector( names( which( R$dn_in_A_VS_B <= getCutoffPvalue( x ) ) ) );
	
	names( R )[ grep( "up_in_A_VS_B", names( R ) )] <- paste( groupA,"[H]_", groupB, sep = "" );
	names( R )[ grep( "dn_in_A_VS_B", names( R ) )] <- paste( groupA,"[L]_", groupB, sep = "" );
	
	options( warn = 1 ); # Turn warnings back on
	
	return( R );
	
}#demi.wilcox.test


#' Cluster probes into higher and lower clusters based on their differential signalling
#' 
#' Performs t.test on normalized expression value matrix defined in 'DEMIClust'
#' object.
#' 
#' @param x A \code{DEMIClust} object. The \code{DEMIClust} object containing normalized expression
#' 		  values used for statistical significance test on differential signalling
#' 		  of probes. The object contains the column indexes of groups (e.g. 'test'
#' 		  and 'control') used in the analysis.
#' @return A \code{list}. Returns a \code{list} containing different sets of probes that behave
#' 		   similarly under current statistical test (e.g. up- or down-regulated probes).
#' @seealso \code{\link{t.test}} which this function wraps.
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
#' # Now we can continue the example of the function demi.t.test
#' 
#' # Basic experiment set up.
#' demiexp <- DEMIExperiment(analysis = 'gene', celpath = destfolder,
#' 			experiment = 'myexperiment', organism = 'homo_sapiens')
#'  
#' # Create clusters with default behaviour
#' demiclust <- DEMIClust( demiexp, group = c( "BRAIN", "UHR" ) )
#' 
#' # Retrieve probes whose differential signalling was statistically significant
#' sigprobes <- demi.t.test( demiclust )
#' 
#' # However it makes more sense to incorporate the method straight into \code{DEMIClust} object
#' demiclust <- DEMIClust( demiexp, group = c( "BRAIN", "UHR" ), clust.method = demi.t.test )
#' 
#' # Retrieve the probes whose differential signalling was statistically significant
#' sigprobes <- getCluster( demiclust )
#' 
#' # Retrieve the cluster names since we have both up-regulated and down-regulated probe clusters
#' names( sigprobes )
#' 
#' # Retrieve the up-regulated probes whose cluster names contain the sign '[H]'
#' head( sigprobes[[grep("\\[H\\]", names( sigprobes ))]] )
#' 
#' # Retrieve the down-regulated probes whose cluster names contain the sign '[L]'
#' head( sigprobes[[grep("\\[L\\]", names( sigprobes ))]] )
#' 
#' }
#' 
#' @export
#' @docType methods
#' @rdname demi.t.test-methods 
"demi.t.test" <-
function( x = "DEMIClust" )
{
	
	#cat( "*Clustering probes into 'higher' and 'lower' clusters based on differential probe signal using function 't.test'\n" );
	cat( DEMIMessages$demiExperiment_t.test$main );
	
	ndata <- getNormMatrix( getExperiment( x ) );
	rows <- nrow( ndata );
	
	#cat( paste( "\tIt can take some time for there are a total of", rows, "probes to be analyzed.\n" ) );
	cat( DEMIMessages$takestime( rows ) );
	
	group1 <- getGroup( x )@indexA;
	group2 <- getGroup( x )@indexB;
	groupA <- getGroup( x )@groupA;
	groupB <- getGroup( x )@groupB;
	
	H <- numeric( rows );
	L <- numeric( rows );
	
	for (i in 1:rows ) {
		#calculate t-test p-value
		if ( i %% 100000 == 0) {
			#cat( paste( "\t\t", i, " probes done\n", sep = "" ) );
			cat( DEMIMessages$probesdone( i ) );
		}
		
		pval <- t.test( ndata[i, group1], ndata[i, group2])$p.value;
		if ( mean( ndata[i, group1]) > mean( ndata[i, group2] ) ) {
			H[i] <- pval;
			L[i] <- 1;
		} else if ( mean( ndata[i, group1] ) < mean( ndata[i, group2] ) ) {
			H[i] <- 1;
			L[i] <- pval;
		} else if ( mean( ndata[i, group1] ) == mean( ndata[i, group2] ) ) {
			H[i] <- 1;
			L[i] <- 1;
		}
		
	}
	
	#	results object - contains both H and L
	R <- list( up_in_A_VS_B = H, dn_in_A_VS_B = L );
	names(R$up_in_A_VS_B) <- rownames( ndata );
	names(R$dn_in_A_VS_B) <- rownames( ndata );
	
	R$up_in_A_VS_B <- as.vector( names( which( R$up_in_A_VS_B <= getCutoffPvalue( x ) ) ) );
	R$dn_in_A_VS_B <- as.vector( names( which( R$dn_in_A_VS_B <= getCutoffPvalue( x ) ) ) );
	
	names( R )[ grep( "up_in_A_VS_B", names( R ) )] <- paste( groupA,"[H]_", groupB, sep = "" );
	names( R )[ grep( "dn_in_A_VS_B", names( R ) )] <- paste( groupA,"[L]_", groupB, sep = "" );
	
	options( warn = 1 ); # Turn warnings back on
	
	return( R );
	
}#demi.t.test

#' Cluster probes into higher and lower clusters based on their differential signalling
#' 
#' Performs higher or lower comparison test on normalized expression matrix defined
#' in the \code{DEMIClust} object. Only probes whose expression values in one group
#' are all either bigger or smaller then the expression values in the comparative group
#' are termed with significant differential expression. 
#' 
#' @param x A \code{DEMIClust} object. The \code{DEMIClust} object containing normalized expression
#' 		  values used for statistical significance test on differential signalling
#' 		  of probes. The object contains the column indexes of groups (e.g. 'test'
#' 		  and 'control') used in the analysis.
#' @return A \code{list}. Returns a \code{list} containing different sets of probes that behave
#' 		   similarly under current statistical test (e.g. up- or down-regulated probes).
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
#' # Now we can continue the example of the function demi.comp.test
#' 
#' # Basic experiment set up.
#' demiexp <- DEMIExperiment(analysis = 'gene', celpath = destfolder,
#' 			experiment = 'myexperiment', organism = 'homo_sapiens')
#' #' # Create clusters with default behaviour
#' demiclust <- DEMIClust( demiexp, group = c( "BRAIN", "UHR" ) )
#' 
#' # Retrieve probes whose differential signalling was statistically significant
#' sigprobes <- demi.comp.test( demiclust )
#' 
#' # However it makes more sense to incorporate the method straight into \code{DEMIClust} object
#' demiclust <- DEMIClust( demiexp, group = c( "BRAIN", "UHR" ), clust.method = demi.comp.test )
#' 
#' # Retrieve the probes whose differential signalling was statistically significant
#' sigprobes <- getCluster( demiclust )
#' 
#' # Retrieve the cluster names since we have both up-regulated and down-regulated probe clusters
#' names( sigprobes )
#' 
#' # Retrieve the up-regulated probes whose cluster names contain the sign '[H]'
#' head( sigprobes[[grep("\\[H\\]", names( sigprobes ))]] )
#' 
#' # Retrieve the down-regulated probes whose cluster names contain the sign '[L]'
#' head( sigprobes[[grep("\\[L\\]", names( sigprobes ))]] )
#' 
#' }
#' 
#' @export
#' @docType methods
#' @rdname demi.comp.test-methods 
"demi.comp.test" <-
function( x = "DEMIClust" )
{
	
	#cat( "*Clustering probes into 'higher' and 'lower' clusters\n" );
	cat( DEMIMessages$demi.comp.test$main );
	
	ndata <- getNormMatrix( getExperiment( x ) );
	rows <- nrow( ndata );
	
	#cat( paste( "\tIt can take some time for there are a total of", rows, "probes to be analyzed.\n" ) );
	cat( DEMIMessages$takestime( rows ) );
	
	group1 <- getGroup( x )@indexA;
	group2 <- getGroup( x )@indexB;
	groupA <- getGroup( x )@groupA;
	groupB <- getGroup( x )@groupB;
	
	H <- numeric( nrow( ndata ) );
	L <- numeric( nrow( ndata ) );
	
	for ( i in 1:nrow( ndata ) ) {
		if ( i %% 100000 == 0) {
			#cat( paste ( "\t", i, " probes done\n", sep = "" ) );
			cat( DEMIMessages$probesdone( i ) );
		}
		if ( max( ndata[i, group1] ) < min( ndata[i, group2] ) ) {
			L[i] <- 0;
			H[i] <- 1;
		} else if ( min( ndata[i, group1] ) > max( ndata[i, group2] ) ) {
			H[i] <- 0;
			L[i] <- 1;
		} else {
			L[i] <- 1;
			H[i] <- 1;
		}
	}
	
	
	#	results object - contains both H and L
	R <- list( up_in_A_VS_B = H, dn_in_A_VS_B = L );
	names(R$up_in_A_VS_B) <- rownames( ndata );
	names(R$dn_in_A_VS_B) <- rownames( ndata );
	
	R$up_in_A_VS_B <- as.vector( names( which( R$up_in_A_VS_B <= getCutoffPvalue( x ) ) ) );
	R$dn_in_A_VS_B <- as.vector( names( which( R$dn_in_A_VS_B <= getCutoffPvalue( x ) ) ) );
	
	names( R )[ grep( "up_in_A_VS_B", names( R ) )] <- paste( groupA,"[H]_", groupB, sep = "" );
	names( R )[ grep( "dn_in_A_VS_B", names( R ) )] <- paste( groupA,"[L]_", groupB, sep = "" );
	
	return( R );
	
}#demi.comp.test
