# Script comments and history
# 2014
# 5:35:25 PM

# Author: Federico Comoglio @ D-BSSE, ETH Zurich
###############################################################################

loadProtein <- function (pdbID, join.gaps = FALSE, cutoff = 7, ...) {
#	pkg <- require(bio3d)
#	if(!pkg) stop('The package bio3d is missing. bio3d is now available on CRAN. Please see ?install.packages for more information.')
	if (missing(pdbID)) 
		stop("fileImport: argument 'filename' missing, with no default\n")
	if (!is.character(pdbID)) 
		stop("fileImport: argument 'filename' must be string\n")    
	pdb <- read.pdb(file = pdbID, ...)
	tab.selection <- as.matrix(table(pdb$atom[,"chain"], pdb$atom[,"elety"]))
	tab.selection <- as.matrix(tab.selection [,"CA"])
	colnames(tab.selection) <- "#aminoacids"
	chain <- unique(pdb$atom[, "chain"])
	
	if (length(chain) !=1 ) {
		cat("Loading chains: \n")
		print(tab.selection)
	}
	
	coordinates <- vector( 'list', length(chain) )
	chain.names <- chain
	
	for(i in 1 : length(chain)) {
		subset <- pdb$atom[pdb$atom[, "chain"] == chain[i], ]
		alphatrace <- subset[subset[, "elety"] == "CA", ]
		tmp.coord <- as.matrix( alphatrace[, c("x","y","z")] )
		#check for gaps
		if(!join.gaps) { #execute default
			tmp <- findGaps(tmp.coord, cutoff = cutoff)
			coordinates[[i]] <- tmp #is a list of lists if gaps are found
			n.split <- length(tmp)
			split.names <- if(n.split == 1)
								chain[i]
							else 
								paste(chain[i], 1 : length(tmp), sep = '')
			names(coordinates[[i]]) <- split.names		
		}
		else
			coordinates[[i]] <- tmp.coord
	}
	if(join.gaps) 
		names(coordinates) <- chain.names
	else 
		coordinates <- unlist(coordinates, recursive = FALSE)
	
	#last control, single points (C or N-ter usually, esp. after split) -> matrix
	single.aa <- which( sapply(coordinates, class) == 'numeric' )	
	if( length(single.aa) >= 1) {
		for(i in 1 : length(single.aa)) {
			coordinates[[ single.aa[i] ]] <- matrix(coordinates[[ single.aa[i] ]], ncol = 3)
		}	
	}
	
	return(coordinates)
}


