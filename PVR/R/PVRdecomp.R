PVRdecomp <- function(phy = NULL, type = "newick", dist = NULL, scale = FALSE, ...){
	##############---------The function and its arguments
	
	##############---------Loading and preparing data

#---------loading the phylogeny
	
	if(is.null(phy) & is.null(dist)){
		
		cat("Enter the phylogeny (in newick format). Remember that the last string must be blank", "\n")			
		phy <- read.tree()	
	} else
	
	if(!is.null(phy) & is.null(dist)){
		
		if(is.character(phy)){
			
			if(type == "newick"){
				phy = read.tree(phy)
			} else
			
			if(type == "nexus"){
				
				phy = read.nexus(phy)
			}
		} else {
			
			if(class(phy) == "phylo"){
				
				if (is.null(phy$edge.length)){ 
					stop("the phylogeny has no branch lengths")
				}
			} 
			
		}

	} else {
	
		if(!is.null(dist)){
			
			if(!is.null(phy)){
				
				warning("Both phylogeny and phylogenetic distance matrix have data.")
				if(is.character(dist)){
					
					pD <- read.table(dist)
				} else {
					
					if(is.matrix(dist)){
						
						pD <- dist
					}
				}
				
			} else {
			
				if(is.null(phy)){
					
					if(is.character(dist)){
						
						pD <- read.table(dist)
					} else {
						
						if(is.matrix(dist)){
							
							pD <- dist
						}
					}
				}
			}	
		}		
	}
	
#---------Creating the phylogenetic distances matrix	
	if(is.null(dist)){
		
		pD <- cophenetic.phylo(phy)
	}
	
#---------scaling phylogenetic distances into range 0 to 1
	if(scale){
		
		pD <- pD/max(pD)
	}
	
#---------Double centering pD (without squaring it)
	A <- as.matrix(-0.5*(pD))
	l <- matrix(1, nrow = nrow(pD), ncol = 1)
	L <- l%*%t(l)
	pD <- (diag(ncol(pD)) - ((1/ncol(pD)) * L))%*%A%*%(diag(ncol(pD)) - ((1/ncol(pD)) * L))
	
	
#	#pD <- -0.5*(pD)
#	di_hat <- matrix(rowSums(pD)/(nrow(pD)-1), ncol = ncol(pD), nrow = nrow(pD))
#	dj_hat <- matrix(colSums(pD)/(ncol(pD)-1), ncol = ncol(pD), nrow = nrow(pD))
#	pD <- pD-di_hat-dj_hat + (sum(pD)/((nrow(pD)*ncol(pD))-nrow(pD)))	
#	pD <- -0.5*(pD)
	
	##############---------Computing PVR
	####---------Computig pvr (spectral decomposition of matrix pD using LAPACK routines through base::eigen function)
	pvr <- (eigen(pD, symmetric = TRUE))
	Naxis <- ncol(pD) - sum((pvr$values <= 0)*1)
	pvr$values <- pvr$values[1:Naxis]
	pvr$vectors <- pvr$vectors[ ,1:Naxis]
	v <- "c1"
	for(i in 2:ncol(pvr$vectors)){
		
		v <- c(v, paste("c", i, sep = ""))
	}
	colnames(pvr$vectors) <- v
	names(pvr$values) <- v
	
	results <- new("PVR")
	
	results@Eigen <- pvr
	results@phyDist <- pD
	results@phylo <- phy
	
	return(results)
}