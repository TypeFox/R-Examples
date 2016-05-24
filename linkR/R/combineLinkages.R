combineLinkages <- function(linkage1, linkage2, linkage3=NULL, linkage4=NULL){

	if(!is.null(linkage3) && is.null(linkage4)){
		combine_linkages <- combineLinkages(linkage1, linkage2)
		return(combineLinkages(combine_linkages, linkage3))
	}

	if(!is.null(linkage3) && !is.null(linkage4)){
		combine_linkages <- combineLinkages(linkage1, linkage2)
		combine_linkages <- combineLinkages(combine_linkages, linkage3)
		return(combineLinkages(combine_linkages, linkage4))
	}

	# CHECK THAT LINK NAMES ARE DEFINED FOR LINKAGE INPUTS
	if(class(linkage1) == 'linkage' && is.null(linkage1$link.names)) stop("'link.names' is NULL for linkage1. 'link.names' must be defined for linkage inputs.")
	if(class(linkage2) == 'linkage' && is.null(linkage2$link.names)) stop("'link.names' is NULL for linkage2. 'link.names' must be defined for linkage inputs.")

	# CHECK THAT THERE IS OVERLAP IN LINK NAMES BETWEEN LINKAGE(S)/LINKAGE SYSTEM

	if(!is.null(linkage1$points) || !is.null(linkage2$points)){
		
		# START COMBINE POINT MATRICES
		points <- linkage1$points

		# COMBINE LINK ASSOCIATIONS
		link.assoc <- linkage1$link.assoc
		link.assoc <- c(link.assoc, linkage2$link.assoc[!rownames(linkage2$points) %in% rownames(points)])
		
		# COMBINE LINK NAMES
		link.names <- unique(c(linkage1$link.names, linkage2$link.names))

		# FINISH COMBINE POINT MATRICES
		points <- rbind(points, linkage2$points[!rownames(linkage2$points) %in% rownames(points), ])
		
		# RE-SET THE POINTS ASSOCIATED WITH EACH LINK
		points.assoc <- setNames(vector("list", length(link.names)), link.names)			
		
		# IF LINK.ASSOC ARE NUMERIC INTEGERS
		if(is.numeric(link.assoc[1])){
			for(i in 1:length(link.assoc))
				points.assoc[[names(points.assoc)[link.assoc[i]]]] <- c(points.assoc[[names(points.assoc)[link.assoc[i]]]], i)
		}else{
			for(i in 1:length(link.assoc)) points.assoc[[link.assoc[i]]] <- c(points.assoc[[link.assoc[i]]], i)
		}
		
	}

	# COMBINE LINKAGES
	if(class(linkage1) == 'linkage' && class(linkage2) == 'linkage'){

		linkage_system <- list()
		linkage_system[[1]] <- linkage1
		linkage_system[[2]] <- linkage2

	}else if(class(linkage1) == 'linkage_system' && class(linkage2) == 'linkage'){

		linkage_system <- list()		
		linkage_system[1:sum(names(linkage1) == "")] <- linkage1[1:sum(names(linkage1) == "")]
		linkage_system[[length(linkage_system)+1]] <- linkage2		

	}else if(class(linkage1) == 'linkage' && class(linkage2) == 'linkage_system'){

		linkage_system <- list()		
		linkage_system[1:sum(names(linkage2) == "")] <- linkage1[1:sum(names(linkage2) == "")]
		linkage_system[[length(linkage_system)+1]] <- linkage1		

	}else if(class(linkage1) == 'linkage_system' && class(linkage2) == 'linkage_system'){

		linkage_system <- linkage1
		linkage_ct1 <- sum(names(linkage_system) == "")
		linkage_ct2 <- sum(names(linkage2) == "")

		for(i in 1:linkage_ct2) linkage_system[[linkage_ct1+i]] <- linkage2[[i]]

	}

	linkage_system$link.names <- link.names
	linkage_system$points <- points
	linkage_system$link.assoc <- link.assoc
	linkage_system$points.assoc <- points.assoc

	#if(class(linkage1) == 'linkage_system') print(linkage_system[[3]][['joints']])

	class(linkage_system) <- 'linkage_system'

	linkage_system
}