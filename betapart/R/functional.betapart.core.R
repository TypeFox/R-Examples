functional.betapart.core<-function(x, traits, multi=TRUE, warning.time=TRUE, return.details=FALSE) 
{

	 


    if (!is.matrix(x)) {
        x <- as.matrix(x)
    }
    if (!is.numeric(x)) 
        stop("The data in 'x' is not numeric.", call. = TRUE)
    
    xvals <- unique(as.vector(x))
    if (any(!is.element(xvals, c(0, 1)))) 
        stop("The 'x' table contains values other than 0 and 1: data should be presence/absence.", call. = TRUE)
            
    if (!is.numeric(traits)) 
        stop("The data in 'traits' is not numeric.", call. = TRUE)    
         
    if (any(is.na(traits)))
    	stop("NA are not allowed in 'traits'", call. = TRUE)
    	
    if(ncol(x)!=nrow(traits)) 
    	stop("Number of species in 'x' and 'traits' must be identical", call. = TRUE) 
     
    D<-ncol(traits) 
    Si<-apply(x,1,sum)
    if (any(Si<=D ))
    	stop(paste("'community ",row.names(x)[which(Si<=D)]," must contain at least ",D+1, " species",sep=""))
                  
   	N<-nrow(x)
    if(N<2) 
    	stop("Computing dissimilairty requires at least 2 communities", call. = TRUE) 
   	
   	# dataframe to return computation step done
   	nb.step<-2 ; if(multi==T) nb.step<-N
   	step.fbc<-as.data.frame(matrix("",nb.step,1, dimnames=list( c("           FRi", paste("intersection",2:nb.step,sep="_")), c("iteration")  ) ) )
   	step.fbc[,1]<-as.character(step.fbc[,1])
   	step.fbc[1,1]<-paste("0/",N,sep="")
   	for (k in 2:nb.step)
	   	step.fbc[k,1]<-paste("0/",choose(N,k),sep="")	
   	
   		
   	# computing vertices coordinates and convex hull volume of each polytope
    FRi<-rep(NA,N) ; names(FRi)<-row.names(x)
    coord_vert_i<-list()
    
    for (i in 1:N)
    {tr_i<-traits[which(x[i,]==1),] # traits values present
    
    vert0<-convhulln(tr_i,"Fx TO 'vert.txt'") # vertices indentity
	vert1<-scan("vert.txt",quiet=T)
	verti<-(vert1+1)[-1]
		
    coord_vert_i[[i]]<-tr_i[verti,] # vertices coordinates
    
    FRi[i]<-convhulln(tr_i[verti,],"FA")$vol # convex hull volume
	
	# writing flag
	step.fbc["           FRi",1]<-paste(i,"/",N,sep="") ;   step.fbc[,1]<-as.character(step.fbc[,1])
	write.table(step.fbc, file="step.fbc.txt", row.names=T, col.names=F, sep="\t")
	} # end of i
	
	sumFRi <- sum(FRi)


	# function to compute coordinates of vertices and volume of intersection between two convex polytopes
	intersect<-function(set1, set2) # set= points coordinates (Si*D)
	{
	
		# tranforming species coordinates to true rational number written as character strings
		set1rep <- d2q(cbind(0, cbind(1, set1)))
		set2rep <- d2q(cbind(0, cbind(1, set2)))

		# keeping only vertices
		polytope1 <- redundant(set1rep, representation = "V")$output
		polytope2 <- redundant(set2rep, representation = "V")$output

		# changing polytope representation: vertices coordinates -> inequality constraints
		H_chset1 <- scdd(polytope1, representation = "V")$output
		H_chset2 <- scdd(polytope2, representation = "V")$output

		# intersection between the two polytopes
		H_inter <- rbind(H_chset1, H_chset2)
		V_inter <- scdd(H_inter, representation = "H")$output
		vert_1n2 <- q2d(V_inter[ , - c(1, 2)])

    	# volume and vertices coordinates of intersection if it exists
    	coord_vert_inter<-rep(NA,ncol(set1))
    	vol_inter<-0
    	if (is.matrix(vert_1n2)==T) # vector if one vertex in common
      		if( nrow(vert_1n2)>ncol(vert_1n2) ) 
      			{ 
      				coord_vert_inter<-vert_1n2
      				vol_inter<-convhulln(vert_1n2,"FA")$vol
      				} # end of if intersection exists
		
		res<-list( coord_vert_inter=coord_vert_inter, vol_inter=vol_inter)
		return(res)
	} # end of intersect
	


	# pariwise comparisons	
	comb2<-combn(1:N, 2, simplify=T) # table with all pairs
	
	# matrix and vector to store volumes of intersection and list to store vertices coordinates
	vol_inter2_mat<-matrix(0,N,N,dimnames=list(row.names(x), row.names(x) ) )
	vol_inter2<-rep(0,ncol(comb2))
	coord_vert_inter2<-list()
	
	
	# loop on pairs
		for (k in 1:ncol(comb2))
		{
		# communities i and j for pair k
		i<-comb2[1,k]
		j<-comb2[2,k]
		
		# species coordinates in the multidimensional functional space 
		seti<-traits[which(x[i,]==1),]
		setj<-traits[which(x[j,]==1),]	
	
		# computing intersection
		interij<-intersect(seti, setj)
		vol_inter2_mat[j,i] <-interij$vol_inter
		vol_inter2[k]<-interij$vol_inter
		coord_vert_inter2[[k]]<-interij$coord_vert_inter
		
		# writing flag
		step.fbc["intersection_2",1]<-paste(k,"/",ncol(comb2),sep="")
		write.table(step.fbc, file="step.fbc.txt", row.names=T, col.names=F, sep="\t")
		} # end of k


    # functional richness shared or unique between pairs
	matNN<-matrix(0,N,N,dimnames=list(row.names(x), row.names(x) ) )
	shared<-matNN
	not.shared<-matNN
	
	for (i in 1:(N-1))
	for (j in (i+1):N)
	{
	shared [j,i] <- vol_inter2_mat[j,i]
	not.shared[i,j] <- FRi[i]-vol_inter2_mat[j,i]
	not.shared[j,i] <- FRi[j]-vol_inter2_mat[j,i]
	} # end of i,j
		
    # sum, max, min of functional richness not shared
	sum.not.shared <- not.shared + t(not.shared)
    max.not.shared <- pmax(not.shared, t(not.shared))
    min.not.shared <- pmin(not.shared, t(not.shared))
    
    
    # general lists to store  combinations, vertices coordinates and volumes of intersections
    # and filling them with results for pariwise intersections
    comb_inter<-list() ; comb_inter[[1]]<-comb2
    coord_vert_inter<-list() ; coord_vert_inter[[1]]<-coord_vert_inter2
    vol_inter<-list() ; vol_inter[[1]]<-vol_inter2
    
    # by default multiple functional a and total richness set to NA (i.e. meaningless for pairwise dissimilarity)
    FRt <- NA
    a <- NA
    
    # computing multiple dissimilarity only if it is required
    
    if( N>2 & multi==T  ) {
     	
     	if (warning.time==T & N>10 )
    		stop(paste("Computing mulitple functional dissimilarity on more than 10 communities may take a long time. 
    									Set 'multi' or 'warning.time' to FALSE")) 
    	
     	if (warning.time==T & D>4 )
    		stop(paste("Computing mulitple functional dissimilarity in a",D,"-dimensions functional space may take a long time. 
    									Set 'multi' or 'warning.time' to FALSE"))            


    	# loop for computing intersection from z=3 to z=N communities based on intersections at the z-1 level
    	for (z in 3:N)
    	{
		
		# table of combinations of z communities
    	comb_z<-combn(1:N, z, simplify=T) 
    	
		# matrix and vector to store volume of intersection and vertices coordinates
		vol_inter_z<-rep(0,ncol(comb_z))
		coord_vert_inter_z<-list()	

    
    		# loop on all z-combinations
			for (k in 1:ncol(comb_z))
			{
		
			# extracting vertices coordinates of the 2 intersections at z-1 level involved for the k-th intersection at z-level
			seti<-coord_vert_inter[[z-2]] [[ which( apply( comb_inter[[z-2]]  , 2, identical, comb_z[1:(z-1),k]  )==T) ]]
			setj<-coord_vert_inter[[z-2]] [[ which( apply( comb_inter[[z-2]]  , 2, identical, comb_z[2:z,k]  )==T) ]]	
	
			# k-th intersection at z level if both intersections at z-1 level exist, else NA
			coord_vert_inter_z[[k]]<-rep(NA,D)
				if( is.na( sum(seti)+sum(setj) )==F ) {
							interij<-intersect(seti, setj)
							vol_inter_z[k]<-interij$vol_inter
							coord_vert_inter_z[[k]]<-interij$coord_vert_inter } # end of intersections z-1 exists
			
			# writing flag	
			step.fbc[paste("intersection",z,sep="_"),1]<-paste(k,"/",ncol(comb_z),sep="")
			write.table(step.fbc, file="step.fbc.txt", row.names=T, col.names=F, sep="\t")

    		} # end of k
    	
    	# storing results
    	comb_inter[[z-1]]<-comb_z
    	coord_vert_inter[[z-1]]<-coord_vert_inter_z
    	vol_inter[[z-1]]<-vol_inter_z	
    
    	} # end of z
    
    
    	# inclusion-exclusion principle for computing volume of the union of the N communities (functional equivalent to St)
  		
  		# sum of volumes of intersections from 2 to N with sign (-1)^(k-1)
  		sumvol_sign<-rep(NA,N-1)
  			for (k in 2:N)
  				{   sumvol_sign[k-1]<- (-1)^(k-1) * sum(vol_inter[[k-1]])     } # end of k

			# volume of union
			FRt<- sumFRi +  sum(sumvol_sign)
  
  			# volume shared by at least two communities
    		a <- sumFRi - FRt
    
    } # end of if more than two communities
  
    
    # grouping details about convex hulls and their intersections in a list if they have to be returned
    details<-NA
   	if( return.details==T) 
   	   {  
   	   	CH<-list(FRi=FRi, coord_vertices=coord_vert_i)
   	   	intersections<-list(combinations=comb_inter, volumes=vol_inter, coord_vertices=coord_vert_inter)
   	   	details<-list(CH=CH, intersections=intersections)
   	   	} # end of if details	

   	
   	
  	# returning results of functional.betapart.core
    functional.computations<-list( sumFRi=sumFRi, FRt=FRt, a=a, 
    								shared = shared, not.shared = not.shared, sum.not.shared = sum.not.shared, 
    								max.not.shared = max.not.shared, min.not.shared = min.not.shared,
    								details=details)
    
    class(functional.computations) <- "functional.betapart"
    return(functional.computations)
	

} # end of function