##################################################################################
#' Sorting by NDS-rank and Hypervolume Contribution
#'
#' Sorts the large design for the purpose of multi objective optimization with SPOT.
#' First non dominated sorting rank (NDS) is used. If the choice of points for the next
#' sequential step is not clear by NDS rank, the hypervolume contribution of the 
#' competing points is recalculated sequentially to remove those with the smallest 
#' contribution.
#'
#' @param largeDesign the design matrix in the parameter space, to be sorted by the associated y-values for each objective
#' @param designY objective value matrix. Contains objective values associated to largeDesign
#' @param newsize this is the number of points that need to be selected, i.e. the seq.design.new.size
#' @return largeDesign \cr 
#' - The sorted large design
#' @export
#' @keywords internal
###################################################################################
spotMcoSort <- function (largeDesign, designY, newsize){
	lhdY=t(as.matrix(designY))
	ndsR<-nds_rank(lhdY) #First: Sort by nds rank
	largeDesign <-  as.data.frame(largeDesign[order(ndsR,decreasing=FALSE),]);
	lhdY<-lhdY[,order(ndsR,decreasing=FALSE)]
	ndsR<-ndsR[order(ndsR)]
	summe<-0
	if(newsize<nrow(largeDesign)){
		for(i in 1:max(ndsR)){ #Look for front that will not be used completely, sort it by hypervol contrib
			index<-which(ndsR==i)
			summe<-summe+length(index)
			if(summe==newsize)break;
			if((summe>newsize) && (length(index)>1)){
				set<-largeDesign[index,]
				removeN=summe-newsize
				sortVec=rep(0,length(index))
				frontY<-lhdY[,index]
				for(jj in 1:removeN){ #repeated selection by hypervolume contribution, selected individuals will be last in order
					iREM <- nds_hv_selection(frontY[,sortVec==0])
					iREM <- which(frontY==frontY[,sortVec==0][,iREM],arr.ind=TRUE)[1,2] #ugly hack to select which column to remove by comparing with original front
					sortVec[iREM]=1						
				}
				set<-set[order(sortVec),]
				largeDesign[index,]<-set
				break;
			}
		}
	}
	largeDesign
}


##################################################################################
#' Sorting by NDS-rank and Hypervolume Contribution, with known points
#'
#' Sorts the large design for the purpose of multi objective optimization with SPOT.
#' First non dominated sorting (NDS) rank is used. If the choice of points for the next
#' sequential step is not clear by NDS rank, the hypervolume contribution of the 
#' competing points is recalculated sequentially to remove those with the smallest 
#' contribution. 
#'
#' In contrast to \code{\link{spotMcoSort}}, this function considers the known points in \code{mergedX} and \code{mergedY}
#' so that new points will rather be chosen in between known points, thus producing a better Pareto front.
#' To do so, the known points are added to the set of solutions. To ensure that they are not removed, they receive infinite hypervolume contribution,
#' and are not counted when determining the number of NDS ranks to be considered.
#'
#' @param largeDesign the design matrix in the parameter space, to be sorted by the associated y-values for each objective
#' @param designY objective value matrix. Contains objective values associated to largeDesign
#' @param newsize this is the number of points that need to be selected, i.e. the seq.design.new.size
#' @param mergedX position of the already known points in parameter space (vector of parameter values)
#' @param mergedY y-values of the already known points (vector of objective values)
#' @return largeDesign \cr 
#' - The sorted large design
#' @export
#' @keywords internal
###################################################################################
spotMcoSelectionHypervol <- function (largeDesign, designY, newsize, mergedX, mergedY,ref=NULL){
	info<-c(rep(0,nrow(designY)),rep(1,nrow(mergedY))) # create info of point origin (0=largedesign, 1=allready evaluated)
	colnames(designY)<-colnames(mergedY)
	designY<-rbind(designY,mergedY)#merge points and known values	
	largeDesign<-rbind(largeDesign,mergedX)#merge large design and known values	
	lhdY=t(as.matrix(designY))
	ndsR<-nds_rank(lhdY) #determine nds rank
	info<-info[order(ndsR,decreasing=FALSE)] #sort info about points by nds rank
	largeDesign <-  as.data.frame(largeDesign[order(ndsR,decreasing=FALSE),]) #sort large design by nds rank
	lhdY<-lhdY[,order(ndsR,decreasing=FALSE)] #sort y values by nds rank
	ndsR<-ndsR[order(ndsR)]#sort ndsR values by nds rank
	summe<-0 #initial value for sum
	if(newsize<nrow(largeDesign)){
		for(i in 1:max(ndsR)){ #Look for front that will not be used completely, sort it by hypervol contrib
			index<-which(ndsR==i)
			summe<-summe+length(index)- sum(info[index]) #sum of points from this point (which will be added): all points, minus known points
			if(summe==newsize)break;
			if((summe>newsize) && (length(index)>1)){
				set<-largeDesign[index,]
				removeN=summe-newsize
				sortVec=rep(0,length(index))
				frontY<-lhdY[,index]
				frontInfo<-info[index]
				for(jj in 1:removeN){ #repeated selection by hypervolume contribution, selected individuals will be last in order
					if(is.null(ref)){
						contrib<-hypervolume_contribution(frontY[,sortVec==0]) #todo: check if this algorithm works robustly
					}else{
						contrib<-hypervolume_contribution(frontY[,sortVec==0],ref)
					}
					contrib[frontInfo==1]<-Inf#max(contrib)+1 #make sure known points can not be removed!    
					iREM=which.min(contrib)    
					iREM= which(frontY==frontY[,sortVec==0][,iREM],arr.ind=TRUE)[1,2] #ugly hack to select which column to remove by comparing with original front
					sortVec[iREM]=1						
				}
				set<-set[order(sortVec),]
				info[index]<-frontInfo[order(sortVec)]
				largeDesign[index,]<-set #fill sorted pareto set into large design
				break;
			}
		}
	}	
	largeDesign<-largeDesign[info==0,]#remove known points from largedesign
	largeDesign
}

##################################################################################
#' Multi criteria sorting by a tournament scheme
#'
#' Experimental, do not use
#'
#' @param largeDesign the design matrix in the parameter space, to be sorted by the associated y-values for each objective
#' @param designY objective value matrix. Contains objective values associated to largeDesign
#' @param tsize tournament size (percentage of large design in values between zero and one)
#' @param ssize selection size number of points that will be selected in SPOT, lower bound for absolute tsize
#' @return largeDesign \cr 
#' - The sorted large design
#' @keywords internal
###################################################################################
spotMcoCrowdTournament <- function (largeDesign, designY, tsize,ssize){			
	##### load rgp (needed for nondeterministicRanking orderByParetoCrowdingDistance)
	spotInstAndLoadPackages("rgp");
	#####calculate number of competitors, repair
	tsize <- tsize*nrow(largeDesign)
	tsize <- max(2,ssize,tsize)
	tsize <- min(nrow(largeDesign),tsize)
	winners<-NULL
	i<-0
	while(i<ssize){
		selection <- sample(nrow(largeDesign),tsize)
		competitors <- largeDesign[selection,]
		competitors <- competitors[rgp::orderByParetoCrowdingDistance(t(as.matrix(designY[selection,]))),]
		winners<-rbind(competitors[1,],winners)
		winners<-unique(winners)
		i<-nrow(winners)
	}
	winners
}

##################################################################################
#' Pareto Tournament
#'
#' Experimental, do not use
#'
#' @param largeDesign the design matrix in the parameter space, to be sorted by the associated y-values for each objective
#' @param designY objective value matrix. Contains objective values associated to largeDesign
#' @param tsize tournament size (percentage of large design in values between zero and one)
#' @param ssize selection size number of points that will be selected in SPOT, lower bound for absolute tsize
#' @return largeDesign \cr 
#' - The sorted large design
#' @keywords internal
###################################################################################
spotMcoTournament <- function (largeDesign, designY, tsize, ssize){
	spotInstAndLoadPackages("emoa");
	#####calculate number of competitors, repair
	tsize <- tsize*nrow(largeDesign)
	tsize <- max(2,ssize,tsize)
	tsize <- min(nrow(largeDesign),tsize)
	#####sample competitors
	i<-0
	winners<-NULL
	while(i<ssize){   #Scheme: all non dominated points win the tournament. Problem: usually those will be more than ssize...
						#that means: all that the tournament really does is randomly disregard solutions in one single step. and the
						#rest is not even sorted in any way
						
						#big problem: due to previous sms-emoa or nsga2 most solutions in large design ARE non dominated... 
						#not much to do for this type of tournament
			
						# smaller problem: new parameter introduced.... tuning necessary.
		selection <- sample(nrow(largeDesign),tsize)
		competitors <- largeDesign[selection,]
		competitorsY <- designY[selection,]
		winidx<-!is_dominated(t(competitorsY))
		winner<-competitors[winidx,]
		winners<-rbind(winner,winners)
		winners<-unique(winners)
		i<-nrow(winners)
	}
	winners[1:ssize,]
}