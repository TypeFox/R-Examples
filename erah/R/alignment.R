
# QualityFactorControl <- function(Experiment)
# {
	
	# for(Smp in 1:length(Experiment@Data@FactorList))
	# {
		# Peak.Heights <- as.numeric(as.vector(unlist(Experiment@Data@FactorList[[Smp]]$"Peak Height")))
		# delete.factors <- which(Peak.Heights<Experiment@Data@Parameters$min.peak.height)
		
		# if(length(delete.factors)!=0)
		# {
			# FactorList.new <- as.data.frame(lapply(Experiment@Data@FactorList[[Smp]], function(x) x[-delete.factors]))
			# if(any(colnames(FactorList.new)=="Peak.Height")) colnames(FactorList.new)[which(colnames(FactorList.new)=="Peak.Height")] <- "Peak Height"

		# Experiment@Data@FactorList[[Smp]] <- FactorList.new	
		# }
	# }
	# Experiment 	
# }



align.factors <- function(factors.list, min.spectra.cor, max.time.dist, max.mz, mz.range)
{				
	
	stopifnot(min.spectra.cor<1,min.spectra.cor>0)
	
	empty.samples <- which(lapply(factors.list,nrow)==0)
	factors.list.original <- NULL
	if(length(empty.samples)!=0)
	{
		factors.list.original <- factors.list
		factors.list <- factors.list[-empty.samples] 
	}
	if(length(factors.list)==1) stop("Only one sample has been processed. No alignment needed")
	
	N.samples <- length(factors.list)
	
	factors.assignment.matrix <- apply(as.matrix(1:length(factors.list)),1,function(x) {
		as.data.frame(matrix(c(rep(x,length(factors.list[[x]]$ID)), 1:length(factors.list[[x]]$ID)), ncol=2))
		}) 
	
	factors.assignment.matrix <- do.call(rbind, factors.assignment.matrix)
	colnames(factors.assignment.matrix) <- c("Sample","Element")
	
	
	retention.time.vector <- lapply(factors.list,function(x){as.numeric(as.vector(x[,"RT"]))})
	retention.time.vector <- as.vector(unlist(retention.time.vector))
	ret.iterator <- as.matrix(1:length(retention.time.vector))
	
	## New Method:
	
	order.vector <- order(retention.time.vector)
	retention.time.vector.o <- retention.time.vector[order.vector]
	
	res.vector <- retention.time.vector.o - c(retention.time.vector.o[-1],retention.time.vector.o[length(retention.time.vector.o)])
	group.flags <- which(abs(res.vector)>max.time.dist)
	
	#Fixed Bug #AL_290316:
	if(length(group.flags)==0)
    {
    	time.dist.clustlist <- list(order.vector[1:length(order.vector)])
    }else{
	    if (group.flags[1] != 1) group.flags <- c(0, group.flags)
	    if (group.flags[length(group.flags)] != length(order.vector))  group.flags <- c(group.flags, length(order.vector))
	    time.dist.clustlist <- sapply(1:(length(group.flags) - 1), function(i) order.vector[(group.flags[i] + 1):group.flags[(i + 1)]])
	}
	
	#if(group.flags[1]!=1) group.flags <- c(0,group.flags)
	#if(group.flags[length(group.flags)]!=length(order.vector)) group.flags <- c(group.flags,length(order.vector))
		
	#time.dist.clustlist	<- sapply(1:(length(group.flags)-1), function(i) order.vector[(group.flags[i]+1):group.flags[(i+1)]])

	id.vector <- as.vector(unlist(lapply(factors.list,function(x){as.vector(x[,"ID"])})))
	global.class.vector <- unlist(apply(as.matrix(1:length(factors.list)),1,function(x){rep(x,length(factors.list[[x]][,1]))}))

	#uu <<- time.dist.clustlist

	del.inds <- which(unlist(lapply(time.dist.clustlist, length))==1)
	if(length(del.inds)!=0) time.dist.clustlist <- time.dist.clustlist[-del.inds]

	###################
	#k <- 1
	pb <- txtProgressBar(min=0,max=length(time.dist.clustlist), width=50, style=3)
	global.aligned.factors <- list()
	for(k in 1:length(time.dist.clustlist))		
	{
		local.clust <- time.dist.clustlist[k][[1]]
		class.vector <- global.class.vector[local.clust]
			
		time.dist <- as.matrix(dist(retention.time.vector[local.clust]), upper=T, diag=F)
				
		local.spectra.matrix <- sapply(local.clust, function(x){
			convertMSPspectra(factors.list[[factors.assignment.matrix[x,1]]]$Spectra[[factors.assignment.matrix[x,2]]],max.mz)
		})
		
		local.spectra.matrix[-mz.range,] <- 0  
		cor.dist <- suppressWarnings(cor(local.spectra.matrix))
		cor.dist[cor.dist<0] <- 0
		cor.dist <- 1 -  cor.dist

		cor.dist[abs(cor.dist)>(1-min.spectra.cor)] <- NA  ##Eliminar upper limits
		
		time.dist[abs(time.dist)>max.time.dist] <- NA
		max.eu.dist <- sqrt((1-min.spectra.cor)^2+max.time.dist^2)

		eu.dist <- sqrt(time.dist^2+cor.dist^2)
		eu.dist[eu.dist==0] <- NA
		eu.dist.vector <- which(eu.dist<max.eu.dist, arr.ind=T)
	
		eu.dist.graph <- graph.data.frame(eu.dist.vector, directed = FALSE)
		eu.dist.clustlist <- split(unique(as.vector(eu.dist.vector)), clusters(eu.dist.graph)$membership)

		clustlist.unit.length <- which(as.vector(unlist(lapply(eu.dist.clustlist,length)))==1)
		if(length(clustlist.unit.length)!=0) eu.dist.clustlist <- eu.dist.clustlist[-clustlist.unit.length]

		aligned.factors <-  lapply(eu.dist.clustlist,function(clust){
			#clust <- eu.dist.clustlist[[2]]
			
			inter.distance.matrix <- eu.dist[clust,clust]
			
			forbidden.combinations <- which((class.vector[clust] %*% t(class.vector[clust]))== class.vector[clust]^2, arr.ind=F)
			inter.distance.matrix[forbidden.combinations] <- NA
			inter.distance.matrix[inter.distance.matrix==0] <- NA
			
			clusts <- comp.clusters(inter.distance.matrix, class.vector[clust])
			clusts <- lapply(clusts,function(x) {clust[as.vector(x$elements)]})
			clusts
		})
		
		aligned.factors <- unlist(aligned.factors, recursive = FALSE)
		aligned.factors <- lapply(aligned.factors, function(x) local.clust[x])

		global.aligned.factors <- c(global.aligned.factors,aligned.factors)

	setTxtProgressBar(pb, getTxtProgressBar(pb)+1)
	}
		
	aligned.factors <- global.aligned.factors
			
	if(!(any(unlist(lapply(factors.list,function(x) {is.null(x$AlignID)}))==FALSE)))
	{	
		factors.list <- lapply(factors.list, function(x){
			outp <- cbind(x,matrix(0,nrow=length(x$ID)))
			colnames(outp)[ncol(outp)] <- "AlignID"
			outp
			})
	}else{
		factors.list <- lapply(factors.list, function(x){
			x$AlignID <- rep(0,nrow=length(x$ID))
			x
			})
	}
	
	for(i in 1:length(aligned.factors))
	{
		x <- aligned.factors[[i]] #[[1]]	
		loc <- factors.assignment.matrix[x,]
		for(j in 1:length(loc[,1])) factors.list[[loc[j,1]]]$AlignID[loc[j,2]] <- i	
	}
	if(!is.null(factors.list.original))
	{
		factors.list.original[-empty.samples] <- factors.list
		factors.list <- factors.list.original
	}
	
	factors.list
}

create.factorlist.table <- function(object)
{

	empty.samples <- which(lapply(object@Data@FactorList,nrow)==0)
	if(length(empty.samples)!=0) object@Data@FactorList <- object@Data@FactorList[-empty.samples]
	
	factors.list <- object@Data@FactorList
	
	#Spectra table
		
	alignId <- lapply(factors.list,function(x){x$AlignID})
	N.groups <- unique(unlist(alignId))
	N.groups <- N.groups[-which(N.groups==0)]
	N.samples <- length(object@Data@FactorList)
	samples.name <- names(object@Data@FactorList)		
	max.mz <- max(object@Results@Parameters@Alignment$mz.range)

	spectra.list <- apply(as.matrix(N.groups),1,function(Ng){
		#cat(Ng, " \n")
		align.iterator <- as.vector(which(lapply(alignId,function(x){
			if(length(intersect(Ng,x))!=0) {
				return(T)
				}else{
					return(F)
				}
			})==T))
			group.spectra <- lapply(as.matrix(align.iterator), function(x){ 
				list(spectra=convertMSPspectra(factors.list[[x]][which(alignId[[x]]==Ng),"Spectra"],max.mz),time=factors.list[[x]][which(alignId[[x]]==Ng),"RT"]) 
				})
			if(length(align.iterator)!=0)
			{	
				common.spectra <- normalize(rowSums(do.call(cbind,lapply(group.spectra,function(x){x$spectra})))) 
				common.time <- mean(do.call(cbind,lapply(group.spectra,function(x){x$time})))
				return(list(spectra=common.spectra,time=common.time, nsamples=length(align.iterator)))
			}else{return(NULL)}
	})
	
	#list.delete <- unlist(lapply(spectra.list, is.null))
	#if(any(list.delete)) spectra.list <- spectra.list[-which(list.delete==T)]
	#N.groups <- length(spectra.list)	

	spectra.matrix <- do.call(cbind,lapply(spectra.list,function(x){x$spectra}))
	#time.vector <- unlist(lapply(spectra.list,function(x){x$time}))
	#foundin.vector <- unlist(lapply(spectra.list,function(x){x$nsamples}))	
	

	for(i in 1:N.samples) samples.name[i] <- strsplit(as.character(samples.name[i]), split="\\.")[[1]][1]
	
	align.matrix <- t(apply(as.matrix(N.groups),1,function(i){
		local.index <- as.vector(unlist(lapply(object@Data@FactorList,function(x) {
			outp <- as.character(x[which(x$AlignID==i),"Peak Height"])
			if(length(outp)==0) outp <- 0
			as.numeric(as.vector(outp))
			})))
	}))
	
	foundIn.vector <- apply(align.matrix,1,function(x) length(which(x!=0)))
	
	time.align.matrix <- t(apply(as.matrix(N.groups),1,function(i){
		local.index <- as.vector(unlist(lapply(object@Data@FactorList,function(x) {
			outp <- as.character(x[which(x$AlignID==i),"RT"])
			if(length(outp)==0) outp <- 0
			outp
			})))
	}))
	
	class(time.align.matrix) <- "numeric"
	rt.mean <- rowSums(time.align.matrix)/apply(time.align.matrix,1,function(x) length(which(x!=0)))
	
	align.List <- matrix(0,nrow=nrow(align.matrix),ncol=(5+N.samples))
	colnames(align.List) <- c("AlignID","Factor","tmean","FoundIn","Spectra",as.character(samples.name))

	align.List[,6:(5+N.samples)] <- align.matrix
	
	# for(i in 1:N.groups) 
	# {
	#	local.pos <- which(object@Results@Identification$AlignID==i)
	#	align.List[i,"AlignID"] <- object@Results@Identification[local.pos,"AlignID"]
	#	align.List[i,"Name"] <- as.character(object@Results@Identification[local.pos,"Name"])
	#  	align.List[i,"tmean"] <- object@Results@Identification[local.pos,"tmean"]
	# 	#alignList[i,c(1:3)] <- object@Results@Identification[which(object@Results@Identification$AlignID==i),c("AlignID","Name","tmean")]
	# }
	
	align.List[,"AlignID"] <- N.groups
	align.List[,"tmean"] <- round(rt.mean, digits=4)
	align.List[,"Factor"] <- apply(align.List,1,function(x) paste("Factor #", x["AlignID"] , sep="")) 
	align.List[,"FoundIn"] <- foundIn.vector
	align.List[,"Spectra"] <- apply(as.matrix(1:ncol(spectra.matrix)),1,function(j){
			spectra <- spectra.matrix[,j]
			spectra.index <- which(spectra!=0) 
			spectra.pos <- spectra.index
			spectra.int <- round(spectra[spectra.index]*1000)
			spectra.text <- paste(sweep(as.matrix(spectra.pos),1,as.matrix(spectra.int),"paste.sp"), collapse=" ")
			spectra.text
		})	
	
	align.List <- as.data.frame(align.List[order(as.vector(as.numeric(align.List[,"tmean"]))),], row.names=1:nrow(align.List)) 
	align.List[,"tmean"] <- as.numeric(as.vector(align.List[,"tmean"]))
	align.List[,"FoundIn"] <- as.numeric(as.vector(align.List[,"FoundIn"]))
	align.List[,"AlignID"] <- as.integer(as.vector(align.List[,"AlignID"]))
	align.List[,"Spectra"] <- as.character(align.List[,"Spectra"])
	align.List[,"Factor"] <- as.character(align.List[,"Factor"])

	for(i in 6:ncol(align.List)) align.List[,i] <- as.numeric(as.vector(align.List[,i]))
	
	align.List
}


comp.clusters <- function(hdist, classes)
{
	k <- nrow(hdist)
	n.class <- unique(classes)
	
	group.list <- list()
	it <- 1
	
	while(1)
	{
		min.list <- list()
		min.index <- 0
		min.previous <- NA
		for(j in 1:k)
		{
			minn <- rep(NA,length(n.class))
			w.minn <- minn
			for(i in n.class)
			{
				D <- hdist[j,which(classes==i)]
				if(!all(is.na(D)))
				{ 
					minn[i] <- min(D, na.rm=T)
					w.minn[i] <- which(classes==i)[which.min(D)]
				}
			}
			group.distance <- sqrt(sum(minn^2, na.rm=T))/length(which(is.na(minn)==F))
			if(is.na(group.distance)) group.distance <- 0
			min.list[[j]] <- list(dist = group.distance, elements = w.minn)
			if(which.min(c(group.distance,min.previous))==1 && group.distance!=0) {min.previous <- group.distance; min.index <- j}
			#cat(j,": ", sqrt(sum(minn^2, na.rm=T)), " with N: ", length(which(is.na(minn)==T)) ,"\n", sep="")
		}
		if(min.index==0) break
		g.elements <- as.vector(na.omit(min.list[[min.index]]$elements))
		g.elements <- c(g.elements,min.index)
		if(!any(is.null(unlist(group.list)))) 
		{
			remove.i <- which((g.elements %in% unique(unlist(group.list)))==T)
			if(length(remove.i)!=0) g.elements <- g.elements[-remove.i]
		}
		#if(length(g.elements)==0) break;
		group.list[[it]] <- list(elements=g.elements)
		it <- it + 1
		#hdist[g.elements,g.elements] <- NA
		hdist[,g.elements] <- NA
		hdist[g.elements,] <- NA

	}
	group.list	
}





