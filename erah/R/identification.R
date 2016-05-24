cor.sinus <- function(s1.mat, s2.mat)
{
	if(is.matrix(s1.mat)==F) s1.mat <- matrix(s1.mat, ncol=1)
	if(is.matrix(s2.mat)==F) s2.mat <- matrix(s2.mat, ncol=1)
	
	Z <- matrix(diag(t(s1.mat)%*%s1.mat),ncol=1)
	H <- matrix(diag(t(s2.mat)%*%s2.mat),ncol=1)
	
	cos.angle <- (t(s1.mat)%*%s2.mat)/sqrt(Z%*%t(H))
	cos.angle
}

convertMSPspectra <- function(MSPchain, maxMZ)
{
	out.sp <- rep(0,maxMZ) 
	splitted.spectra.list <- strsplit(as.character(MSPchain),split=" ")[[1]]
	splitted.spectra.matrix <- apply(as.matrix(splitted.spectra.list),1,function(c.bin){strsplit(c.bin,split=",")[[1]]})
	
	if(all(is.na(splitted.spectra.matrix))) return(out.sp)	
	if(max(as.numeric(splitted.spectra.matrix[1,]))>maxMZ)
	{
		delete.mz <- which(as.numeric(splitted.spectra.matrix[1,])>maxMZ)
		splitted.spectra.matrix[,-delete.mz]
	}	
		
	out.sp[as.numeric(splitted.spectra.matrix[1,])] <- as.numeric(splitted.spectra.matrix[2,])
	if(length(out.sp)>maxMZ) out.sp <- out.sp[-((maxMZ+1):length(out.sp))] 
	out.sp
}

convertMSPspectra.dot <- function(MSPchain, maxMZ)
{
	out.sp <- rep(0,maxMZ) 
	splitted.spectra.list <- strsplit(as.character(MSPchain),split=" ")[[1]]
	splitted.spectra.matrix <- apply(as.matrix(splitted.spectra.list),1,function(c.bin){strsplit(c.bin,split=":")[[1]]})
		
	if(max(as.numeric(splitted.spectra.matrix[1,]))>maxMZ)
	{
		delete.mz <- which(as.numeric(splitted.spectra.matrix[1,])>maxMZ)
		splitted.spectra.matrix[,-delete.mz]
	}	
		
	out.sp[as.numeric(splitted.spectra.matrix[1,])] <- as.numeric(splitted.spectra.matrix[2,])
	if(length(out.sp)>maxMZ) out.sp <- out.sp[-((maxMZ+1):length(out.sp))] 
	out.sp
}

identify.factors <- function(Experiment, maxMZ, compare.only.mz, avoid.processing.mz, id.database, n.putative)
{	
	
	factors.list <- Experiment@Data@FactorList
	empty.samples <- which(lapply(factors.list,nrow)==0)
	if(length(empty.samples)!=0) factors.list <- factors.list[-empty.samples]
	
	if(length(factors.list)==1) {
		factors.list[[1]] <- cbind(factors.list[[1]],matrix(c(1:(length(factors.list[[1]]$ID))),nrow=length(factors.list[[1]]$ID)))
		colnames(factors.list[[1]])[ncol(factors.list[[1]])] <- "AlignID"	
	#}
	
		alignId <- lapply(factors.list,function(x){x$AlignID})
		N.groups <- max(unique(unlist(alignId)))	
	
		spectra.list <- apply(as.matrix(1:N.groups),1,function(Ng){
			align.iterator <- as.vector(which(lapply(alignId,function(x){
				if(length(intersect(Ng,x))!=0) {
					return(T)
					}else{
						return(F)
					}
				})==T))
				group.spectra <- lapply(as.matrix(align.iterator), function(x){ 
					list(spectra=convertMSPspectra(factors.list[[x]][which(alignId[[x]]==Ng),"Spectra"],maxMZ),time=factors.list[[x]][which(alignId[[x]]==Ng),"RT"]) 
					})
				common.spectra <- normalize(rowSums(do.call(cbind,lapply(group.spectra,function(x){x$spectra})))) 
				common.time <- mean(do.call(cbind,lapply(group.spectra,function(x){x$time})))
				list(spectra=common.spectra,time=common.time, nsamples=length(align.iterator))
		})
		
		spectra.matrix <- do.call(cbind,lapply(spectra.list,function(x){x$spectra}))
		time.vector <- unlist(lapply(spectra.list,function(x){x$time}))
		appear.vector <- unlist(lapply(spectra.list,function(x){x$nsamples}))
	}else{
		spectra.matrix <- apply(Experiment@Results@Alignment, 1, function(x) convertMSPspectra(x["Spectra"], maxMZ))
	}

	delete.mask <- rep(0,maxMZ)
	delete.mask[compare.only.mz] <- 1
	delete.mask[avoid.processing.mz] <- 0
	delete.mz <- which(delete.mask==0)
	spectra.matrix.compare <- spectra.matrix
	spectra.matrix.compare[delete.mz,] <- 0		

	cat("Constructing matrix database... \n")	
	spects <- get.spectra.matrix(id.database, maxMZ)
	spects[delete.mz,] <- 0		
	cat("Comparing spectra... \n")	
	mat.correlation <- cor.sinus(spectra.matrix.compare,spects)	

	identification.results <- NULL
	for(n.put in 1:n.putative)
	{
		identification.results.p <- t(apply(mat.correlation,1,function(x){
			out.mat <- matrix(0,ncol=2,nrow=1)
			if(!all(is.na(x)))
			{
				out.mat[,1] <- order(x, decreasing=T)[n.put]
				out.mat[,2] <- round(as.numeric(x[order(x, decreasing=T)][n.put])*100,digits=2)
			}else{
				out.mat[,1] <- NA
				out.mat[,2] <- NA
				}
			out.mat
		}))
		identification.results <- cbind(identification.results, identification.results.p)
		i.start <- (n.put-1)*2 + 1
		col.names <- colnames(identification.results)
		col.names[i.start:(i.start+1)] <- c(paste("DB.Id.",n.put,sep=""),paste("MatchFactor.",n.put,sep=""))	
		colnames(identification.results) <- col.names
	}	
	
	spectra.msp <- apply(as.matrix(1:ncol(spectra.matrix)),1,function(j){
			spectra <- spectra.matrix[,j]
			spectra.index <- which(spectra!=0) 
			spectra.pos <- spectra.index
			spectra.int <- round(spectra[spectra.index]*1000)
			spectra.text <- paste(sweep(as.matrix(spectra.pos),1,as.matrix(spectra.int),"paste.sp"), collapse=" ")
			spectra.text
		})	
		
	spectra.msp <- as.data.frame(spectra.msp)
	
	if(!is.null(Experiment@Results@Alignment$"tmean"))
	{
		time.vector <- Experiment@Results@Alignment$"tmean"
		foundIn.vector <- Experiment@Results@Alignment$"FoundIn"
		alignID.vector <- Experiment@Results@Alignment$"AlignID"
	}else{
		time.vector <- Experiment@Data@FactorList[[1]]$"RT"
		foundIn.vector <- rep(1, length(time.vector))
		alignID.vector <- 1:length(time.vector)
	}
		
	identification.metadata <- cbind(matrix(round(time.vector, digits=4),ncol=1), matrix(alignID.vector, ncol=1), matrix(foundIn.vector,ncol=1))
	identification.list <- cbind(spectra.msp,identification.metadata,identification.results)
	colnames(identification.list) <- c("Spectra","tmean","AlignID","FoundIn",colnames(identification.results))
	cat("Done! \n")
	id.list <- as.data.frame(identification.list, row.names=1:nrow(identification.list))
	id.list[,"Spectra"] <- as.character(id.list[,"Spectra"])
	id.list
}

get.spectra.matrix <- function(id.database, maxMZ)
{
	spects <- lapply(id.database, function(x) convertMSPspectra.dot(x$Spectra,maxMZ))
	spects <- do.call(cbind, spects)
	spects
}

setDBinMemory <- function(id.database, maxMZ)
{
	cat("Constructing matrix database... \n")	
	spects <- get.spectra.matrix(id.database@database, maxMZ)
	colnames(spects) <- unlist(lapply(id.database@database,function(x){x$Name}))
	cat("Done! \n")	
	spects
}

# seekSimilar <- function(Experiment, AlignId, id.database=mslib, n.similar=10)
# {
	# if(class(id.database)=="GCRMetDatabase") 
	# {
		# maxMZ <- max(Experiment@Results@Parameters@Alignment$mz.range)
		# spects <- get.spectra.matrix(id.database@database, maxMZ)
	# }else if(class(id.database)=="matrix"){
		# spects <- id.database
	# }else{
		# stop("Database format not recognized")
	# }
	
	# maxMZ <- max(Experiment@Results@Parameters@Alignment$mz.range)
	# index <- which(Experiment@Results@Identification[,"AlignID"]==AlignId)
	# MSP.spect.emp <- Experiment@Results@Identification[index,"Spectra"]
	# empiric.spectra <- convertMSPspectra(MSP.spect.emp,maxMZ)
	# empiric.spectra <- empiric.spectra/max(empiric.spectra)*1000
	
	# ##S'ha de multiplicar la matriu de spects i empiric.spectra per el mz.range etc.
	# mat.correlation <- cor.sinus(empiric.spectra,spects)
	
	# corr.vector <- as.vector(mat.correlation)
	# sim.index <- order(corr.vector, decreasing=T)[1:n.similar]
	# sim.corr <- round(corr.vector[sim.index]*100, digits=2)
	
	# sim.list <- as.data.frame(NULL)
	# if(class(id.database)=="GCRMetDatabase") 
	# {
		# sim.names <- unlist(lapply(id.database@database[sim.index], function(x){x$Name}))
		# sim.cas <- unlist(lapply(id.database@database[sim.index], function(x){x$CAS}))
		# sim.form <- unlist(lapply(id.database@database[sim.index], function(x){x$Formula}))
		
		# sim.list <- as.data.frame(matrix(c(sim.index,sim.names,sim.corr,sim.cas,sim.form), nrow=n.similar))
		# colnames(sim.list) <- c("DB.Id","Name","MatchFactor","CAS","Formula")
		
	# }else if(class(id.database)=="matrix"){
		# sim.names <- colnames(spects)[sim.index]
		# sim.list <- as.data.frame(matrix(c(sim.index,sim.names,sim.corr), nrow=n.similar))
		# colnames(sim.list) <- c("DB.Id","Name","MatchFactor")
	
	# }

	# sim.list
# }
