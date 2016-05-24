# msap - Statistical analysis for Methilation-Sensitive Amplification Polimorphism data
# Author: Andrés Pérez-Figueroa (anpefi@uvigo.es)




msap <- function(datafile, name=datafile, no.bands="u", nDec=4, meth=TRUE, rm.redundant=TRUE, rm.monomorphic=TRUE, do.pcoa=TRUE, do.shannon=TRUE, do.amova=TRUE, do.pairwisePhiST=FALSE, do.cluster=TRUE, use.groups=NULL, do.mantel=FALSE, np.mantel=1000, loci.per.primer=NULL, error.rate.primer=NULL, uninformative=TRUE){
	msapVer <- as.character(packageVersion("msap"))
	cat(paste("\nmsap ", msapVer, " - Statistical analysis for Methylation-Sensitive Amplification Polimorphism data\n"),sep="")
	 
	 ########## CHECKING PARAMETERS ############
	
	#check no.bands parameter
	#uninformative (default)
	if (is.character(no.bands)){
		if(no.bands=="u") uninformative <- TRUE #uninformative (default)
		else if(no.bands=="h") uninformative <- FALSE #assuming hypermethilated target
		else{
			cat("\nWARNING: no.bands argument wrong. Using default value (no.bands=\"u\")\n\n")
			uninformative <- TRUE #uninformative (default)
		}
	}
	else cat("\nWARNING: no.bands argument wrong. Using default value (no.bands=\"u\")\n\n")
	
	#check loci.per.primer
	n.primer <- 1 #default
	if(is.vector(loci.per.primer)){
		if(is.numeric(loci.per.primer)){
			n.primer <- length(loci.per.primer)
			
		}
		else cat("\nWARNING: fragments.per.primer argument wrong. Combining all fragments.\n\n")
	}
	
	
	#check error.rate.primer
	if(is.vector(error.rate.primer)){
		if(is.numeric(error.rate.primer)){
			if(length(error.rate.primer)!=n.primer){
				cat("\nWARNING: error.rate.primer argument doesn't match number of primers in loci.per.primer. Using default value (0.05) for all primer combinations\n\n")
				error.rate.primer<-rep(0.05, n.primer)			
			}
			
		}
		else{
			 cat("\nWARNING: error.rate.primer argument not valid. Using default value (0.05) for all primer combinations\n\n")
			error.rate.primer<-rep(0.05, n.primer)		
		}
	}

  ## CHECKING AND READING DATAFILE ##
  
	#File exists?
	if(!file.exists(datafile)) stop(datafile," unreachable. File does not exist or you have no permissions to open it\n", call.=F)
	#File has commas?
	con<-file(datafile,"r")
	str <- readLines(con)
  close(con)
	if(sum(unlist(strsplit(str[1], "")) == ",")==0) stop("Doesn't find \",\" in ",datafile,": Please use a .csv file with fields sepparated by commas")
	else {
	  # have all files the same number of fields?
	  if(!(sum(unlist(strsplit(str[1], "")) == ",") == sum(unlist(strsplit(str[1:length(str)-1], "")) == ",")/(length(str)-1))) 
	    stop ("Bad formatting. There are not the same number of columns in all the rows of the file, please check ", datafile," carefully (it should have the same number of commas in all its rows)")
	}
  
	cat("\nReading ", datafile,".....")
	data <- read.csv(datafile, header=TRUE)
	
	if(meth){
    #Check enzymes
	  temp <- factor(data[,3])
	
	  if(!is.element("HPA",levels(temp))) stop("HPA missing. The third column should include both HPA and MSP (capital letters)")
	  if(!is.element("MSP",levels(temp))) stop("MSP missing. The third column should include both HPA and MSP (capital letters)")
	  if(!table(temp)["MSP"]==table(temp)["HPA"]) stop("There are not equal number of rows for HPA and MSP, please check ", datafile)
	}
	#Check other values than 0/1
	vals<-levels(factor(as.matrix(data[,4:length(data[1,])])))
	if(!(length(vals)==2 && is.element("0", vals) && is.element("1",vals))) stop("There is one (or more) unusual values in the data matrix. Remember only 0 and 1 are allowed. Please, check ",datafile)
	
	cat(" Ok!\n")
  
  ## SETTING DATA ##
  #get the loci names
	locus <- colnames(data)[4:length(data[1,])]
	cat("Number of loci: ",length(locus),"\n")
	#sorting
	data <- data[with(data, order(data[,1],data[,2],data[,3])),]
	#set primer mask
	if(is.null(loci.per.primer)) loci.per.primer <- c(length(locus)) #Avoids crash
	if(is.null(error.rate.primer)) error.rate.primer <- rep(0.05, length(loci.per.primer))
	primer.mask <- rep (1:length(loci.per.primer), loci.per.primer)
	
	
	#If user just want to use some of the groups in the datafile
	#Should be provided with the same name that in file
	if(is.vector(use.groups)){
		if(is.character(use.groups)){
			data <- data[which(is.element(data[,1],use.groups)),]
			cat("\nUsing a subset of the data:\n")
			print(levels(factor(data[,1]))) 
		}
		else cat("\nWARNING: use.groups argument wrong. Using all groups.\n\n")
	}
	
	
	if(meth){ #meth=TRUE -> MSAP
	
		############ FRAGMENT CLASSIFICATION ##############
		inds <- as.character(data[data[,3]=="HPA",2])
		groups <- factor(data[data[,3]=="HPA",1])  #Needed to refactor here 
		ntt <- length(levels(groups))
		cat("Number of samples/individuals: ",length(data[,1])/2,"\n")
		cat("Number of groups/populations: ",ntt,"\n")
		cat("Number of primer combinations: ",n.primer,"\n")
		cat("Loci per primer combinations ",loci.per.primer,"\n")
		cat("Error rates per primer combination: ",error.rate.primer,"\n")
		dataMSP <- data[data[,3]=="MSP",4:length(data[1,])]
		dataHPA <- data[data[,3]=="HPA",4:length(data[1,])]
		dataMIX <- dataHPA*10+dataMSP
		#check methStatus within different primer combinations (since 1.1.0)
		Met <- NULL
		for (i in 1:n.primer){
			temp <- apply(dataMIX[,which(primer.mask==i)], 2, methStatusEval, error=error.rate.primer[i], uninformative=uninformative)
			Met <- c(Met,temp)
			cat("Primer: ",i,"\n")
			cat("--Number of Methylation-Susceptible Loci (MSL): ",length(which(temp)),"\n")
			cat("--Number of No Methylated Loci (NML): ",length(which(!temp)),"\n\n")
		}
		#Met <- apply(dataMIX, 2, methStatusEval, type=uninformative)
		MSL <- which(Met)
		NML <- which(!Met)
		MSL.nloci <- length(Met[MSL])
		NML.nloci <- length(Met[NML])
		cat("All combinations: \n")
		cat("Number of Methylation-Susceptible Loci (MSL): ",MSL.nloci,"\n")
		cat("Number of No Methylated Loci (NML): ",NML.nloci,"\n\n")
				
		
		if (uninformative){
		#This transformation assumes that HPA-/MSP- pattern is uninformative as it could represent full methylation of cytosines in the target or that target is missing by mutation. 
		#Herrera & Bazaga 2010
			
			if(MSL.nloci>0) dataMSL <- dataMIX[,MSL]
			if(NML.nloci>0) dataNML <- dataMIX[,NML]
			if(MSL.nloci>0){
				dataMSL[dataMSL==0] <- NA
				dataMSL[dataMSL==1] <- 2
				dataMSL[dataMSL==11] <- 1
				dataMSL[dataMSL==10] <- 2
			}
			if(NML.nloci>0) dataNML <- ifelse(dataNML==0, 1, 2)
		
			if(MSL.nloci>0) matM <- as.matrix(dataMSL)
			if(NML.nloci>0) matN <- as.matrix(dataNML)
		}
		else 
		{
			#This transformation assumes that HPA-/MSP- pattern represents full methylation (hypermethilation) of cytosines in the target ignoring possible genetic differences. Thus, all bands are scored as methylation.
			#Lu et al. 2008; Gupta et al. 2012

			dataMIXb <- ifelse(dataMIX==11, 1, 2)  #changed in 1.1.2 for consistency to 1:unmmethylated 2:methylated. Previous values were not affect to the results from the analysis present at the moment
			if(MSL.nloci>0) matM <- as.matrix(dataMIXb[,MSL])
			if(NML.nloci>0) matN <- as.matrix(dataMIXb[,NML])
		}
		
		if(MSL.nloci>0) PolyM <- apply(matM, 2, polymorphic)
		if(NML.nloci>0) PolyN <- apply(matN, 2, polymorphic)
		if(MSL.nloci>0) MSL.ploci <- length(which(PolyM))
		else MSL.ploci <- 0
		if(NML.nloci>0) NML.ploci <- length(which(PolyN))
		else NML.ploci <-0
		if(MSL.nloci>0) cat("Number of polymorphic MSL: ",MSL.ploci," (",format(MSL.ploci/MSL.nloci*100,digits=1),"% of total MSL)\n")
		if(NML.nloci>0) cat("Number of polymorphic NML: ",NML.ploci," (",format(NML.ploci/NML.nloci*100,digits=1),"% of total NML)\n\n")
		
		if(rm.monomorphic){
			if(MSL.ploci>0) matM <- matM[,PolyM]
			if(NML.ploci>0) matN <- matN[,PolyN]
			NML.nloci <- NML.ploci
			MSL.nloci <- MSL.ploci 
		}
		
		#save transformed files (MSL and NML)
		if(MSL.nloci>0){
			cat("- Saving transformed matrix for MSL in file: ",paste(name,"-MSL-transformed.csv",sep=""),"\n")
		 	write.csv(data.frame(groups,inds,matM), file=paste(name,"-MSL-transformed.csv",sep=""), row.names=FALSE)
		}
		if(NML.nloci>0){
			cat("- Saving transformed matrix for NML in file: ",paste(name,"-NML-transformed.csv",sep=""),"\n")
			write.csv(data.frame(groups,inds,matN), file=paste(name,"-NML-transformed.csv",sep=""), row.names=FALSE)
		}
		if(MSL.nloci>1) MSL.I <- apply(matM, 2, shannon)
		if(NML.nloci>1) NML.I <- apply(matN, 2, shannon)
		if(MSL.nloci>1)cat("\nShannon's Diversity Index \n")
		if(MSL.nloci>1) cat("MSL: I = ", mean(MSL.I, na.rm=T),"  (SD: ", sd(MSL.I, na.rm=T),")\n")
		if(NML.nloci>1){
		cat("NML: I = ", mean(NML.I, na.rm=T),"  (SD: ", sd(NML.I, na.rm=T),")\n")
		wt<-wilcox.test(MSL.I,NML.I)
		pval<-ifelse(wt$p.value<0.0001, "P < 0.0001", paste("P = ",wt$p.value))
		cat(wt$method,": ", names(wt$statistic)," = ",wt$statistic[[1]]," (",pval,")\n")
		#png(filename=paste(name,"Shannon.png", sep='-'), width=680, height=680)
		#boxplot(MSL.I, NML.I, names=c("MSL","NSL"), ylab="Shannon's I")
		#dev.off()
		}
	
		### MSL 
		if(MSL.nloci>1){
		
			cat("\n\n*****************************\nAnalysis of MSL\n")
			meth.rep <- repMet(dataMIX[,MSL], groups, nDec)
		  cat("\n\n")


		  DM<- dist(matM, method="euclidean", upper=T, diag=T)
			if(do.cluster){
				# cluster
				DM_copy <- DM
				attr(DM_copy, "Labels") <- inds
				np <- table(groups)[]
				MSL.cluster <-nj(DM_copy) 
				tipCol<-rep(rainbow(length(np)+1)[-2], np)
				ecol<-unlist(lapply(MSL.cluster$edge[,2], edgeCol, length(MSL.cluster$tip.label), tipCol))
				cat("- Saving clustering tree figure for MSL in file: ", paste(name,"MSL-NJ.png", sep='-'),".....")
				png(filename=paste(name,"MSL-NJ.png", sep='-'))
				plot.phylo(MSL.cluster, tip.color=tipCol, use.edge.length=T, edge.color=ecol, edge.width=3, show.tip.label=T, main="MSL")
				legend("bottomright", as.character(levels(groups)), col=rainbow(length(np)+1)[-2], lwd=3)
				dev.off()
        cat("Ok!\n")
			}
		
			#PCoA
			if(do.pcoa)pcoa(DM, groups, inds, name, "MSL")
			if(do.amova){
				#AMOVA
				cat("\nPerforming AMOVA\n")
     		diffAmova(DM, groups, nDec, do.pairwisePhiST)
			}
			
      DM.MSL <- DM
		} #end if MSL.loci>0
		else{
		cat("There are not polymorphic MSL. Diversity Analysis skipped.\n")
	}		
	
	} #end if meth
	else{  #meth=FALSE -> AFLP
		cat("Analysis for standard AFLP loci !!!\n\n")
		inds<-as.character(data[,2])
		groups <- factor(data[,1]) #Here, all rows are indepedent samples
		ntt <- length(levels(groups))
		matN <- as.matrix(data[,4:length(data[1,])])
		NML.nloci <- length(locus)
		matN[matN==1]<-2
		matN[matN==0]<-1 
		cat("Number of samples/individuals: ",length(inds),"\n")
		cat("Number of groups/populations: ",ntt,"\n")
		PolyN <- apply(matN, 2, polymorphic)
		NML.ploci <- length(which(PolyN))
		cat("Number of polymorphic AFLP: ",NML.ploci," (",format(NML.ploci/NML.nloci*100,digits=1),"% of total)\n\n")
		NML.nloci <- NML.ploci
		
		
		DM<- dist(matN, method="euclidean", upper=T, diag=T)
		if(do.cluster){
			# cluster
			DM_copy <- DM
			attr(DM_copy, "Labels") <- inds
			np <- table(groups)[]
			MSL.cluster <-nj(DM_copy) 
			tipCol<-rep(rainbow(length(np)+1)[-2], np)
			ecol<-unlist(lapply(MSL.cluster$edge[,2], edgeCol, length(MSL.cluster$tip.label), tipCol))
			cat("- Saving clustering tree figure for AFLP in file: ", paste(name,"AFLP-NJ.png", sep='-'),".....")
			png(filename=paste(name,"AFLP-NJ.png", sep='-'))
			plot.phylo(MSL.cluster, tip.color=tipCol, use.edge.length=T, edge.color=ecol, edge.width=3, show.tip.label=T, main="AFLP")
			legend("bottomright", as.character(levels(groups)), col=rainbow(length(np)+1)[-2], lwd=3)
			dev.off()
      cat("Ok!\n")
			
		}
		
		#PCoA
		if(do.pcoa) pcoa(DM, groups, inds, name, "AFLP")
		
		DM.AFLP<-DM
		if(do.amova){
			#AMOVA
			cat("\nPerforming AMOVA\n")
			diffAmova(DM, groups, nDec, do.pairwisePhiST)
		}
		
		NML.nloci <- 0  #Added to skip NML analysis
	}
	#### NML
	if(NML.nloci>1){	
		if(meth) cat("\n\n*****************************\nAnalysis of NML\n\n")
		DM<- dist(matN, method="euclidean", upper=T, diag=T)
		if(do.cluster){
			# cluster
			DM_copy <- DM
			attr(DM_copy, "Labels") <- inds
			np <- table(groups)[]
			MSL.cluster <-nj(DM_copy) 
      tipCol<-rep(rainbow(length(np)+1)[-2], np)
			ecol<-unlist(lapply(MSL.cluster$edge[,2], edgeCol, length(MSL.cluster$tip.label), tipCol))
			cat("- Saving clustering tree figure for NML in file: ", paste(name,"NML-NJ.png", sep='-'),".....")
			png(filename=paste(name,"NML-NJ.png", sep='-'))
			plot.phylo(MSL.cluster, tip.color=tipCol, use.edge.length=T, edge.color=ecol, edge.width=3, show.tip.label=T, main="NML")
			legend("bottomright", as.character(levels(groups)), col=rainbow(length(np)+1)[-2], lwd=3)
			dev.off()
      
      DM.NML <- DM
      cat("Ok!\n")
			
		}
		#PCoA
		if(do.pcoa) pcoa(DM, groups, inds, name, "NML")
		
		if(do.amova){
			#AMOVA
			cat("\nPerforming AMOVA\n")
			diffAmova(DM, groups, nDec, do.pairwisePhiST)
		}
		
		if(do.mantel){
			mtl <-mantel.randtest(DM.MSL, DM.NML, np.mantel)
			pval<-ifelse(mtl$pvalue<0.0001, "P < 0.0001", paste("P = ",format(mtl$pvalue, digits=nDec)))
			cat("\nMantel test (MSL/NML): r = ", mtl$obs," (",pval,"; nperm= ",mtl$rep,")\n")
		}
	}
	else{
		if(meth) cat("There are not polymorphic NML. Diversity Analysis skipped.\n")
	}
	
  #Storing some interesting item in a list to be returned (as suggested by C. Herrera)
	res <- list(
    groups = if(exists("groups")) {groups} else {NULL},
    patterns = if(meth && exists("meth.rep")) {meth.rep} else {NULL},
    transformed.MSL = if(meth && MSL.nloci>0) {data.frame(groups,inds,matM)} else {NULL},
    transformed.NML =  if(meth && NML.nloci>0) {data.frame(groups,inds,matN)} else {NULL},
    DM.MSL = if(exists("DM.MSL")) {DM.MSL} else {NULL},
    DM.NML = if(exists("DM.NML")) {DM.NML} else {NULL},
    DM.AFLP = if(exists("DM.AFLP")) {DM.AFLP} else {NULL}
     )
  
	 cat("Done!\n")
   invisible(res)
}


edgeCol <- function(x, n, pal){
if(x > n) col <- "black"
else col <-pal[x]
return(col)
}





