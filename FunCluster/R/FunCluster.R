.packageName <- "FunCluster"

ref.list <- NULL
up <- NULL
down <- NULL
genes.lst <- NULL
f.version <- "1.09"

.First.lib <- function(lib, pkg, ...)
{
  verbose <- .Options$Hverbose
  if(!length(verbose) || verbose)

cat(paste("\nThis is FunCluster package version ",f.version," maintained by Corneliu Henegar.\n\n",
	"FunCluster(wd='', org='HS', go.direct=FALSE, clusterm='cc', compare='common.correl.genes',\n",
		      "\tcorr.met='greedy', corr.th=0.85, two.lists=TRUE, restrict=FALSE, alpha=0.05,\n",
		      "\tlocation=FALSE, details=FALSE)\n\n",sep=''))
  invisible()
}


#############################################################################################
#
# 1. Function FunCluster() -> Assures the global control of the treatment process
#
#############################################################################################


FunCluster <- function(wd="", org="HS", go.direct=FALSE, clusterm="cc", compare="common.correl.genes", 
			corr.met="greedy", corr.th=0.85, two.lists=TRUE, restrict=FALSE, alpha=0.05, 
			location=FALSE, details=FALSE){

	# "two.lists" means 2 lists of genes instead of one (i.e. UP & DOWN)
	# "restrict" means restrict the enrichment calculus to a list of genes (in oposition with the whole genome by default)

	cat(paste("\n\tFunCluster started at: ",date(),sep=""))
	cat(paste("\n\t\tUsing annotations updated on: ",annot.date,sep=""))

	if(wd != ""){setwd(wd)} # set working directory
	# create results directory
	results.dir <- paste("Results_",format(Sys.time(), "%Y_%b_%d_%H-%M-%S"),sep="")
	dir.create(paste(getwd(),"/",results.dir,sep=""))

# GO annotations: direct vs indirect
	if(go.direct == FALSE){
		godir <- "IND"
	}else{
		godir <- "DIR"
	}
# Load data files	

	wd <- getwd()	# working directory

	# EntrezGene ID's	
	if(org == "HS"){locus.name <- HS.locus.name}
	if(org == "MM"){locus.name <- MM.locus.name}
	if(org == "SC"){locus.name <- SC.locus.name}

	if(restrict == TRUE){	# specifing a reference list of genes for the calculation of annotations enrichment
		if(is.null(ref.list)){ref.list <- read.table(paste(wd,"/","ref.txt",sep=""),sep="\t")}

		ref.lst <- filter.genes(restrict=TRUE,ref.list=ref.list)
		ref.list <- ref.lst$ref.list

		rm(ref.lst)

	}else{
		ref.list <- NULL
	}

	if(two.lists == TRUE){	# when 2 lists of genes are analysed in oposition (UP & DOWN)

		if(is.null(up)){up <- read.table(paste(wd,"/","up.txt",sep=""),sep="\t")}	# genes UP
		if(is.null(down)){down <- read.table(paste(wd,"/","down.txt",sep=""),sep="\t")}	# genes DOWN

		up.down <- filter.genes(up=up,down=down,two.lists=TRUE)
		up <- up.down$up
		down <- up.down$down
		up.frame <- up.down$up.frame

		down.frame <- up.down$down.frame


		if(compare == "common.correl.genes"){
			up.correl <- gene.correl(up.frame,corr.th,method=corr.met,"UP")

			down.correl <- gene.correl(down.frame,corr.th,method=corr.met,"DOWN")

		}else{
			up.correl <- NULL
			down.correl <- NULL
		}

		if(location == TRUE){
			if(org == "HS"){locus.cyto <- HS.locus.cyto}
			if(org == "MM"){locus.cyto <- MM.locus.cyto}
			if(org == "SC"){locus.cyto <- SC.locus.cyto}
			# LocusLink cytoband

			location.analysis(list.anal=up,ref.list=ref.list,nom="UP",cyto=locus.cyto,locus.name=locus.name)	# cytoband analysis
			location.analysis(list.anal=down,ref.list=ref.list,nom="DOWN",cyto=locus.cyto,locus.name=locus.name)
		}


	}else{
		if(is.null(genes.lst)){genes.lst <- read.table(paste(wd,"/","genes.txt",sep=""),sep="\t")}	# a unique list of genes

		genes <- filter.genes(genes.lst=genes.lst,two.lists=FALSE)
		genes.lst <- genes$genes.lst
		genes.frame <- genes$genes.frame


		if(compare == "common.correl.genes"){
			genes.correl <- gene.correl(genes.frame,corr.th,"LIST")

		}else{
			genes.correl <- NULL
		}

		if(location == TRUE){
			if(org == "HS"){locus.cyto <- HS.locus.cyto}
			if(org == "MM"){locus.cyto <- MM.locus.cyto}
			if(org == "SC"){locus.cyto <- SC.locus.cyto}
			# LocusLink cytoband
			location.analysis(list.anal=genes.lst,ref.list=ref.list,nom="list",cyto=locus.cyto,locus.name=locus.name)	# cytoband analysis	
		}
	}



# KEGG Annotations

	terms.name <- KEGG.terms.name	# KEGG structure
	# KEGG annotations file
	if(org == "HS"){file.annot <- HS.KEGG.file.annot}
	if(org == "MM"){file.annot <- MM.KEGG.file.annot}
	if(org == "SC"){file.annot <- SC.KEGG.file.annot}
	taxoname <- "KEGG"


	if(two.lists == TRUE){
		main.loop(file.annot=file.annot,taxoname=taxoname,terms.name=terms.name,up=up,down=down,
		results.dir=results.dir,clusterm=clusterm,alpha=alpha,locus.name=locus.name,compare=compare,
		up.frame=up.frame,down.frame=down.frame,up.correl=up.correl,down.correl=down.correl,restrict=restrict,ref.list=ref.list,details=details)
	}else{
		main.loop(file.annot=file.annot,taxoname=taxoname,terms.name=terms.name,results.dir=results.dir,
		clusterm=clusterm,alpha=alpha,locus.name=locus.name,compare=compare,genes.lst=genes.lst,
		genes.frame=genes.frame,genes.correl=genes.correl,restrict=restrict,ref.list=ref.list,details=details)
	}


# Gene Ontology Annotations

	go.name <- c("GO Biological Process","GO Cellular Component","GO Molecular Function") # names of GO branches
	terms.name <- GO.terms.name # GO structure  

	for(i in 1:3){

	# GO annotations file HS
		if(org == "HS" && godir == "DIR" && i == 1){file.annot <- HS.GO.DIR.BP.file.annot}
		if(org == "HS" && godir == "DIR" && i == 2){file.annot <- HS.GO.DIR.CC.file.annot}
		if(org == "HS" && godir == "DIR" && i == 3){file.annot <- HS.GO.DIR.MF.file.annot}
		if(org == "HS" && godir == "IND" && i == 1){file.annot <- HS.GO.IND.BP.file.annot}
		if(org == "HS" && godir == "IND" && i == 2){file.annot <- HS.GO.IND.CC.file.annot}
		if(org == "HS" && godir == "IND" && i == 3){file.annot <- HS.GO.IND.MF.file.annot}
	# GO annotations file MM
		if(org == "MM" && godir == "DIR" && i == 1){file.annot <- MM.GO.DIR.BP.file.annot}
		if(org == "MM" && godir == "DIR" && i == 2){file.annot <- MM.GO.DIR.CC.file.annot}
		if(org == "MM" && godir == "DIR" && i == 3){file.annot <- MM.GO.DIR.MF.file.annot}
		if(org == "MM" && godir == "IND" && i == 1){file.annot <- MM.GO.IND.BP.file.annot}
		if(org == "MM" && godir == "IND" && i == 2){file.annot <- MM.GO.IND.CC.file.annot}
		if(org == "MM" && godir == "IND" && i == 3){file.annot <- MM.GO.IND.MF.file.annot}
	# GO annotations file SC
		if(org == "SC" && godir == "DIR" && i == 1){file.annot <- SC.GO.DIR.BP.file.annot}
		if(org == "SC" && godir == "DIR" && i == 2){file.annot <- SC.GO.DIR.CC.file.annot}
		if(org == "SC" && godir == "DIR" && i == 3){file.annot <- SC.GO.DIR.MF.file.annot}
		if(org == "SC" && godir == "IND" && i == 1){file.annot <- SC.GO.IND.BP.file.annot}
		if(org == "SC" && godir == "IND" && i == 2){file.annot <- SC.GO.IND.CC.file.annot}
		if(org == "SC" && godir == "IND" && i == 3){file.annot <- SC.GO.IND.MF.file.annot}



		taxoname <- go.name[i]	# name of the GO branch considered for gene annotations

		if(two.lists == TRUE){
			main.loop(file.annot=file.annot,taxoname=taxoname,terms.name=terms.name,up=up,down=down,
			results.dir=results.dir,clusterm=clusterm,alpha=alpha,locus.name=locus.name,compare=compare,
			up.frame=up.frame,down.frame=down.frame,up.correl=up.correl,down.correl=down.correl,restrict=restrict,ref.list=ref.list,details=details)
		}else{
			main.loop(file.annot=file.annot,taxoname=taxoname,terms.name=terms.name,results.dir=results.dir,clusterm=clusterm,alpha=alpha,locus.name=locus.name,compare=compare,genes.lst=genes.lst,
			genes.frame=genes.frame,genes.correl=genes.correl,restrict=restrict,ref.list=ref.list,details=details)
		}		
	}


	cat(paste("\n\tEnd  of treatment at: ",date(),"\n",sep=""))	
	rm()


}


#############################################################################################
#
# 2. Function MAIN.LOOP() -> Recursive treatment of data
#
#############################################################################################

main.loop <- function(file.annot,taxoname,terms.name,up=NULL,down=NULL,results.dir,clusterm,alpha,
			locus.name,compare,up.frame=NULL,down.frame=NULL,up.correl=NULL,down.correl=NULL,
			genes.lst=NULL,genes.frame=NULL,genes.correl=NULL,restrict=FALSE,ref.list=NULL,details){

	annotations <- preliminary(file.annot=file.annot,taxoname=taxoname,restrict=restrict,ref.list=ref.list)	# preliminary treatment of annotations file


if(!is.null(up) && !is.null(down) && is.null(genes.lst) && !is.null(annotations)){
	up.data <- annotate(annotations,up,"UP",terms.name,taxoname)
	down.data <- annotate(annotations,down,"DOWN",terms.name,taxoname)


	up.data <- pterm(up.data,taxoname,"UP")	
	down.data <- pterm(down.data,taxoname,"DOWN")		

	if(details == TRUE){
		resterm(up.data,taxoname,"UP",locus.name,results.dir=results.dir)
		resterm(down.data,taxoname,"DOWN",locus.name,results.dir=results.dir)
	}

	if(clusterm == "cs"){
		c.up <- clusterm.cs(up.data,taxoname,"UP",alpha,compare,up.frame,up.correl)
		c.down <- clusterm.cs(down.data,taxoname,"DOWN",alpha,compare,down.frame,down.correl)

		c.up <- cluster.c(c.up,up.data,annotations,taxoname,"UP")
		c.down <- cluster.c(c.down,down.data,annotations,taxoname,"DOWN")

		c.up <- pclust(up.data,c.up,annotations,taxoname,"UP")		
		c.down <- pclust(down.data,c.down,annotations,taxoname,"DOWN")		

		resclust(up.data,c.up,taxoname,"UP","CS",locus.name,results.dir=results.dir)
		resclust(down.data,c.down,taxoname,"DOWN","CS",locus.name,results.dir=results.dir)
	}else if(clusterm == "cc"){
		c.up <- clusterm.cc(up.data,taxoname,"UP",alpha,compare,up.frame,up.correl)
		c.down <- clusterm.cc(down.data,taxoname,"DOWN",alpha,compare,down.frame,down.correl)

		c.up <- cluster.c(c.up,up.data,annotations,taxoname,"UP")
		c.down <- cluster.c(c.down,down.data,annotations,taxoname,"DOWN")

		c.up <- pclust(up.data,c.up,annotations,taxoname,"UP")		
		c.down <- pclust(down.data,c.down,annotations,taxoname,"DOWN")		

		resclust(up.data,c.up,taxoname,"UP","CC",locus.name,results.dir=results.dir)
		resclust(down.data,c.down,taxoname,"DOWN","CC",locus.name,results.dir=results.dir)
	}

}else if(!is.null(genes.lst) && !is.null(annotations)){


	genes.data <- annotate(annotations,genes.lst,"LIST",terms.name,taxoname)
	genes.data <- pterm(genes.data,taxoname,"LIST")



	if(details == TRUE){
		resterm(genes.data,taxoname,"LIST",locus.name,results.dir=results.dir)
	}

	if(clusterm == "cs"){
		c.genes <- clusterm.cs(genes.data,taxoname,"LIST",alpha,compare,genes.frame,genes.correl)

		c.genes <- cluster.c(c.genes,genes.data,annotations,taxoname,"LIST")

		c.genes <- pclust(genes.data,c.genes,annotations,taxoname,"LIST")		

		resclust(genes.data,c.genes,taxoname,"LIST","CS",locus.name,results.dir=results.dir)


	}else if(clusterm == "cc"){
		c.genes <- clusterm.cc(genes.data,taxoname,"LIST",alpha,compare,genes.frame,genes.correl)

		c.genes <- cluster.c(c.genes,genes.data,annotations,taxoname,"LIST")

		c.genes <- pclust(genes.data,c.genes,annotations,taxoname,"LIST")		

		resclust(genes.data,c.genes,taxoname,"LIST","CC",locus.name,results.dir=results.dir)


	}
}


	cat(paste("\n\t",taxoname," annotations treatment finished... ",date(),sep=""))
	rm()

}




#############################################################################################
#
# 3. Function FILTER.GENES() -> Filtering genes UP vs DOWN
#
#############################################################################################


filter.genes <- function(up=NULL,down=NULL,genes.lst=NULL,two.lists=TRUE,restrict=FALSE,ref.list=NULL){
	
	if(restrict == TRUE && !is.null(ref.list)){
		
		cat(paste("\n\tFiltering reference list started... ",format(Sys.time(), "%X"),sep=""))
		
		genes.locus <- NULL
		
		for(i in 1:length(ref.list[,1])){
			if(regexpr(";",as.character(ref.list[i,1])) != -1){
				x <- strsplit(gsub(" ","",as.character(ref.list[i,1])),";")[[1]]
							
				for(j in 1:length(x)){
					genes.locus <- c(genes.locus,x[j])				
				}
			}else{
				genes.locus <- c(genes.locus,as.character(ref.list[i,1]))				
			}		
		}
		
		genes.locus <- levels(as.factor(genes.locus))
		genes.locus <- genes.locus[genes.locus != ""]
		names(genes.locus) <- genes.locus
		
		return(list(ref.list=genes.locus))
	}



if(two.lists == FALSE && !is.null(genes.lst)){	

	cat(paste("\n\tFiltering genes list started... ",format(Sys.time(), "%X"),sep=""))

	
	genes.matrix <- NULL
	genes.locus <- NULL
	
	for(i in 1:length(genes.lst[,1])){
		if(regexpr(";",as.character(genes.lst[i,1])) != -1){
			x <- strsplit(gsub(" ","",as.character(genes.lst[i,1])),";")[[1]]
						
			for(j in 1:length(x)){
				genes.locus <- c(genes.locus,x[j])				
				genes.matrix <- rbind(genes.matrix,as.double(genes.lst[i,2:length(genes.lst)]))	
			}
		}else{
			genes.locus <- c(genes.locus,as.character(genes.lst[i,1]))				
			genes.matrix <- rbind(genes.matrix,as.double(genes.lst[i,2:length(genes.lst)]))
		}		
	}
	

	names(genes.matrix) <- 1:length(genes.matrix)

	genes.frame <- data.frame(genes.locus,genes.matrix)
	
	genes.locus <- levels(as.factor(genes.locus))
	
	genes.frame.new <- NULL
	
	for(i in 1:length(genes.locus)){
		
		datas <- genes.frame[genes.frame[,1] == genes.locus[i],][,2:length(genes.frame)]
		
		if(is.vector(datas) != TRUE){
			datas <- apply(datas,2,function(x) mean(as.double(x),na.rm=TRUE))
			datas[is.nan(datas)] <- NA
		}else{
			datas <- as.double(datas)
		}
		
		genes.frame.new <- rbind(genes.frame.new,c(genes.locus[i],datas))
		
	}
	
	genes.frame <- as.data.frame(genes.frame.new)
	
	return(list(genes.lst=genes.locus,genes.frame=genes.frame))



}else if(two.lists == TRUE && !is.null(up) && !is.null(down)){
	
	cat(paste("\n\tFiltering genes UP/DOWN started... ",format(Sys.time(), "%X"),sep=""))
	
	up.matrix <- NULL
	up.locus <- NULL
	
	for(i in 1:length(up[,1])){
		if(regexpr(";",as.character(up[i,1])) != -1){
			x <- strsplit(gsub(" ","",as.character(up[i,1])),";")[[1]]
						
			for(j in 1:length(x)){
				up.locus <- c(up.locus,x[j])				
				up.matrix <- rbind(up.matrix,as.double(up[i,2:length(up)]))	
			}
		}else{
			up.locus <- c(up.locus,as.character(up[i,1]))				
			up.matrix <- rbind(up.matrix,as.double(up[i,2:length(up)]))
		}		
	}
	

	names(up.matrix) <- 1:length(up.matrix)

	up.frame <- data.frame(up.locus,up.matrix)

	
	down.matrix <- NULL
	down.locus <- NULL
		
	for(i in 1:length(down[,1])){
		if(regexpr(";",as.character(down[i,1])) != -1){
			x <- strsplit(gsub(" ","",as.character(down[i,1])),";")[[1]]
							
			for(j in 1:length(x)){
				down.locus <- c(down.locus,x[j])				
				down.matrix <- rbind(down.matrix,as.double(down[i,2:length(down)]))	
			}
		}else{
			down.locus <- c(down.locus,as.character(down[i,1]))				
			down.matrix <- rbind(down.matrix,as.double(down[i,2:length(down)]))
		}		
	}
		
	names(down.matrix) <- 1:length(down.matrix)	

	down.frame <- data.frame(down.locus,down.matrix)
	

	

	up.locus <- levels(as.factor(up.locus))
	names(up.locus) <- up.locus
	down.locus <- levels(as.factor(down.locus))
	names(down.locus) <- down.locus
	
	x <- sort(up.locus[down.locus])
        
	if(length(x) > 0){
	for(i in 1:length(x)){
		up.frame <- up.frame[up.frame[,1] != x[i],]
		up.locus <- up.locus[up.locus != x[i]]
		down.frame <- down.frame[down.frame[,1] != x[i],]
		down.locus <- down.locus[down.locus != x[i]]
	}}
	
	rm(x)
	
	up.matrix <- NULL
	
	for(i in 1:length(up.locus)){
		x <- up.frame[up.frame[,1] == up.locus[i],2:ncol(up.frame)]
		x <- mean(x)
		up.matrix <- rbind(up.matrix,x)		
	}
	
	up.frame <- data.frame(up.locus,up.matrix)

	up.locus <- levels(as.factor(up.locus))
		
		up.frame.new <- NULL
		
		for(i in 1:length(up.locus)){
			
			datas <- up.frame[up.frame[,1] == up.locus[i],][,2:ncol(up.frame)]
			
			if(is.vector(datas) != TRUE){
				datas <- apply(datas,2,function(x) mean(as.double(x),na.rm=TRUE))
				datas[is.nan(datas)] <- NA
			}else{
				datas <- as.double(datas)
			}			
			
			up.frame.new <- rbind(up.frame.new,c(up.locus[i],datas))
			
		}
		
	up.frame <- as.data.frame(up.frame.new)
	
	

	down.matrix <- NULL
		
	for(i in 1:length(down.locus)){
		x <- down.frame[down.frame[,1] == down.locus[i],2:ncol(down.frame)]
		x <- mean(x)
		down.matrix <- rbind(down.matrix,x)		
	}
		
	down.frame <- data.frame(down.locus,down.matrix)
	
	down.locus <- levels(as.factor(down.locus))
		
		down.frame.new <- NULL
		
		for(i in 1:length(down.locus)){
			
			datas <- down.frame[down.frame[,1] == down.locus[i],][,2:ncol(down.frame)]
			
			if(is.vector(datas) != TRUE){
				datas <- apply(datas,2,function(x) mean(as.double(x),na.rm=TRUE))
				datas[is.nan(datas)] <- NA
			}else{
				datas <- as.double(datas)
			}
			
			down.frame.new <- rbind(down.frame.new,c(down.locus[i],datas))
			
		}
		
	down.frame <- as.data.frame(down.frame.new)
	

	return(list(up=up.locus,down=down.locus,up.frame=up.frame,down.frame=down.frame))
}


	rm()

}




#############################################################################################
#
# 4. Function LOCATION.ANALYSIS() -> Cytoband analysis of genes data
#
#############################################################################################

location.analysis <- function(list.anal,ref.list=NULL,nom,cyto,locus.name){

	cat(paste("\n\tLocation enrichment analysis of genes ",nom," started... ",format(Sys.time(), "%X"),sep=""))

	if(!is.null(ref.list)){	# filtering location data by reference list
		x <- NULL

		for(i in 1:length(ref.list)){
			if(length(cyto[cyto[,1] == as.vector(ref.list[i]),][1]) > 0){
				x <- rbind(x,cyto[cyto[,1] == as.vector(ref.list[i]),])
			}
		}

		cyto <- x
		rm(x)	
	}

	if(length(cyto[,1]) > 0){

		a <- as.factor(as.character(cyto[,1]))
		genes.cyto <- as.matrix(levels(a))	# list of all genes with cytoband information available
		genes.total <- length(genes.cyto[,1])	# number of all genes with cytoband information available
		rm(a)	



		a <- as.factor(as.character(cyto[,2]))
		chromo.id <- as.matrix(levels(a))	# list of IDs for the chromosomes
		rm(a)

		a <- as.factor(as.character(cyto[,3]))
		cyto.id <- as.matrix(levels(a))	# list of IDs for the cytobands
		rm(a)

		chromo.genes.nr <- matrix(0,length(chromo.id[,1]),1)	# list of number of genes located on each chromosome
		chromo.genes.list <- as.list(matrix(NA,length(chromo.id[,1]),1))	# list of vectors of genes located on each chromosome


		for(i in 1:length(chromo.id)){
			chromo.genes.list[[i]] <- cyto[,1][cyto[,2] == chromo.id[i]]
			chromo.genes.nr[i] <- length(chromo.genes.list[[i]])
		}

		cyto.genes.nr <- matrix(0,length(cyto.id[,1]),1)	# list of number of genes located on each cytoband
		cyto.genes.list <- as.list(matrix(NA,length(cyto.id[,1]),1))	# list of vectors of genes located on each cytoband


		for(i in 1:length(cyto.id)){
			cyto.genes.list[[i]] <- cyto[,1][cyto[,3] == cyto.id[i]]
			cyto.genes.nr[i] <- length(cyto.genes.list[[i]])
		}

		# returning chromosome and cytoband location data
		chromo.data <- list(annot.id=chromo.id, genes.nr=chromo.genes.nr, genes.total=genes.total, genes.list=chromo.genes.list, genes.annot=genes.cyto)
		cyto.data <- list(annot.id=cyto.id, genes.nr=cyto.genes.nr, genes.total=genes.total, genes.list=cyto.genes.list, genes.annot=genes.cyto)

		chromo.results <- annotate(annotations=chromo.data,exp.genes=list.anal,nom=nom,terms.name=cbind(chromo.id,chromo.id),taxoname="Chromosome")
		cyto.results <- annotate(annotations=cyto.data,exp.genes=list.anal,nom=nom,terms.name=cbind(cyto.id,cyto.id),taxoname="Cytoband")

		chromo.results <- pterm(exp.data=chromo.results,taxoname="Chromosome",nom=nom)
		cyto.results <- pterm(exp.data=cyto.results,taxoname="Cytoband",nom=nom)

		resterm(exp.data=chromo.results,taxoname="Chromosome",nom=nom,locus.name=locus.name)
		resterm(exp.data=cyto.results,taxoname="Cytoband",nom=nom,locus.name=locus.name)


	}else{
		cat(paste("\n\t\tCytoband enrichment analysis impossible...",sep=""))

	}

}





#############################################################################################
#
# 5. Function PRELIMINARY() -> Preliminary treatment of genes annotation data
#
#############################################################################################

preliminary <- function(file.annot,taxoname,restrict=FALSE,ref.list=NULL){

	cat(paste("\n\tPreliminary treatment of ",taxoname," annotations started... ",format(Sys.time(), "%X"),sep=""))

	if(restrict == TRUE && !is.null(ref.list)){

		x <- NULL

		for(i in 1:length(ref.list)){
			if(length(file.annot[file.annot[,1] == as.vector(ref.list[i]),][,1])){
				x <- rbind(x,file.annot[file.annot[,1] == as.vector(ref.list[i]),])
			}
		}

		file.annot <- x
		rm(x)
	}	

	if(length(file.annot[,1]) > 0){

		a <- as.factor(as.character(file.annot[,1]))
		genes.annot <- as.matrix(levels(a))	# list of all annotated genes within the considered taxonomical system
		genes.total <- length(genes.annot[,1])	# number of all annotated genes
		rm(a)	



		a <- as.factor(as.character(file.annot[,2]))
		annot.id <- as.matrix(levels(a))	# list of IDs for the terms annotating genes
		rm(a)

		genes.nr <- matrix(0,length(annot.id[,1]),1)	# list of number of genes annotated by each term
		genes.list <- as.list(matrix(NA,length(annot.id[,1]),1))	# list of vectors of genes annotated by each term


		for(i in 1:length(annot.id)){
			genes.list[[i]] <- file.annot[,1][file.annot[,2] == annot.id[i]]
			genes.nr[i] <- length(genes.list[[i]])
		}

	# returning annotations data
		annotations <- list(annot.id=annot.id, genes.nr=genes.nr, genes.total=genes.total, genes.list=genes.list, genes.annot=genes.annot)
	}else{
		annotations <- NULL
		cat(paste("\n\t\tNo ",taxoname," annotations are available...",sep=""))
	}


	return(annotations)
	rm()

}




#############################################################################################
#
# 6. Function ANNOTATE() -> Annotation of genes UP/DOWN
#
#############################################################################################

annotate <- function(annotations,exp.genes,nom,terms.name,taxoname){

	cat(paste("\n\t",taxoname," annotation of genes ", nom," started... ",format(Sys.time(), "%X"),sep=""))

	# getting annotation data
	annot.id <- annotations$annot.id
	genes.nr <- annotations$genes.nr
	genes.total <- annotations$genes.total
	genes.list <- annotations$genes.list
	genes.annot <- annotations$genes.annot	

	exp.genes <- as.matrix(as.character(as.matrix(exp.genes)))	# getting experience genes
	a <- matrix(0,length(exp.genes),1)	# calculating the number of annotated genes among experience genes

	for(i in 1:length(exp.genes)){
		a[i] <- length(genes.annot[,1][genes.annot[,1] == exp.genes[i]])
	}

	exp.total <- sum(a)	# number of experience genes annotated
	rm(a)


	exp.id <- as.matrix("")	# term IDs annotating experience genes
	exp.nr <- as.matrix(0)	# number of annotated genes by term (among experience ones)
	exp.list <- list(NULL)	# list of annotated genes by term (among experience ones)

	k <- 0

	for(i in 1:length(annot.id)){	# selecting a term ID

		# getting annotated genes list for the selected term
		a <- data.frame(as.character(genes.list[[i]]),matrix(0,length(genes.list[[i]]),1))

		# testing for annotated genes among experience ones
		for(j in 1:length(a[,1])){
			a[j,2] <- length(exp.genes[exp.genes == a[j,1]])		
		}

		if(sum(a[,2]) > 0){	# if experience genes were founded
			k + 1 -> k
			if(k == 1){	# for the first term
				exp.id[k] <- annot.id[i]
				exp.nr[k] <- sum(a[,2])	# number of annotated genes by the term (among experience ones)
				exp.list[[k]] <- as.vector(a[,1][a[,2] > 0])	# list of annotated genes by the term (among experience ones)
			}else{	# for the next terms
				exp.id <- rbind(exp.id,as.matrix(annot.id[i]))	# adding the term ID
				exp.nr <- rbind(exp.nr,as.matrix(sum(a[,2])))	# number of annotated genes by the term (among experience ones)
				exp.list[[length(exp.list)+1]]<- as.vector(a[,1][a[,2] > 0])	# list of annotated genes by the term (among experience ones)
			}
		}
	}

	rm(k,a)

	if(exp.id[1,1] != ""){	# if GO annotations were found...
		exp.term <- matrix("",length(exp.id),1)	# list of terms corresponding to IDs selected previously (annotating experience genes)

		for(i in 1:length(exp.id)){
			if(length(terms.name[,2][terms.name[,1] == exp.id[i]]) > 0){
				exp.term[i,1] <- as.character(terms.name[,2][terms.name[,1] == exp.id[i]])
			}else{
				exp.term[i,1] <- "---"	# avoiding taxonomy versions incompatibility
			}
		}

		exp.total <- matrix(exp.total,length(exp.id),1)

		x <- data.frame(as.character(annot.id),genes.nr)

		pop.hits <- matrix(0,length(exp.id),1)	# population genes annotated by each term
		pop.total <- matrix(genes.total,length(exp.id),1)	# total number of genes annotated among population genes

		for(i in 1:length(exp.id)){
			pop.hits[i] <- x[,2][x[,1] == exp.id[i]]
		}


		exp.data <- list(exp.id=exp.id, exp.term=exp.term, exp.nr=exp.nr, exp.total=exp.total, pop.hits=pop.hits, pop.total=pop.total, exp.list=exp.list)
	}else{	# if no GO annotations were founs...
		exp.data <- NULL
		cat(paste("\n\t\tNo annotations were detected for the genes ", nom, "...", sep=""))
	}


	return(exp.data)
	rm()
}




#############################################################################################
#
# 7. Function ANNOT.LIST() -> Annotate a list of genes for displaying purposes only (LocusLink)
#
#############################################################################################


annot.list <- function(ll,file.annot,taxoname,terms.name,locus.name){
	
	cat(paste("\n\t",taxoname," annotation of genes started...\n",sep=""))
	
	ll <- as.character(ll)
	annot.matrix <- matrix("",length(ll),2)
	annot.matrix <- cbind(ll,annot.matrix)
	rownames(locus.name) <- locus.name[,1]
	rownames(annot.matrix) <- ll
	
	for(i in 1:length(ll)){
		annot.gene <- file.annot[file.annot[,1] == ll[i],]
		
		if(length(annot.gene[,1]) > 0){
			annot.string <- ""
			
			for(j in 1:length(annot.gene[,2])){
				annot.string <- paste(annot.string,terms.name[terms.name[,1]==as.character(annot.gene[j,2]),2],"; ",sep="")
			}
			
			
			annot.matrix[i,3] <- annot.string
			
		}
		if(length(as.character(locus.name[locus.name[,1] == ll[i],2])) > 0){
			annot.matrix[i,2] <- as.character(locus.name[locus.name[,1] == ll[i],2])
		}
	
	}
	
	return(annot.matrix)

}




#############################################################################################
#
# 8. Function FDR.ADJUST() -> P-values adjustment using FDR of Benjamini & Hochberg (1995) modified by Paciorek (2004)
#
#############################################################################################



fdr.adjust <- function(pvalues, qlevel=0.05, method="original", adjust=NULL){	


	x.fdr <- fdr(pvals=pvalues,qlevel=qlevel,method=method,adjustment.method=adjust)	 
	y.fdr <- as.vector(matrix(1,length(pvalues),1))

	if(!is.null(x.fdr)){
		for(i in 1:length(x.fdr)){ y.fdr[x.fdr[i]]<-pvalues[x.fdr[i]]}
	}

	return(y.fdr)

}




#############################################################################################
#
# 9. Function PTERM() -> P-values calculation and adjustment
#
#############################################################################################

pterm <- function(exp.data,taxoname,nom){

	cat(paste("\n\t",taxoname," terms P-values calculation and adjustment for genes ",nom, " started... ",format(Sys.time(), "%X"),sep=""))
	if(!is.null(exp.data)){
	# getting experience data

		exp.id <- exp.data$exp.id
		exp.term <- exp.data$exp.term
		exp.nr <- exp.data$exp.nr
		exp.total <- exp.data$exp.total
		pop.hits <- exp.data$pop.hits
		pop.total <- exp.data$pop.total
		exp.list <- exp.data$exp.list


		term.p <- pvalues(data.frame(exp.nr,exp.total,pop.hits,pop.total))	# calculating p-values
		term.hommel <- p.adjust(term.p,method="hommel",n=length(term.p))	# Hommel adjusted p-values
		#term.fdr <- p.adjust(term.p,method="fdr",n=length(term.p))	# FDR (Benjamini & Hochberg) adjusted 
		term.fdr <- fdr.adjust(pvalues=term.p)	# FDR Paciorek


		exp.data <- list(exp.id=exp.id, exp.term=exp.term, exp.nr=exp.nr, exp.total=exp.total, pop.hits=pop.hits, pop.total=pop.total, exp.list=exp.list, term.p=term.p, term.hommel=term.hommel, term.fdr=term.fdr)

	}else{
		exp.data <-NULL
		cat("\n\t\tNo p-values were calculated...")
	}
	return(exp.data)
	rm()

}




#############################################################################################
#
# 10. Function GENE.CORREL() -> Grouping genes in clusters based on expression profiles correlation
#
#############################################################################################

gene.correl <- function(gene.frame,corr.th,method,nom){

	cat(paste("\n\tClustering genes ", nom," by correlating expression profiles started... ",format(Sys.time(), "%X"),sep=""))


	genes.correl <- as.list(NULL)

	w <- gene.frame


	if(method == "hierarchical"){

		x <- as.character(w[,1])

		y <- w[,2:ncol(w)]
		rownames(y) <- x

		y.corr <- rcorr(t(y),type="spearman")$r
		#cat(paste("\n\t\tCorrelation matrix computed... ",format(Sys.time(), "%X"),sep=""))
		y.dist <- as.dist(1 - y.corr)
		#hc <- agnes(y.dist,method = "ave")
		hc <- hclust(y.dist, method = "ave")
		sil.hc <- as.vector(NULL)

		for(i in 2:(nrow(y.corr)-1)){ 
    		# cut the tree 
    			memb <- cutree(hc, k = i)
    			sil <- silhouette(memb, y.dist)
    			sil.hc <- c(sil.hc, mean(summary(sil)$clus.avg.width))
  		}

		names(sil.hc) <- 2:(length(x)-1)

		best.index <- as.numeric(names(sil.hc[sil.hc == max(sil.hc)]))

		best.partition <- cutree(hc, k = best.index)

		for(i in 1:best.index){		
			genes.correl[[length(genes.correl) + 1]] <- as.character(names(best.partition[best.partition == i]))
		}
	}else if(method == "hgreedy"){

		x <- as.character(w[,1])

		y <- w[,2:ncol(w)]
		rownames(y) <- x

		y.corr <- rcorr(t(y),type="spearman")$r
		#cat(paste("\n\t\tCorrelation matrix computed... ",format(Sys.time(), "%X"),sep=""))

		genes.ok <- as.vector(NULL)
		for(i in 1:ncol(y.corr)){
			z <- y.corr[,i][names(y.corr[,i])!= x[i]]
			if(max(z)>= corr.th){genes.ok <- c(genes.ok,x[i])}
		}
		y.corr <- y.corr[rownames(y.corr) %in% genes.ok,colnames(y.corr) %in% genes.ok]

		y.dist <- as.dist(1 - y.corr)
		#hc <- agnes(y.dist,method = "ave")
		hc <- hclust(y.dist, method = "ave")
		sil.hc <- as.vector(NULL)

		for(i in 2:(nrow(y.corr)-1)){ 
    		# cut the tree 
    			memb <- cutree(hc, k = i)
    			sil <- silhouette(memb, y.dist)
    			sil.hc <- c(sil.hc, mean(summary(sil)$clus.avg.width))
  		}

		names(sil.hc) <- 2:(length(x)-1)

		best.index <- as.numeric(names(sil.hc[sil.hc == max(sil.hc)]))

		best.partition <- cutree(hc, k = best.index)

		for(i in 1:best.index){		
			genes.correl[[length(genes.correl) + 1]] <- as.character(names(best.partition[best.partition == i]))
		}

	}else if(method == "greedy"){

		w.init <- w
		x <- as.character(w[,1])

		y <- w[,2:ncol(w)]
		rownames(y) <- x

		y.corr <- rcorr(t(y),type="spearman")$r

		rownames(w) <- x
		rm(x)

		while(!is.null(w) && nrow(w) > 0){

		if(nrow(w)>1){
			x <- rownames(w)[1]
			y <- y.corr[,colnames(y.corr) == x]
		}else{
			x <- rownames(w)[1]
			y <- y.corr
		}


		if(length(y[y >= corr.th]) > 1){
			y <- y[y >= corr.th]
			z <- names(y)


			genes.correl[[length(genes.correl) + 1]] <- z
			w <- w[!(rownames(w) %in%  z),]
			y.corr <- y.corr[!(rownames(y.corr) %in%  z),!(colnames(y.corr) %in%  z)]




		}else{
			genes.correl[[length(genes.correl) + 1]] <- x


			if(nrow(w) > 1){
				w <- w[rownames(w) != x,]
				y.corr <- y.corr[!(rownames(y.corr) %in% x),!(colnames(y.corr) %in% x)]
			}else{
				w <- NULL
				y.corr <- NULL
			}

		}}


	}


	return(genes.correl)
	rm()

}





#############################################################################################
#
# 11. Function TERM.COMPARE() -> Comparing terms and clusters using a parameterised method
#
#############################################################################################

term.compare <- function(term1,term2,gene.frame,exp.total,compare,genes.correl){

	result <- NA	# the p-value to return

	if(compare == "common.genes"){

		names(term1) <- term1
		names(term2) <- term2

		z <- length(sort(term1[term2]))		

		if(z > 0){
			result <- pvalues(data.frame(z,length(term1),length(term2),exp.total))
		}

		rm(z)		

	}else if(compare == "correl.mean.exp"){

		x <- mean(gene.frame[term1,2:length(gene.frame)],na.rm=TRUE)
		y <- mean(gene.frame[term2,2:length(gene.frame)],na.rm=TRUE)

		if(suppressWarnings(cor.test(x,y,alternative="greater",method="spearman")[[4]] >= 0.8)){

			suppressWarnings(result <- cor.test(x,y,alternative="greater",method="spearman")[[3]])

		}

		rm(x,y)

	}else if(compare == "common.correl.genes"){

		n <- 0
		a <- length(term1)
		b <- length(term2)

		x <- term1
		names(x) <- x
		y <- term2
		names(y) <- y

		if(length(sort(x[y])) > 0){
			n <- length(sort(x[y]))
			z <- sort(x[y])
			for(i in 1:length(z)){
				term1 <- term1[term1 != z[i]]
				term2 <- term2[term2 != z[i]]
			}
		}

		rm(x,y)	#


		if(length(genes.correl) > 0){	# avoid incorrelable genes with 0 gene clusters
			for(i in 1:length(genes.correl)){
				z <- genes.correl[[i]]
				names(z) <- z
				x <- sort(z[term1])
				y <- sort(z[term2])

				if(length(x) > 0 & length(y) > 0){
					m <- min(length(x),length(y))
					n <- n + m
				}			
			}
		}

		if(n > 0){
			result <- pvalues(data.frame(n,a,b,exp.total))
		}

	}

	return(result)
	rm()
}




#############################################################################################
#
# 12. Function CLUSTERM.CC() -> Grouping terms in clusters based on common annotations in a concomitant way
#
#############################################################################################

clusterm.cc <- function(exp.data,taxoname,nom,alpha,compare,gene.frame,genes.correl){

	cat(paste("\n\t",taxoname," annotation concomitent clustering of genes ", nom," started... ",format(Sys.time(), "%X"),sep=""))

if(!is.null(exp.data)){	# 

	v <- as.list(NULL)

	for(i in 1:length(genes.correl)){
		if(length(genes.correl[[i]]) > 1){
			v[[length(v) + 1]] <- genes.correl[[i]]
		}
	}

	genes.correl <- v
	rm(v)

	exp.id <- exp.data$exp.id
	names(exp.id) <- exp.id
	exp.term <- exp.data$exp.term
	names(exp.term) <- exp.id
	term.p <- exp.data$term.p
	names(term.p) <- exp.id
	exp.list <- exp.data$exp.list
	names(exp.list) <- exp.id


	exp.id <- exp.id[order(term.p)]
	exp.term <- exp.term[order(term.p)]
	exp.list <- exp.list[order(term.p)]
	term.p <- sort(term.p)

	exp.total <- exp.data$exp.total[1]

# p-values term selection

	if(alpha != 1){
		term.p <- term.p[term.p <= alpha]
		exp.id <- exp.id[names(term.p)]
		exp.term <- exp.term[names(term.p)]
		exp.list <- exp.list[names(term.p)]		
	}

# starting term clustering

	id.cluster <- as.list(NULL)
	gene.cluster <- as.list(NULL)


	w <- exp.id	# the list of terms not already clustered (initially all the terms)


	while(length(w) > 0){
		x <- matrix("",length(w),1)
		names(x) <- w


		for(i in 1:length(w)){

		# selecting the first term among those not already clustered		
			a <- as.matrix(exp.list[w[i]][[1]])

		# comparing the significance of association with the rest of not already clustered terms
			y <- matrix(NA,length(w),1)			
			names(y) <- w

			for(j in 1:length(w)){
				if(w[j] != w[i]){
					b <- as.matrix(exp.list[w[j]][[1]])

		# using a parameterised method for comparing terms		
					y[j] <- term.compare(a,b,gene.frame,exp.total,compare,genes.correl)

					rm(b)
				}		
			}

			if(length(sort(y)) > 0){

				y <- sort(y)
				names.y <- names(y)
				#y <- p.adjust(y,method="fdr",n=length(y))
				y <- fdr.adjust(pvalues=y)	# FDR Paciorek

				names(y) <- names.y			
				y <- y[y <= 0.05]
			}else{
				y <- NULL
			}

		# comparing the significance of association with existing clusters
			u <- NULL

			if(length(id.cluster) > 0){

				u <- matrix(2,length(id.cluster),1)
				names(u) <- as.character(1:length(u))

				for(j in 1:length(id.cluster)){
					b <- as.matrix(gene.cluster[[j]])

		# using a parameterised method for comparing terms and clusters
					u[j] <- term.compare(a,b,gene.frame,exp.total,compare,genes.correl)

					rm(b)
				}				
			}


			if(length(sort(u)) > 0){

				u <- sort(u)
				names.u <- names(u)
				#u <- p.adjust(u,method="fdr",n=length(u))
				u <- fdr.adjust(pvalues=u)	# FDR Paciorek

				names(u) <- names.u
				u <- u[u <= 0.05]
			}else{
				u <- NULL
			}

			yu <- c(y,u)


			if(length(yu) > 0){
				yu <- yu[yu == sort(yu)[1]]

				if(yu[1] <= 0.05){

					if(length(yu) > 1){
						cluster.found <- FALSE

						for(j in 1:length(yu)){
							z <- as.character(1:length(id.cluster))
							if(length(z[z == names(yu[j])]) > 0 & cluster.found == FALSE){
								x[i] <- names(yu[j])
								cluster.found <- TRUE
							}
						}

						if(cluster.found == FALSE){
							x[i] <- names(yu[1])
						}
					}else{
						x[i] <- names(yu[1])
					}					

				}
			}

			rm(a,y,u,yu)		
		}

		if(length(x[x == ""]) == length(x)){
			n <- length(id.cluster)

			for(i in 1:length(x)){
				n <- n + 1
				id.cluster[[n]] <- names(x[i])
				gene.cluster[[n]] <- exp.list[[names(x[i])]]
			}
			w <- NULL

		}else{

			x <- x[x != ""]

			if(length(id.cluster) > 0){

				modif.cluster <- FALSE

				for(i in 1:length(id.cluster)){
					if(length(x[x == as.character(i)]) > 0){
						modif.cluster <- TRUE
						y <- x[x == as.character(i)]
						for(j in 1:length(y)){
							id.cluster[[i]] <- c(id.cluster[[i]],names(y[j]))
							gene.cluster[[i]] <- levels(as.factor(c(gene.cluster[[i]],exp.list[[names(y[j])]])))
							w <- w[w != names(y[j])]
						}
					}
				}

				if(modif.cluster == FALSE){
					id.cluster[[length(id.cluster) + 1]] <- c(x[[1]],names(x[1]))
					gene.cluster[[length(gene.cluster) + 1]] <- levels(as.factor(c(exp.list[[names(x[1])]],exp.list[[x[[1]]]])))
					w <- w[w != x[[1]] & w != names(x[1])]
				}

			}else{
				id.cluster[[length(id.cluster) + 1]] <- c(x[[1]],names(x[1]))
				gene.cluster[[length(gene.cluster) + 1]] <- levels(as.factor(c(exp.list[[names(x[1])]],exp.list[[x[[1]]]])))
				w <- w[w != x[[1]] & w != names(x[1])]
			}

			rm(x)		

		}


	}


	if(length(id.cluster) > 0){

		term.cluster <- as.list(NULL)

		for(i in 1:length(id.cluster)){
			a <- as.vector(NULL)
			for(j in 1:length(id.cluster[[i]])){
				a <- c(a,exp.term[as.character(id.cluster[[i]][j])])
			}
			term.cluster[[i]] <- a
		}

		clusters <- list(id.cluster=id.cluster,term.cluster=term.cluster,gene.cluster=gene.cluster)
	}else{
		clusters <- NULL
		cat("\n\t\tNo clusters were created...")
	}
}else{
	clusters <- NULL
	cat("\n\t\tNo clusters were created...")
}
	return(clusters)
	rm()

}



#############################################################################################
#
# 13. Function CLUSTERM.CS() -> Grouping terms in clusters based on common annotations in a consecutive way
#
#############################################################################################

clusterm.cs <- function(exp.data,taxoname,nom,alpha,compare,gene.frame,genes.correl){

	cat(paste("\n\t",taxoname," annotation consecutive clustering of genes ", nom," started... ",format(Sys.time(), "%X"),sep=""))

if(!is.null(exp.data)){

	v <- as.list(NULL)

	for(i in 1:length(genes.correl)){
		if(length(genes.correl[[i]]) > 1){
			v[[length(v) + 1]] <- genes.correl[[i]]
		}
	}

	genes.correl <- v
	rm(v)

	exp.id <- exp.data$exp.id
	names(exp.id) <- exp.id
	exp.term <- exp.data$exp.term
	names(exp.term) <- exp.id
	term.p <- exp.data$term.p
	names(term.p) <- exp.id
	exp.list <- exp.data$exp.list
	names(exp.list) <- exp.id


	exp.id <- exp.id[order(term.p)]
	exp.term <- exp.term[order(term.p)]
	exp.list <- exp.list[order(term.p)]
	term.p <- sort(term.p)

	exp.total <- exp.data$exp.total[1]

# p-values term selection

	if(alpha != 1){
		term.p <- term.p[term.p <= alpha]
		exp.id <- exp.id[names(term.p)]
		exp.term <- exp.term[names(term.p)]
		exp.list <- exp.list[names(term.p)]		
	}

# starting term clustering

	id.cluster <- as.list(NULL)

	id.cluster[[1]] <- exp.id[1]
	gene.cluster <- as.list(NULL)
	gene.cluster[[1]] <- exp.list[as.character(exp.id[1])][[1]]


	x <- exp.id[2:length(exp.id)]
	names(x) <- x
	y <- matrix(NA,length(x),1)
	names(y) <- x

	k <- 1

	while(length(x) > 0){


		for(i in 1:length(x)){
			a <- x[i]
			b <- as.matrix(exp.list[a][[1]])

		# using a parameterised method for comparing terms and clusters
			y[i] <- term.compare(b,gene.cluster[[k]],gene.frame,exp.total,compare,genes.correl)

		}

		if(length(sort(y)) > 0){
			y <- sort(y)			
			z <- data.frame(x[names(y)],y)
			z <- z[order(z[,2]),]
			#z[,2] <- p.adjust(z[,2],method="fdr",n=length(z[,2]))
			z[,2] <- fdr.adjust(pvalues=z[,2])	# FDR Paciorek


			if(z[1,2] <= 0.05){
				id.cluster[[k]] <- c(id.cluster[[k]],as.character(z[1,1]))
				gene.cluster[[k]] <- levels(as.factor(c(gene.cluster[[k]],exp.list[as.character(z[1,1])][[1]])))
				x <- x[x != as.character(z[1,1])]
				names(x) <- x
				y <- matrix(NA,length(x),1)
				names(y) <- x
			}else if(length(x) > 1){
				k <- k + 1
				id.cluster[[k]] <- x[1]
				gene.cluster[[k]] <- exp.list[as.character(x[1])][[1]]
				x <- x[2:length(x)]
				names(x) <- x
				y <- matrix(NA,length(x),1)
				names(y) <- x
			}else{
				k <- k + 1
				id.cluster[[k]] <- x[1]
				gene.cluster[[k]] <- exp.list[as.character(x[1])][[1]]
				x <- as.vector(NULL)
			}

		}else if(length(x) >= 2){
			k <- k + 1
			id.cluster[[k]] <- x[1]
			gene.cluster[[k]] <- exp.list[as.character(x[1])][[1]]
			x <- x[2:length(x)]
			names(x) <- x
			y <- matrix(NA,length(x),1)
			names(y) <- x

		}else{
			k <- k + 1
			id.cluster[[k]] <- x[1]
			gene.cluster[[k]] <- exp.list[as.character(x[1])][[1]]
			x <- as.vector(NULL)
		}

	}

	if(length(id.cluster) > 0){


		term.cluster <- as.list(NULL)

		for(i in 1:length(id.cluster)){
			a <- as.vector(NULL)
			for(j in 1:length(id.cluster[[i]])){
				a <- c(a,exp.term[as.character(id.cluster[[i]][j])])
			}
			term.cluster[[i]] <- a
		}

		clusters <- list(id.cluster=id.cluster,term.cluster=term.cluster,gene.cluster=gene.cluster)
	}else{
		clusters <- NULL
		cat("\n\t\tNo clusters were created...")
	}

}else{
	clusters <- NULL
	cat("\n\t\tNo clusters were created...")
}
	return(clusters)	
	rm()	
}



#############################################################################################
#
# 14. Function CLUSTER.C() -> Weight center calculation for the clusters of terms
#
#############################################################################################

cluster.c <- function(clusters,exp.data,annotations,taxoname,nom){

	cat(paste("\n\t",taxoname," clusters weight center calculation for genes ",nom, " started... ",format(Sys.time(), "%X"),sep=""))

if(!is.null(clusters)){

	id.cluster <- clusters$id.cluster
	term.cluster <- clusters$term.cluster
	gene.cluster <- clusters$gene.cluster

	annot.id <- annotations$annot.id
	genes.nr <- annotations$genes.nr
	names(genes.nr) <- annot.id


	exp.id <- exp.data$exp.id
	names(exp.id) <- exp.id
	exp.term <- exp.data$exp.term
	names(exp.term) <- exp.id
	exp.list <- exp.data$exp.list
	names(exp.list) <- exp.id

	center.cluster <- matrix("NA",length(id.cluster),1)

	for(i in 1:length(id.cluster)){
		x <- id.cluster[[i]]
		y <- gene.cluster[[i]]
		names(y) <- y

		for(j in 1:length(x)){
			z <- exp.list[[x[j]]]
			y <- sort(y[z])		
		}

		if(length(y) > 0){

			z <- matrix(2,length(x),1)
			names(z) <- x

			for(j in 1:length(x)){
				z[j] <- pvalues(data.frame(length(y),length(y),length(exp.list[[x[j]]]),length(gene.cluster[[i]])))
			}


			if(length(z[z != 2]) > 0){
				names.z <- names(z[z != 2])
				z <- z[z != 2]
				#z <- p.adjust(z,method="fdr",n=length(z))
				z <- fdr.adjust(pvalues=z)	# FDR Paciorek

				names(z) <- names.z

				z <- sort(z)

				if(length(z[z == z[1]]) > 1){
					z <- z[z == z[1]]
					w <- matrix(0,length(z),1)
					names(w) <- names(z)

					for(j in 1:length(w)){
						w[j] <- genes.nr[[names(w[j])]]
					}

					center.cluster[i] <- names(sort(w)[1])

				}else{
					center.cluster[i] <- names(z[1])
				}								
			}			
		}
		rm(x,y,z)
	}


	clusters <- list(id.cluster=id.cluster,term.cluster=term.cluster,gene.cluster=gene.cluster,center.cluster=center.cluster)
}else{
	clusters <- NULL
	cat("\n\t\tNo weight center was calculated...")
}


	return(clusters)	
	rm()

}




#############################################################################################
#
# 15. Function PCLUST() -> P-values calculation and adjustment for the clusters of terms
#
#############################################################################################


pclust <- function(exp.data,clusters,annotations,taxoname,nom){

	cat(paste("\n\t",taxoname," clusters P-values calculation and adjustment for genes ",nom, " started... ",format(Sys.time(), "%X"),sep=""))

if(!is.null(clusters)){

	id.cluster <- clusters$id.cluster
	term.cluster <- clusters$term.cluster
	gene.cluster <- clusters$gene.cluster
	center.cluster <- clusters$center.cluster

	exp.id <- exp.data$exp.id
	names(exp.id) <- exp.id
	exp.list <- exp.data$exp.list
	names(exp.list) <- exp.id
	exp.total <- exp.data$exp.total[1]

	annot.id <- annotations$annot.id 
	genes.total <- annotations$genes.total
	genes.list <- annotations$genes.list
	names(genes.list) <- annot.id

	list.hits <- matrix(0,length(id.cluster),1)
	list.total <- matrix(exp.total,length(id.cluster),1)
	pop.hits <- matrix(0,length(id.cluster),1)
	pop.total <- matrix(genes.total,length(id.cluster),1)

	for(i in 1:length(id.cluster)){
		list.hits[i] <- length(gene.cluster[[i]])
		a <- genes.list[id.cluster[[i]]]
		b <- NULL

		for(j in 1:length(a)){
			b <- c(b,a[[j]])
		}

		b <- levels(as.factor(b))
		pop.hits[i] <- length(b)
		rm(a,b)		
	}

	data.cluster <- data.frame(list.hits,list.total,pop.hits,pop.total)

	cluster.p <- pvalues(data.cluster)
	cluster.hommel <- p.adjust(cluster.p,method="hommel",n=length(cluster.p))
	#cluster.fdr <- p.adjust(cluster.p,method="fdr",n=length(cluster.p))
	cluster.fdr <- fdr.adjust(pvalues=cluster.p)	# FDR Paciorek

	id.cluster <- id.cluster[order(cluster.p)]
	term.cluster <- term.cluster[order(cluster.p)]
	gene.cluster <- gene.cluster[order(cluster.p)]
	center.cluster <- center.cluster[order(cluster.p)]

	data.cluster <- data.cluster[order(cluster.p),]

	cluster.hommel <- cluster.hommel[order(cluster.p)]
	cluster.fdr <- cluster.fdr[order(cluster.p)]
	cluster.p <- sort(cluster.p)

	data.cluster <- data.frame(data.cluster,cluster.p,cluster.hommel,cluster.fdr)

	clusters <- list(id.cluster=id.cluster,term.cluster=term.cluster,gene.cluster=gene.cluster,center.cluster=center.cluster,data.cluster=data.cluster)
}else{
	clusters <- NULL
	cat("\n\t\tNo p-values were calculated...")
}


	return(clusters)
	rm()

}



#############################################################################################
#
# 16. Function RESCLUST() -> Organizing and saving of the results (clusters of terms and genes)
#
#############################################################################################

resclust <- function(exp.data,clusters,taxoname,nom,pref,locus.name,results.dir){

	cat(paste("\n\tSaving ",taxoname," ",nom," clustering results... ",format(Sys.time(), "%X"),sep=""))

if(!is.null(clusters)){

	id.cluster <- clusters$id.cluster
	term.cluster <- clusters$term.cluster
	gene.cluster <- clusters$gene.cluster
	center.cluster <- clusters$center.cluster
	data.cluster <- clusters$data.cluster

	term.p <- exp.data$term.p
	term.hommel <- exp.data$term.hommel[order(term.p)]
	term.fdr <- exp.data$term.fdr[order(term.p)]

	exp.id <- exp.data$exp.id[order(term.p)]
	exp.term <- exp.data$exp.term[order(term.p)]
	exp.nr <- exp.data$exp.nr[order(term.p)]

	exp.total <- exp.data$exp.total[order(term.p)]
	pop.hits <- exp.data$pop.hits[order(term.p)]
	pop.total <- exp.data$pop.total[order(term.p)]
	exp.list <- exp.data$exp.list[order(term.p)]

	term.p <- sort(term.p)
	index <- row(as.matrix(term.p))


	names(exp.term) <- exp.id
	names(exp.nr) <- exp.id

	names(exp.total) <- exp.id
	names(pop.hits) <- exp.id

	names(pop.total) <- exp.id
	names(exp.list) <- exp.id
	names(term.p) <- exp.id

	names(term.hommel) <- exp.id
	names(term.fdr) <- exp.id
	names(index) <- exp.id


	str.genes.cluster <- matrix("",length(gene.cluster),1)
	str.names.cluster <- matrix("",length(gene.cluster),1)

	for(i in 1:length(gene.cluster)){
		gene.cluster[[i]] -> a

		"" -> b
		"" -> d
		for(j in 1:length(a)){
			paste(b,a[j],"\n",sep="") -> b
			paste(d,"<option value='",locus.name[locus.name[,1] == a[j],2],"'>",locus.name[locus.name[,1] == a[j],2],"</option>",sep="") -> d
		}
		b -> str.genes.cluster[i]
		d -> str.names.cluster[i]
	}

	rm(a,b,d)

	wd <- getwd()	
	options(warn=-1)

	head <- paste("<html>

<head>
<meta http-equiv='Content-Type' content='text/html; charset=windows-1252'>
<title>Genes ",nom," - ",taxoname," Annotation Clusters</title>
</head>

<body link='#0000FF' vlink='#0000FF' alink='#0000FF'>

<table id='table0' style='width: 700px; border-collapse: collapse' cellPadding='0' border='0'>
	<tr>
		<td style='font-weight: 700; font-size: 10pt; vertical-align: middle; color: navy; font-style: normal; font-family: Arial, sans-serif; white-space: nowrap; height: 38pt; text-decoration: none; text-align: general; border: medium none; padding-left: 1px; padding-right: 1px; padding-top: 1px' align='justify'>
		<span lang='en-us' style='LETTER-SPACING: 2pt'><font size='3'>Genes ",nom," -
		",taxoname," Annotation Clusters</font></span></td>
	</tr>
	<tr>
		<td style='font-weight: 400; font-size: 10pt; vertical-align: middle; width: 645pt; color: windowtext; font-style: normal; font-family: Arial; white-space: nowrap; text-decoration: none; text-align: general; border: medium none; padding-left: 1px; padding-right: 1px; padding-top: 1px' align='justify'>
		&nbsp;</td>
	</tr>
	<tr>
		<td style='font-weight: 400; font-size: 10pt; vertical-align: middle; width: 645pt; color: windowtext; font-style: normal; font-family: Arial; white-space: nowrap; text-decoration: none; text-align: general; border: medium none; padding-left: 1px; padding-right: 1px; padding-top: 1px' align='justify'>
		<font face='Times New Roman' size='3'>
		<b><span lang='en-us'>Used notations:</span></b></font></td>
	</tr>
	<tr>
		<td style='font-weight: 400; font-size: 10pt; vertical-align: middle; width: 645pt; color: windowtext; font-style: normal; font-family: Arial; white-space: nowrap; text-decoration: none; text-align: general; border: medium none; padding-left: 1px; padding-right: 1px; padding-top: 1px' align='justify'>
		&nbsp;</td>
	</tr>
	<tr>
		<td style='font-weight: 400; font-size: 10pt; vertical-align: middle; width: 645pt; color: windowtext; font-style: normal; font-family: Arial; white-space: nowrap; text-decoration: none; text-align: general; border: medium none; padding-left: 1px; padding-right: 1px; padding-top: 1px' align='justify'>
		<ul>
			<li><font face='Times New Roman' size='3'><span lang='en-us'><b>List Hits</b> - the number of genes
			annotated by the considered ",taxoname," category or annotation cluster within
			the analyzed list of target genes</span> </font> </li>
		</ul>
		</td>
	</tr>
	<tr>
		<td style='font-weight: 400; font-size: 10pt; vertical-align: middle; width: 645pt; color: windowtext; font-style: normal; font-family: Arial; white-space: nowrap; text-decoration: none; text-align: general; border: medium none; padding-left: 1px; padding-right: 1px; padding-top: 1px' align='justify'>
		<ul>
			<li><font face='Times New Roman' size='3'><span lang='en-us'><b>List </b></span><b>Total</b><span lang='en-us'>
			- the number of genes within the analyzed list of target genes
			having at least one ",taxoname," annotation</span> </font> </li>
		</ul>
		</td>
	</tr>
	<tr>
		<td style='font-weight: 400; font-size: 10pt; vertical-align: middle; width: 645pt; color: windowtext; font-style: normal; font-family: Arial; white-space: nowrap; text-decoration: none; text-align: general; border: medium none; padding-left: 1px; padding-right: 1px; padding-top: 1px' align='justify'>
		<ul>
			<li><font face='Times New Roman' size='3'><b>Population<span lang='en-us'> Hits</span></b><span lang='en-us'>
			- the number of genes, available on the entire microarray, annotated
			by the considered ",taxoname," category or annotation cluster</span>
			</font>
			</li>
		</ul>
		</td>
	</tr>
	<tr>
		<td style='font-weight: 400; font-size: 10pt; vertical-align: middle; width: 645pt; color: windowtext; font-style: normal; font-family: Arial; white-space: nowrap; text-decoration: none; text-align: general; border: medium none; padding-left: 1px; padding-right: 1px; padding-top: 1px' align='justify'>
		<ul>
			<li><font face='Times New Roman' size='3'><b>Population Total</b> - the number of genes available on the
			entire microarray and having at least one ",taxoname," annotation </font>
			</li>
		</ul>
		</td>
	</tr>
	<tr>
		<td style='font-weight: 400; font-size: 10pt; vertical-align: middle; width: 645pt; color: windowtext; font-style: normal; font-family: Arial; white-space: nowrap; text-decoration: none; text-align: general; border: medium none; padding-left: 1px; padding-right: 1px; padding-top: 1px' align='justify'>
		<ul>
			<li><font face='Times New Roman' size='3'><span lang='en-us'><b>P-value</b> - the significance p-value of
			the gene enrichment of the considered ",taxoname," category or annotation
			cluster, calculated with a unilateral Fisher exact test</span>
			</font> </li>
		</ul>
		</td>
	</tr>
	<tr>
		<td style='font-weight: 400; font-size: 10pt; vertical-align: middle; width: 645pt; color: windowtext; font-style: normal; font-family: Arial; white-space: nowrap; text-decoration: none; text-align: general; border: medium none; padding-left: 1px; padding-right: 1px; padding-top: 1px' align='justify'>
		<ul>
			<li><font face='Times New Roman' size='3'><span lang='en-us'><b>FDR (q-value)</b> - the estimated value of
			the false discovery rate (FDR) by using the Benjamini &amp; Hochberg
			(2001) method</span> </font> </li>
		</ul>
		</td>
	</tr>
</table>
<p>",sep="")


	write(head, file = paste(wd,"/",results.dir,"/Genes ",nom," - ",taxoname," Annotation Clusters.htm",sep=""), append = TRUE, sep = " ")


	for(i in 1:length(id.cluster)){

		cluster.head <- paste("
		<p>
..........................................................................................................................</p>
<p><b><span style='LETTER-SPACING: 2pt; color: navy; font-family: Arial, sans-serif'><font size='3'>Cluster n ",i,"</font></span></b></p>
<table border='1' id='table",i,"' style='border-collapse: collapse' bordercolor='#808080'>
	<tr>
		<td align='center' width='45'><b>Terms</b></td>
		<td align='center' width='45'>
		<p style='margin-top: 0; margin-bottom: 0'><b>List </b></p>
		<p style='margin-top: 0; margin-bottom: 0'><b>Hits</b></td>
		<td align='center' width='45'>
		<p style='margin-top: 0; margin-bottom: 0'><b>List </b></p>
		<p style='margin-top: 0; margin-bottom: 0'><b>Total</b></td>
		<td align='center' width='70'>
		<p style='margin-top: 0; margin-bottom: 0'><b>Population </b></p>
		<p style='margin-top: 0; margin-bottom: 0'><b>Hits</b></td>
		<td align='center' width='70'>
		<p style='margin-top: 0; margin-bottom: 0'><b>Population </b></p>
		<p style='margin-top: 0; margin-bottom: 0'><b>Total</b></td>
		<td align='center' width='100'>
		<p style='margin-top: 0; margin-bottom: 0'><b>FDR </b></p>
		<p style='margin-top: 0; margin-bottom: 0'><b>(q-value)</b></td>
		<td align='center' width='100'><b>Genes ID's</b></td>
		<td align='center' width='200'><b>Genes Names</b></td>
	</tr>
	<tr>
		<td align='center' width='45'><font size='3'>",length(id.cluster[[i]]),"</font></td>
		<td align='center' width='45'><font size='3'>",data.cluster[i,1],"</font></td>
		<td align='center' width='45'><font size='3'>",data.cluster[i,2],"</font></td>
		<td align='center' width='70'><font size='3'>",data.cluster[i,3],"</font></td>
		<td align='center' width='70'><font size='3'>",data.cluster[i,4],"</font></td>
		<td align='center' width='70'><font size='3'>",format(data.cluster[i,7],digits=3,scientific=TRUE),"</font></td>
		<td align='center' width='100'><textarea rows='0' name='Genes",i,"1' cols='7'>",str.genes.cluster[i],"</textarea></td>
		<td align='center' width='200'><select size='1' name='Genes",i,"2'>
	",str.names.cluster[i],"
	</select></td>
	</tr>
</table>
<p><b>Terms</b></p>
<table border='1' id='table2' style='border-collapse: collapse' bordercolor='#808080'>
	<tr>
		<td align='center' width='70'><b>ID</b></td>
		<td align='center' width='200'><b>Term</b></td>
		<td align='center' width='45'><b>Rank</b></td>
		<td align='center' width='45'>
		<p style='margin-top: 0; margin-bottom: 0'><b>List </b></p>
		<p style='margin-top: 0; margin-bottom: 0'><b>Hits</b></td>
		<td align='center' width='45'>
		<p style='margin-top: 0; margin-bottom: 0'><b>List </b></p>
		<p style='margin-top: 0; margin-bottom: 0'><b>Total</b></td>
		<td align='center' width='70'>
		<p style='margin-top: 0; margin-bottom: 0'><b>Population </b></p>
		<p style='margin-top: 0; margin-bottom: 0'><b>Hits</b></td>
		<td align='center' width='70'>
		<p style='margin-top: 0; margin-bottom: 0'><b>Population </b></p>
		<p style='margin-top: 0; margin-bottom: 0'><b>Total</b></td>
		<td align='center' width='100'><b>P-value</b></td>
		<td align='center' width='100'><b>Genes ID's</b></td>
		<td align='center' width='200'><b>Genes Names</b></td>
	</tr>",sep="")

		write(cluster.head, file = paste(wd,"/",results.dir,"/Genes ",nom," - ",taxoname," Annotation Clusters.htm",sep=""), append = TRUE, sep = " ")


		for(j in 1:length(id.cluster[[i]])){

			exp.term.lnk <- ""
			exp.id.lnk <- ""
			str.genes.term <- ""
			str.names.term <- ""


			exp.list[[id.cluster[[i]][j]]] -> a

			"" -> b
			"" -> d
			for(k in 1:length(a)){
				paste(b,a[k],"\n",sep="") -> b
				paste(d,"<option value='",locus.name[locus.name[,1] == a[k],2],"'>",locus.name[locus.name[,1] == a[k],2],"</option>",sep="") -> d
			}
			b -> str.genes.term
			d -> str.names.term

			rm(a,b)			

			if(taxoname == "GO Biological Process" || taxoname == "GO Cellular Component" || taxoname == "GO Molecular Function"){
				exp.term.lnk <- paste("<a target='_blank' href='http://www.ebi.ac.uk/ego/DisplayGoTerm?id=",id.cluster[[i]][j],"' style='text-decoration: none'>
						<font size='3'>",term.cluster[[i]][j],"</font></a>",sep="")
				exp.id.lnk <- paste("<a target='_blank' href='http://www.ebi.ac.uk/ego/DisplayGoTerm?id=",id.cluster[[i]][j],"' style='text-decoration: none'>
						<font size='3'>",id.cluster[[i]][j],"</font></a>",sep="")
			}else{
				exp.term.lnk <- paste("<font size='3'>",term.cluster[[i]][j],"</font>",sep="")
				exp.id.lnk <- paste("<font size='3'>",id.cluster[[i]][j],"</font>",sep="")
			}




		cluster.terms <- paste("
		<tr>
		<td align='center' width='70'>
		",exp.id.lnk,"</td>
		<td align='center' width='200'>
		",exp.term.lnk,"</td>
		<td align='center' width='45'><font size='3'>",index[id.cluster[[i]][j]],"</font></td>
		<td align='center' width='45'><font size='3'>",exp.nr[id.cluster[[i]][j]],"</font></td>
		<td align='center' width='45'><font size='3'>",data.cluster[i,2],"</font></td>
		<td align='center' width='70'><font size='3'>",pop.hits[id.cluster[[i]][j]],"</font></td>
		<td align='center' width='70'><font size='3'>",data.cluster[i,4],"</font></td>
		<td align='center' width='70'><font size='3'>",format(term.p[id.cluster[[i]][j]], digits = 3, scientific = TRUE),"</font></td>
		<td align='center' width='100'><textarea rows='0' name='Genes",i,"3' cols='7'>",str.genes.term,"</textarea></td>
		<td align='center' width='200'><select size='1' name='Genes",i,"4'>",str.names.term,"</select></td>
		</tr>",sep="")


		write(cluster.terms, file = paste(wd,"/",results.dir,"/Genes ",nom," - ",taxoname," Annotation Clusters.htm",sep=""), append = TRUE, sep = " ")

		}
		write("</table>", file = paste(wd,"/",results.dir,"/Genes ",nom," - ",taxoname," Annotation Clusters.htm",sep=""), append = TRUE, sep = " ")

	}

		end.page <- paste("<p>
		..........................................................................................................................</p>
		<p>
		This file was produced on ",date()," with <a target='_blank' href='http://corneliu.henegar.info/FunCluster.htm' style='text-decoration: none'>
		<font size='3'>FunCluster ",f.version,"</font></a> using GO and KEGG annotations updated on ",annot.date,".</p>
		</body>

		</html>",sep="")

		write(end.page, file = paste(wd,"/",results.dir,"/Genes ",nom," - ",taxoname," Annotation Clusters.htm",sep=""), append = TRUE, sep = " ")



}else{
	cat(paste("\n\t\tNo ",taxoname," ",nom," clustering results to save...",sep=""))
}

	rm()

}




#############################################################################################
#
# 17. Function RESTERM() -> Organizing and saving of the results (isolated terms)
#
#############################################################################################

resterm <- function(exp.data,taxoname,nom,locus.name,results.dir){

	cat(paste("\n\tSaving ",taxoname," ",nom," annotation results... ",format(Sys.time(), "%X"),sep=""))

if(!is.null(exp.data)){

	term.p <- exp.data$term.p
	exp.id <- exp.data$exp.id[order(term.p)]
	exp.term <- exp.data$exp.term[order(term.p)]
	exp.nr <- exp.data$exp.nr[order(term.p)]
	exp.total <- exp.data$exp.total[order(term.p)]
	pop.hits <- exp.data$pop.hits[order(term.p)]
	pop.total <- exp.data$pop.total[order(term.p)]
	exp.list <- exp.data$exp.list[order(term.p)]

	term.p <- sort(term.p)
	term.hommel <- exp.data$term.hommel[order(term.p)]
	term.fdr <- exp.data$term.fdr[order(term.p)]
	#term.qval <- exp.data$term.qval



# list (as a string) of genes for each term

	list.genes <- matrix("",length(exp.list),1)
	list.names <- matrix("",length(exp.list),1)

	for(i in 1:length(exp.list)){
		exp.list[[i]] -> a

		"" -> b
		"" -> d
		for(j in 1:length(a)){
			paste(b,a[j],'\n',sep="") -> b
			paste(d,"<option value='",locus.name[locus.name[,1] == a[j],2],"'>",locus.name[locus.name[,1] == a[j],2],"</option>",sep="") -> d

		}
		b -> list.genes[i]
		d -> list.names[i]

	}

	rm(a,b,d)


# routine for saving the results of treatment for GO annotations as HTML files



	wd <- getwd()
	exp.term.lnk <- ""
	exp.id.lnk <- ""


	head <- paste("<html>

	<head>
	<meta http-equiv='Content-Type' content='text/html; charset=windows-1252'>
	<title>Genes ",nom," - ",taxoname," Annotations</title>
	</head>

	<body link='#0000FF' vlink='#0000FF' alink='#0000FF'>

	<table id='table0' style='width: 700px; border-collapse: collapse' cellPadding='0' border='0'>
		<tr>
			<td style='font-weight: 700; font-size: 10pt; vertical-align: middle; color: navy; font-style: normal; font-family: Arial, sans-serif; white-space: nowrap; height: 38pt; text-decoration: none; text-align: general; border: medium none; padding-left: 1px; padding-right: 1px; padding-top: 1px' align='justify'>
			<span lang='en-us' style='LETTER-SPACING: 2pt'><font size='3'>Genes ",nom," -
			",taxoname," Annotations</font></span></td>
		</tr>
		<tr>
			<td style='font-weight: 400; font-size: 10pt; vertical-align: middle; width: 645pt; color: windowtext; font-style: normal; font-family: Arial; white-space: nowrap; text-decoration: none; text-align: general; border: medium none; padding-left: 1px; padding-right: 1px; padding-top: 1px' align='justify'>
			&nbsp;</td>
		</tr>
		<tr>
			<td style='font-weight: 400; font-size: 10pt; vertical-align: middle; width: 645pt; color: windowtext; font-style: normal; font-family: Arial; white-space: nowrap; text-decoration: none; text-align: general; border: medium none; padding-left: 1px; padding-right: 1px; padding-top: 1px' align='justify'>
			<font face='Times New Roman' size='3'>
			<b><span lang='en-us'>Used notations:</span></b></font></td>
		</tr>
		<tr>
			<td style='font-weight: 400; font-size: 10pt; vertical-align: middle; width: 645pt; color: windowtext; font-style: normal; font-family: Arial; white-space: nowrap; text-decoration: none; text-align: general; border: medium none; padding-left: 1px; padding-right: 1px; padding-top: 1px' align='justify'>
			&nbsp;</td>
		</tr>
		<tr>
			<td style='font-weight: 400; font-size: 10pt; vertical-align: middle; width: 645pt; color: windowtext; font-style: normal; font-family: Arial; white-space: nowrap; text-decoration: none; text-align: general; border: medium none; padding-left: 1px; padding-right: 1px; padding-top: 1px' align='justify'>
			<ul>
				<li><font face='Times New Roman' size='3'><span lang='en-us'><b>List Hits</b> - the number of genes
				annotated by the considered ",taxoname," category or annotation cluster within
				the analyzed list of target genes</span> </font> </li>
			</ul>
			</td>
		</tr>
		<tr>
			<td style='font-weight: 400; font-size: 10pt; vertical-align: middle; width: 645pt; color: windowtext; font-style: normal; font-family: Arial; white-space: nowrap; text-decoration: none; text-align: general; border: medium none; padding-left: 1px; padding-right: 1px; padding-top: 1px' align='justify'>
			<ul>
				<li><font face='Times New Roman' size='3'><span lang='en-us'><b>List </b></span><b>Total</b><span lang='en-us'>
				- the number of genes within the analyzed list of target genes
				having at least one ",taxoname," annotation</span> </font> </li>
			</ul>
			</td>
		</tr>
		<tr>
			<td style='font-weight: 400; font-size: 10pt; vertical-align: middle; width: 645pt; color: windowtext; font-style: normal; font-family: Arial; white-space: nowrap; text-decoration: none; text-align: general; border: medium none; padding-left: 1px; padding-right: 1px; padding-top: 1px' align='justify'>
			<ul>
				<li><font face='Times New Roman' size='3'><b>Population<span lang='en-us'> Hits</span></b><span lang='en-us'>
				- the number of genes, available on the entire microarray, annotated
				by the considered ",taxoname," category or annotation cluster</span>
				</font>
				</li>
			</ul>
			</td>
		</tr>
		<tr>
			<td style='font-weight: 400; font-size: 10pt; vertical-align: middle; width: 645pt; color: windowtext; font-style: normal; font-family: Arial; white-space: nowrap; text-decoration: none; text-align: general; border: medium none; padding-left: 1px; padding-right: 1px; padding-top: 1px' align='justify'>
			<ul>
				<li><font face='Times New Roman' size='3'><b>Population Total</b> - the number of genes available on the
				entire microarray and having at least one ",taxoname," annotation </font>
				</li>
			</ul>
			</td>
		</tr>
		<tr>
			<td style='font-weight: 400; font-size: 10pt; vertical-align: middle; width: 645pt; color: windowtext; font-style: normal; font-family: Arial; white-space: nowrap; text-decoration: none; text-align: general; border: medium none; padding-left: 1px; padding-right: 1px; padding-top: 1px' align='justify'>
			<ul>
				<li><font face='Times New Roman' size='3'><span lang='en-us'><b>P-value</b> - the significance p-value of
				the gene enrichment of the considered ",taxoname," category or annotation
				cluster, calculated with a unilateral Fisher exact test</span>
				</font> </li>
			</ul>
			</td>
		</tr>

	</table>
	<p>",sep="")


	write(head, file = paste(wd,"/",results.dir,"/Genes ",nom," - ",taxoname," Annotations.htm",sep=""), append = TRUE, sep = " ")

	table.head <- paste("
	<p>
..........................................................................................................................</p>
	<p><b>Terms</b></p>
	<table border='1' id='table2' style='border-collapse: collapse' bordercolor='#808080'>
		<tr>
			<td align='center' width='45'><b>Rank</b></td>
			<td align='center' width='70'><b>ID</b></td>
			<td align='center' width='200'><b>Term</b></td>
			<td align='center' width='45'>
			<p style='margin-top: 0; margin-bottom: 0'><b>List </b></p>
			<p style='margin-top: 0; margin-bottom: 0'><b>Hits</b></td>
			<td align='center' width='45'>
			<p style='margin-top: 0; margin-bottom: 0'><b>List </b></p>
			<p style='margin-top: 0; margin-bottom: 0'><b>Total</b></td>
			<td align='center' width='70'>
			<p style='margin-top: 0; margin-bottom: 0'><b>Population </b></p>
			<p style='margin-top: 0; margin-bottom: 0'><b>Hits</b></td>
			<td align='center' width='70'>
			<p style='margin-top: 0; margin-bottom: 0'><b>Population </b></p>
			<p style='margin-top: 0; margin-bottom: 0'><b>Total</b></td>
			<td align='center' width='100'><b>P-value</b></td>
			<td align='center' width='100'><b>Genes ID's</b></td>
			<td align='center' width='200'><b>Genes Names</b></td>
	</tr>",sep="")

	write(table.head, file = paste(wd,"/",results.dir,"/Genes ",nom," - ",taxoname," Annotations.htm",sep=""), append = TRUE, sep = " ")




	for(i in 1:length(exp.id)){

		if(taxoname == "GO Biological Process" || taxoname == "GO Cellular Component" || taxoname == "GO Molecular Function"){
			exp.term.lnk <- paste("<a target='_blank' href='http://www.ebi.ac.uk/ego/DisplayGoTerm?id=",exp.id[i],"' style='text-decoration: none'>
						<font size='3'>",exp.term[i],"</font></a>",sep="")
			exp.id.lnk <- paste("<a target='_blank' href='http://www.ebi.ac.uk/ego/DisplayGoTerm?id=",exp.id[i],"' style='text-decoration: none'>
						<font size='3'>",exp.id[i],"</font></a>",sep="")
		}else{
			exp.term.lnk <- paste("<font size='3'>",exp.term[i],"</font>",sep="")
			exp.id.lnk <- paste("<font size='3'>",exp.id[i],"</font>",sep="")
		}




	list.terms <- paste("
		<tr>
		<td align='center' width='45'><font size='3'>",i,"</font></td>
		<td align='center' width='70'>
		",exp.id.lnk,"</td>
		<td align='center' width='200'>
		",exp.term.lnk,"</td>
		<td align='center' width='45'><font size='3'>",exp.nr[i],"</font></td>
		<td align='center' width='45'><font size='3'>",exp.total[i],"</font></td>
		<td align='center' width='70'><font size='3'>",pop.hits[i],"</font></td>
		<td align='center' width='70'><font size='3'>",pop.total[i],"</font></td>
		<td align='center' width='70'><font size='3'>",format(term.p[i], digits = 3, scientific = TRUE),"</font></td>
		<td align='center' width='100'><textarea rows='0' name='Genes",i,"3' cols='7'>",list.genes[i],"</textarea></td>
		<td align='center' width='200'><select size='1' name='Genes",i,"4'>",list.names[i],"</select></td>
		</tr>",sep="")



	write(list.terms, file = paste(wd,"/",results.dir,"/Genes ",nom," - ",taxoname," Annotations.htm",sep=""), append = TRUE, sep = " ")

	}

	end.page <- paste("</table>
		<p>
		..........................................................................................................................</p>
		<p>
		This file was produced on ",date()," with <a target='_blank' href='http://corneliu.henegar.info/FunCluster.htm' style='text-decoration: none'>
		<font size='3'>FunCluster ",f.version,"</font></a> using GO and KEGG annotations updated on ",annot.date,".</p>
		</body>

		</html>",sep="")

	write(end.page, file = paste(wd,"/",results.dir,"/Genes ",nom," - ",taxoname," Annotations.htm",sep=""), append = TRUE, sep = " ")



}else{
	cat(paste("\n\t\tNo ",taxoname," ",nom," annotation results to save...",sep=""))
}
	rm()
}




#############################################################################################
#
# 18. Function PVALUES() -> Calculate the p-value for each term using an exact Fisher test
#
#############################################################################################

#Array DATAS contains all the data necessary to build the contingency table for calculating the Fisher test

pvalues <- function(datas){
	a <- as.matrix(datas[,1])
	b <- as.matrix(datas[,2]) - as.matrix(datas[,1])
	cc <- as.matrix(datas[,3]) - as.matrix(datas[,1])
	d <- as.matrix(datas[,4]) + as.matrix(datas[,1]) - as.matrix(datas[,2]) - as.matrix(datas[,3])
	p <- matrix(0,length(t(a)),1)
	for(i in 1:length(t(a))){
		p[i] <- as.double(matrix(fisher.test(matrix(c(a[i,1],b[i,1],cc[i,1],d[i,1]),2,2),alternative = "g"),1,1))
	}	

# test for p-values > 1 (to avoid errors during Storey q-values calculation)	

	for(i in 1:length(p)){
		if(p[i] > 1){
			p[i] <- 1
		}
	}

	return(p)
	rm()
}



#############################################################################################
#
# 19. FDR correction routine available from Chris Paciorek website
#
#############################################################################################


# File:      fdr.r 
# Date:      5/10/04
# Version:   0.1.3  
# Author:    Chris Paciorek - please contact the author with bug
#            reports: paciorek AT alumni.cmu.edu
# License:   GPL version 2 or later
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# Purpose:   implement False Discovery Rate (FDR) functions for multiple testing, following the Ventura et al. reference below
# Usage:     source('fdr.r');  fdr(my.pvals)
# References:
#             Ventura, V., C.J. Paciorek, and J.S. Risbey.  2004.  Controlling the proportion of falsely-rejected hypotheses when conducting multiple tests with climatological data.  Journal of Climate, in press.  Also Carnegie Mellon University, Department of Statistics technical report 775 (www.stat.cmu.edu/tr/tr775/tr775.html).
#             Benjamini, Y, and Y. Hochberg. 1995. Controlling the false discovery rate: a practical and powerful approach to multiple testing.  JRSSB 57:289-300.
#             Benjamini, Y. and D. Yekutieli.  2001.  The control of the false discovery rate in multiple testing under dependency. Annals of Statistics 29:1165-1188.
#             Benjamini, Y., A. Krieger, and D. Yekutieli.  2001.  Two staged linear step up FDR controlling procedure.  Technical Report, Department of Statistics and Operations Research, Tel Aviv University.  URL: http://www.math.tau.ac.il/~ybenja/Papers.html
#             Storey, J. 2002.  A direct approach to false discovery rates.  JRSSB 64: 479--498.
#

fdr <- function(pvals,qlevel=0.05,method="original",adjustment.method=NULL,adjustment.args=NULL){
#
# Description:
#
#    This is the main function designed for general usage for determining significance based on the FDR approach.
#
# Arguments:
#
#   pvals (required):  a vector of pvals on which to conduct the multiple testing
#
#   qlevel: the proportion of false positives desired
#
#   method: method for performing the testing.  'original' follows Benjamini & Hochberg (1995); 'general' is much more conservative, requiring no assumptions on the p-values (see Benjamini & Yekutieli (2001)).  We recommend using 'original', and if desired, using 'adjustment.method="mean" ' to increase power.
#
#   adjustment.method: method for increasing the power of the procedure by estimating the proportion of alternative p-values, one of "mean", the modified Storey estimator that we suggest in Ventura et al. (2004), "storey", the method of Storey (2002), or "two-stage", the iterative approach of Benjamini et al. (2001)
#
#   adjustment.args: arguments to adjustment.method; see prop.alt() for description, but note that for "two-stage", qlevel and fdr.method are taken from the qlevel and method arguments to fdr()
#
# Value:
#
#   NULL if no significant tests, or a vector of the indices of the significant tests
#
# Examples:
#
#   signif <- fdr(pvals,method="original",adjustment.method="mean")
#
  n <- length(pvals)

  a <- 0   # initialize proportion of alternative hypotheses
  if(!is.null(adjustment.method)){
    if(adjustment.method=="two-stage"){  # set things up for the "two-stage" estimator
      qlevel <- qlevel/(1+qlevel)  # see Benjamini et al. (2001) for proof that this controls the FDR at level qlevel
      adjustment.args$qlevel <- qlevel
      adjustment.args$fdr.method <- method
      cat(paste('Adjusting cutoff using two-stage method, with method ',method,' and qlevel ',round(qlevel,4),'\n',sep=""))
    }
    if(adjustment.method=="mean" & is.null(adjustment.args)){
      adjustment.args <- list(edf.lower=0.8,num.steps=20)  # default arguments for "mean" method of Ventura et al. (2004)
      cat(paste('Adjusting cutoff using mean method, with edf.lower=0.8 and num.steps=20\n',sep=""))
    }
    a <- prop.alt(pvals,adjustment.method,adjustment.args)
  }
  if(a==1){    # all hypotheses are estimated to be alternatives
    return(1:n)
  } else{      # adjust for estimate of a; default is 0
    qlevel <- qlevel/(1-a)
  }

  return(fdr.master(pvals,qlevel,method))
}

fdr.master <- function(pvals,qlevel=0.05,method="original"){
#
# Description:
#
#    This is an internal function that performs various versions of the FDR procedure, but without the modification described in section 4 of our J of Climate paper.
#
# Arguments:
#
#   pvals (required):  a vector of pvals on which to conduct the multiple testing
#
#   qlevel: the proportion of false positives desired
#
#   method: one of 'original', the original method of Benjamini & Hochberg (1995), or 'general', the method of Benjamini & Yekutieli (2001), which requires no assumptions about the p-values, but which is much more conservative.  We recommend 'original' for climatological data, and suspect it works well generally for spatial data.
#
# Value:
#
#   NULL if no significant tests, or a vector of the indices of the significant tests
#
  n <- length(pvals)
  if(method=="general"){
    qlevel <- qlevel/sum(1/(1:n))  # This requires fewer assumptions but is much more conservative
  } else{
    if(method!="original"){
      stop(paste("No method of type: ",method,sep=""))
    }
  }
  return(fdr.basic(pvals,qlevel))
}


fdr.basic <- function(pvals,qlevel=0.05){
#
# Description:
#
#    This is an internal function that performs the basic FDR of Benjamini & Hochberg (1995).
#
# Arguments:
#
#   pvals (required):  a vector of pvals on which to conduct the multiple testing
#
#   qlevel: the proportion of false positives desired
#
# Value:
#
#   NULL if no significant tests, or a vector of the indices of the significant tests
#
  n <- length(pvals)
  sorted.pvals <- sort(pvals)
  sort.index <- order(pvals)
  indices <- (1:n)*(sorted.pvals<=qlevel*(1:n)/n)
  num.reject <- max(indices)
  if(num.reject){
    indices <- 1:num.reject
    return(sort(sort.index[indices]))  
  } else{
    return(NULL)
  }
}


storey <- function(edf.quantile,pvals){
#
# Description:
#
#    This is an internal function that calculates the basic Storey (2002) estimator of a, the proportion of alternative hypotheses.
#
# Arguments:
#
#   edf.quantile (required):  the quantile of the empirical distribution function at which to estimate a
# 
#   pvals (required):  a vector of pvals on which to estimate a
#
# Value:
#
#   estimate of a, the number of alternative hypotheses
#
#
  if(edf.quantile >=1 | edf.quantile <=0){
    stop('edf.quantile should be between 0 and 1')
  }
  a <- (mean(pvals<=edf.quantile)-edf.quantile)/(1-edf.quantile)
  if(a>0){
    return(a)
  } else{
    return(0)
  }
}


prop.alt <- function(pvals,adjustment.method="mean",adjustment.args=list(edf.lower=0.8,num.steps=20)){
#
# Description:
#
#    This is an internal function that calculates an estimate of a, the proportion of alternative hypotheses, using one of several methods.
#
# Arguments:
#
#   pvals (required):  a vector of pvals from which to estimate a
#
#   adjustment.method: method for  estimating the proportion of alternative p-values, one of "mean", the modified Storey estimator suggested in Ventura et al. (2004); "storey", the method of Storey (2002); or "two-stage", the iterative approach of Benjamini et al. (2001)
#
#   adjustment.args: arguments to adjustment.method;
#      for "mean", specify edf.lower, the smallest quantile at which to estimate a, and num.steps, the number of quantiles to use - the approach uses the average of the Storey (2002) estimator for the num.steps quantiles starting at edf.lower and finishing just less than 1
#      for "storey", specify edf.quantile, the quantile at which to calculate the estimator
#      for "two-stage", the method uses a standard FDR approach to estimate which p-values are significant; this number is the estimate of a; therefore the method requires specification of qlevel, the proportion of false positives and fdr.method ('original' or 'general'), the FDR method to be used.  We do not recommend 'general' as this is very conservative and will underestimate a.
#  
# Value:
#
#   estimate of a, the number of alternative hypotheses
#
# Examples:
#
#   a <- prop.alt(pvals,adjustment.method="mean")
#
  n <- length(pvals)
  if(adjustment.method=="two-stage"){
    if(is.null(adjustment.args$qlevel) | is.null(adjustment.args$fdr.method)){
      stop("adjustment.args$qlevel or adjustment.args$fdr.method not specified.  Two-stage estimation of the number of alternative hypotheses requires specification of the FDR threshold and FDR method ('original' or 'general')")
    }
    return(length(fdr.master(pvals,adjustment.args$qlevel,method=adjustment.args$fdr.method))/n)
  }

  if(adjustment.method=="storey"){
    if(is.null(adjustment.args$edf.quantile)){
      stop("adjustment.args$edf.quantile not specified. Using Storey's method for estimating  the number of alternative hypotheses requires specification of the argument of the p-value EDF at which to do the estimation (a number close to one is recommended)")
    }
    return(storey(adjustment.args$edf.quantile,pvals))
  }

  if(adjustment.method=="mean"){
    if(is.null(adjustment.args$edf.lower) | is.null(adjustment.args$num.steps)){
      stop("adjustment.args$edf.lower or adjustment.args$num.steps is not specified. Using the method of Ventura et al. (2004) for estimating  the number of alternative hypotheses requires specification of the lowest quantile of the p-value EDF at which to do the estimation (a number close to one is recommended) and the number of steps between edf.lower and 1, starting at edf.lower, at which to do the estimation")
    }
    if(adjustment.args$edf.lower >=1 | adjustment.args$edf.lower<=0){
      stop("adjustment.args$edf.lower must be between 0 and 1");
    }
    if(adjustment.args$num.steps<1 | adjustment.args$num.steps%%1!=0){
      stop("adjustment.args$num.steps must be an integer greater than 0")
    }
    stepsize <- (1-adjustment.args$edf.lower)/adjustment.args$num.steps
    edf.quantiles <- matrix(seq(from=adjustment.args$edf.lower,by=stepsize,len=adjustment.args$num.steps),nr=adjustment.args$num.steps,nc=1)
    a.vec <- apply(edf.quantiles,1,storey,pvals)
    return(mean(a.vec))
  }
}

