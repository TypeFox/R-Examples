#############################################
### Code to create function "GetGeneList" ###
### Hanna and Riley                       ###
#############################################

GetGeneList = function(Species,build,featuretype=c("GENE","PSEUDO"),savefiles=FALSE,destfile){ 
	if(missing(Species)){
		stop("ERROR: No species specified")
	}
	if(missing(build)){
		stop("ERROR: No map build specified")
	}
	if(missing(destfile)){
		stop("ERROR: No path was specified for saving temporary and permanent files.")
	}
	cat("Please be patient, this could take a few minutes.","\n")
	Species = sub(" ","_",Species)
	dest = paste(destfile,"seq_gene.md.gz",sep="")
	URL = paste("ftp://ftp.ncbi.nih.gov/genomes/MapView/",Species,"/sequence/",build,"/initial_release/seq_gene.md.gz",sep="")
     	download.file(URL,dest,cacheOK=TRUE)
	cat("Reading in file...","\n")
	NCBIList<-read.table(gzfile(dest),header=FALSE,fill=TRUE) 
	if(ncol(NCBIList) > 15){
		NCBIList<-NCBIList[,1:15]
	}
	feature_type = group_label = NULL
  	colnames(NCBIList) = c("tax_id","chromosome","chr_start","chr_stop","chr_orient",
		"contig","ctg_start","ctg_stop","ctg_orient","feature_name","feature_id","feature_type","group_label","transcript",
		"evidence_code")
	ListF = matrix(0,1,15,dimnames=list(1,colnames(NCBIList)))
	for(j in featuretype){
			List = subset(NCBIList,feature_type==j)
			ListF = rbind(ListF,List)
	}
	ListF = ListF[2:nrow(ListF),]
	if(savefiles == TRUE){
		write.table(NCBIList,paste(destfile,"seq_gene.txt",sep=""),quote=FALSE,sep=" ",row.names=FALSE)
		remove(NCBIList)
	}
	if(savefiles == FALSE){
		unlink(dest)
		remove(NCBIList)
	}
	Assembly = matrix(unique(ListF[,'group_label']),ncol=1) 
	Features = matrix(unique(ListF[,'feature_type']),ncol=1)
	cat("Duplicate gene information may be present due to multiple assemblies and feature types.","\n",
		"The following assembly builds are present in this list:","\n")
		print(Assembly)
		cat("The following feature types are present in this list:","\n")
		print(Features)
		y = readline("Please choose which ASSEMBLY that you want to prioritize feature information from \n 
			(e.g. 1 for the first assembly listed, 2 for the second, etc.). \n 
			Other duplicate gene information (if any) will be removed from the list. \n ")
		y = as.numeric(y)
		if(abs(y) > nrow(Assembly)){
			stop("ERROR: You specified a number outside the range possible for the assemblies.")
		}
		if(nrow(Features) > 1){
			x = readline("Do you want to keep multiple feature type information? y = yes, n = no \n ")
		}else{ x = "y" }
		if(x == "n"){
			z = readline("Please choose which FEATURE TYPE that you want to prioritize information from \n 
				(e.g. 1 for the first feature listed, 2 for the second, etc.). \n
				Other duplicate information (if any) will be removed from the list. \n")
			z = as.numeric(z)
			if(abs(z) > nrow(Features)){
				stop("ERROR: You specified a number outside the range possible for the feature types.")
			}
		}
	##Finding Unique List of Genes Based on Assembly Preference ###
	ListAy = subset(ListF, group_label==Assembly[y,1])	
	ListAn = subset(ListF, !(ListF$feature_name %in% ListAy[,'feature_name']))
	ListUnA = rbind(ListAy,ListAn)

	##Finding Unique List of Genes Based on Whether Multiple Features are kept or not ###
	if(x == "n"){
		ListFez = subset(ListUnA,feature_type==Features[z,1])
		ListFen = subset(ListUnA,!(ListUnA$feature_name %in% ListFez[,'feature_name']))
		GeneList = rbind(ListFez,ListFen)
	}
	if(x == "y"){
		GeneList = ListUnA
	}
	if(x != "y"){
		if(x != "n"){
			stop("ERROR: You did not answer if you wanted duplicate feature type information removed.","\n",
				"Please start over and enter y for yes or n for no when prompted.")
		}
	}					
	GeneList=GeneList[order(GeneList$chromosome,GeneList$chr_start),] ##Sorting by Chr & Start Position##
	if(savefiles==TRUE){
		write.table(GeneList,paste(destfile,"GeneList.txt",sep=""),quote=FALSE,sep=" ",row.names=FALSE)
	}
	cat("Finished processing features and assemblies. The list will now be returned to the user.","\n")
	return(GeneList)
}