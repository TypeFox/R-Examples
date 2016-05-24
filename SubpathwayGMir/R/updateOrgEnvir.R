updateOrgEnvir <-
function(org="hsa",path="http://rest.kegg.jp",verbose=TRUE){
      initializeK2ri()
	  # assign("k2ri",new.env(parent=globalenv()),envir=.GlobalEnv)
	  print(paste("Update the current organism : ",org,sep=""))
	  
	  if(verbose==TRUE){
	      miRNA2Org <-  GetK2riData("miRNA2Org")
	      if(!org%in%levels(miRNA2Org[["Species"]])){
	      print(paste("The organism : ",org,", out of our analysis",collapse=NULL))
	      }
      }
      if(verbose==TRUE){
        print("Note that the programming may be time consumming!")  
      } 	 

	  if(verbose==TRUE){
	    print("Download relations between gene and symbol.")
		###  get ncbi gene id to gene symbol data
		file1 <- paste(path,"/","list","/",org,sep="")
        gene2symbol1<-read.table(file1,header=FALSE,sep = "\t", quote="\"",fill=TRUE,stringsAsFactors=FALSE)
		n1 <- grep("uncharacterized",gene2symbol1[,2])
		n2 <- grep("non-protein coding",gene2symbol1[,2])
		n3 <- grep("LINC",gene2symbol1[,2])
		n4 <- grep("LOC",gene2symbol1[,2])
		n  <- c(n1,n2,n3,n4)
		if(length(n)>0){ gene2symbol2 <- gene2symbol1[-n,] }else{ gene2symbol2 <- gene2symbol1}
	    gene2symbol3 <- gene2symbol2[which(gene2symbol2[,2]!=""),1:2]
        ncbi2symbol  <- data.frame(ncbi=sapply(gene2symbol3[,1], function(x) gsub(paste(org,":",sep=""),"",x),USE.NAMES=FALSE),
		             symbol=sapply(gene2symbol3[,2],function(x) unlist(strsplit(unlist(strsplit(unlist(strsplit(x,","))[[1]],";"))[[1]]," gene product from transcript "))[[1]],USE.NAMES=FALSE))
        
		###  get kegg id to ncbi gene id data
		file2 <- paste(path,"/","conv/ncbi-geneid","/",org,sep="")
        kegg2ncbi1 <- read.table(file2,header=FALSE,sep = "\t", quote="\"",fill=TRUE,stringsAsFactors=FALSE)
        kegg2ncbi  <- data.frame(ncbi=sapply(kegg2ncbi1[,1], function(x) gsub(paste(org,":",sep=""),"",x),USE.NAMES=FALSE),
		              gene=sapply(kegg2ncbi1[,2],function(x) gsub("ncbi-geneid:","",x),USE.NAMES=FALSE))
        ###  get kegg id to gene symbol 
		gene2symbol <- merge(kegg2ncbi,ncbi2symbol,by="ncbi")[,-1]
		k2ri$gene2symbol <- gene2symbol
		
		###  get background gene data 
		BGGene <- unique(sapply(strsplit(as.character(kegg2ncbi1[,2]),":"),function(x) x[2]))
        k2ri$BGGene <- BGGene
	   }
	  	 
	  
	   if(verbose==TRUE){
	    print("Download relations between KEGG gene and pathway")
        file2 <- paste(path,"/","link","/","pathway","/",org,sep="")
	    gene2path <- read.table(file2,header=FALSE,sep = "\t", quote="\"",stringsAsFactors =FALSE)	  
        k2ri$gene2path<-gene2path
	  }
	  
	  if(verbose==TRUE){
	   print("Download background of miRNAs")
	   miRNA2Org <- GetK2riData("miRNA2Org")
       BGMiRNA <- unique(as.character(miRNA2Org[which(as.character(miRNA2Org[["Species"]])==org),1]))
       k2ri$BGMiRNA <- BGMiRNA

	  }
	  
	    if(verbose==TRUE){
	   expMir2Tar <- GetK2riData("expMir2Tar")
       k2ri$expMir2Tar<-expMir2Tar
	  } 
	  
	   if(verbose==TRUE){
	   miRNA2Org <- GetK2riData("miRNA2Org")
       k2ri$miRNA2Org<-miRNA2Org
	  } 
	   
	   if(verbose==TRUE){
	   print("Download background of direct KEGG metabolic pathways")
       MetabolicGEGEEMGraph <- paste(toupper(org),"_MetabolicGEGEEMGraph",sep="")
	   MetabolicGEGEEMGraph <- GetK2riData(MetabolicGEGEEMGraph)
	   k2ri$MetabolicGEGEEMGraph<-MetabolicGEGEEMGraph
	  } 
	  
	   if(verbose==TRUE){
	    print("Download background of undirect KEGG metabolic pathways")
	   MetabolicGEGEUEMGraph <- paste(toupper(org),"_MetabolicGEGEUEMGraph",sep="")
	   MetabolicGEGEUEMGraph <- GetK2riData(MetabolicGEGEUEMGraph)
	   k2ri$MetabolicGEGEUEMGraph<-MetabolicGEGEUEMGraph
	  } 
     # save(k2ri,file="k2ri.rda", envir=.GlobalEnv)
	

}
