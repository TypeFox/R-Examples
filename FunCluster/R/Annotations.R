.packageName <- "FunCluster"

MiseEnFormeKegg <- function ( ligne, Fichier) { 
	vec1 <- ligne
	LL <-as.character( vec1[1]) 
	if (!identical(all.equal( as.numeric(LL), as.numeric(NA)),TRUE) ){   
		Kmap <- as.character(vec1[-1])
		Kmap <- Kmap [ Kmap != ""]  
		NbKmap <- length (Kmap)

		LocusLink  <- gl ( 1, NbKmap, label = LL) 

		LLcorrespondance <- cbind.data.frame( LocusLink, Kmap )  
		write.table ( LLcorrespondance, Fichier ,col.names =FALSE,row.names=FALSE,append = TRUE, sep= "\t",quote = FALSE)  


	}

}


MiseEnFormeKeggSC <- function ( ligne, Fichier) { 
	vec1 <- ligne
	LL <-as.character( vec1[1]) 

	Kmap <- as.character(vec1[-1])
	Kmap <- Kmap [ Kmap != ""]  
	NbKmap <- length (Kmap)

	LocusLink  <- gl ( 1, NbKmap, label = LL) 

	LLcorrespondance <- cbind.data.frame( LocusLink, Kmap )  
	write.table ( LLcorrespondance, Fichier ,col.names =FALSE,row.names=FALSE,append = TRUE, sep= "\t",quote = FALSE)  
}




ChargerKegg <- function(espece = espece) {  


	if ( espece == "hsa" ) { 
		DIResp <- "HS"
	} else {
		if ( espece == "mmu"){   
			DIResp <- "MM"
		}else if(espece == "sce"){
			DIResp <- "SC"
		}
	}


	download.file ( paste("ftp://ftp.genome.ad.jp/pub/kegg/pathway/organisms/",espece,"/",espece,"_gene_map.tab",sep="") , paste(getwd(),"/Annotations/KEGG/",DIResp,"/tmp.txt",sep=""), mode = "w" )

	FichierFinalNonQuote <- file (paste(getwd(),"/Annotations/KEGG/",DIResp,"/ll_keggNonQuote.txt",sep=""), "w+")

	FichierFinal <- file (paste(getwd(),"/Annotations/KEGG/",DIResp,"/ll_kegg.txt",sep=""), "w")

	tempo <- file (paste(getwd(),"/Annotations/KEGG/",DIResp,"/tmp.txt",sep=""), "r") 

	tempo2 <- file (paste(getwd(),"/Annotations/KEGG/",DIResp,"/tmp2.txt",sep=""), "w+")






	NbCol <- count.fields(file=tempo)
	MaxNbCol <- max (NbCol)


	seek ( tempo, where = 0)


	DTtemp <- read.table ( file = tempo ,na.strings = "",fill=TRUE,colClasses="character",sep= "\t",header = FALSE,quote = "",comment.char = "")

	write.table ( DTtemp, tempo2 ,col.names =FALSE,row.names=FALSE,quote =FALSE, sep = "\t")
	DTtemp <- read.table ( file = tempo2 ,fill=TRUE,colClasses="character",col.names =c(1: (MaxNbCol)),sep= "", quote = "",comment.char = "")



	close(tempo)
	close(tempo2)
	unlink(paste(getwd(),"/Annotations/KEGG/",DIResp,"/tmp.txt",sep=""))
	unlink (paste(getwd(),"/Annotations/KEGG/",DIResp,"/tmp2.txt",sep=""))





	if(espece == "hsa" || espece == "mmu"){		
		suppressWarnings (apply (DTtemp, 1, FUN = "MiseEnFormeKegg" ,FichierFinalNonQuote))
	}else if (espece == "sce"){
		suppressWarnings (apply (DTtemp, 1, FUN = "MiseEnFormeKeggSC" ,FichierFinalNonQuote))
	}

	if ( espece == "hsa" ) { 

		HS.KEGG.file.annot <<- read.table ( file = FichierFinalNonQuote ,na.strings = "",fill=TRUE,colClasses="character",sep= "\t",header = FALSE,quote = "",comment.char = "")
		write.table (HS.KEGG.file.annot , FichierFinal ,col.names =FALSE,row.names=FALSE,append = TRUE, sep= "\t") 

	} else if ( espece == "mmu"){

		MM.KEGG.file.annot <<- read.table ( file = FichierFinalNonQuote ,na.strings = "",fill=TRUE,colClasses="character",sep= "\t",header = FALSE,quote = "",comment.char = "")
		write.table (MM.KEGG.file.annot , FichierFinal ,col.names =FALSE,row.names=FALSE,append = TRUE, sep= "\t") 
	
	}else if(espece == "sce"){
			
		SC.KEGG.file.annot <<- read.table ( file = FichierFinalNonQuote ,na.strings = "",fill=TRUE,colClasses="character",sep= "\t",header = FALSE,quote = "",comment.char = "")
			
		DTLL <- GeneInfoSC()
		DTLL <- as.matrix(DTLL[DTLL[,1] == "4932",])
		rownames(DTLL) <- as.character(DTLL[,4])
		SC.KEGG.file.annot <<- as.matrix(SC.KEGG.file.annot)
		rownames(SC.KEGG.file.annot) <<- as.character(SC.KEGG.file.annot[,1])
		x <- DTLL[DTLL[,4] %in% as.character(SC.KEGG.file.annot[,1]),2]
		SC.KEGG.file.annot <<- SC.KEGG.file.annot[SC.KEGG.file.annot[,1] %in% names(x),]
			
		y <- NULL
			
		for(i in 1:nrow(SC.KEGG.file.annot)){				
				
			y <- rbind(y,c(x[SC.KEGG.file.annot[i,1]],SC.KEGG.file.annot[i,2]))
		}
		colnames(y) <- as.character(matrix("",1,ncol(y)))
		SC.KEGG.file.annot <<- y	
			
		write.table (SC.KEGG.file.annot , FichierFinal ,col.names =FALSE, row.names=FALSE,append = TRUE, sep= "\t") 

	}


	close (FichierFinalNonQuote)
	unlink(paste(getwd(),"/Annotations/KEGG/",DIResp,"/ll_keggNonQuote.txt",sep=""))
	close (FichierFinal)


}



Kegg <- function () { 


	dir.create(paste(getwd(),"/Annotations/KEGG",sep=""))


	download.file ( "ftp://ftp.genome.ad.jp/pub/kegg/pathway/map_title.tab" ,paste(getwd(), "/Annotations/KEGG/Temp_kegg_terms.txt",sep=""), mode = "w" )
	temp <- file (paste(getwd(),"/Annotations/KEGG/Temp_kegg_terms.txt",sep=""),"r")
	Fichkeggterms <- file (paste(getwd(),"/Annotations/KEGG/kegg_terms.txt",sep=""),"w")

	KEGG.terms.name <<- read.table (file = temp ,na.strings = "",fill=TRUE,colClasses="character",sep= "\t",header = FALSE,quote = "",comment.char = "")
	write.table ( KEGG.terms.name, Fichkeggterms ,col.names =FALSE,row.names=FALSE, sep = "\t")

	close (temp)
	close (Fichkeggterms)
	unlink ( paste(getwd(),"/Annotations/KEGG/Temp_kegg_terms.txt",sep=""))


	dir.create (paste(getwd(),"/Annotations/KEGG/HS",sep=""))
	ChargerKegg ("hsa")


	dir.create (paste(getwd(),"/Annotations/KEGG/MM",sep=""))
	ChargerKegg ("mmu")


	dir.create (paste(getwd(),"/Annotations/KEGG/SC",sep=""))
	ChargerKegg ("sce")
}



MakeGoTerms <- function (DTGoTermsAndIds) {


	DTGoTermsAndIds <- DTGoTermsAndIds [ substr ( DTGoTermsAndIds[,1],1,1) != "!" ,]
	dimens <- nrow(DTGoTermsAndIds)
	DTGoTermsAndIds <- cbind.data.frame( DTGoTermsAndIds, row.names = c(1: dimens)) 
	GO.terms.name <<- DTGoTermsAndIds [,1:2] 



	FichierFinal <- file (paste(getwd(),"/Annotations/GO/go_terms.txt",sep=""), "w") 
	write.table ( GO.terms.name , FichierFinal ,col.names =FALSE,row.names=FALSE, sep = "\t")


	close ( FichierFinal )

}



MakeGoHierarchy <- function (DATE="") {


	if(DATE == ""){
		DATE <- Sys.time()
		DATE <- substr( as.character(DATE),1,7)  
		DATE <- gsub("-","",DATE)   
	}
	NomFichier  <- paste ( "go_", DATE , "-termdb-tables.tar.gz" ,sep = "")




	download.file(paste ("http://godatabase-archive.stanford.edu/latest/", NomFichier,sep = ""),paste (getwd(),"/go_DATE-termdb-tables.tar.gz",sep=""),mode="wb")

	system("gunzip go_DATE-termdb-tables.tar.gz")
	system("tar -xf go_DATE-termdb-tables.tar")



	FichTerm2term <- file (paste(getwd(),"/go_",DATE,"-termdb-tables/term2term.txt",sep=""), "r") 
	FichTerm <- file (paste(getwd(),"/go_",DATE,"-termdb-tables/term.txt",sep=""), "r") 
	FichFinal <- file (paste(getwd(),"/Annotations/GO/go_hierarchy.txt",sep=""), "w") 


	DTterm2term <- read.table ( file = FichTerm2term ,na.strings = "",fill=TRUE,colClasses="character",col.names =c(1:5),sep= "\t",header = FALSE, quote = "",comment.char = "")
	DTterm <- read.table ( file = FichTerm ,na.strings = "",fill=TRUE,colClasses="character",col.names =c(1:6),sep= "\t",header = FALSE, quote = "",comment.char = "")
	GO.terms.hierarchy <<- cbind ( DTterm [ DTterm2term[,4],4],DTterm [ DTterm2term[,3],4] )


	write.table ( GO.terms.hierarchy, FichFinal ,col.names =FALSE,row.names=FALSE, sep= "\t")


	close (	FichTerm2term )
	close(FichTerm )
	close (FichFinal)
	unlink (paste (getwd(),"/go_DATE-termdb-tables.tar",sep=""))
	unlink (paste(getwd(),"/go_",DATE,"-termdb-tables",sep=""), recursive = TRUE)

	return (GO.terms.hierarchy)

}



AnnotationIndirecte <- function (DTGOOntology,DTGOHierarchy , DIResp, NomFichFinal ) {  


	FichFinal <- file (paste(getwd(),"/Annotations/GO/",DIResp,"/IND/",NomFichFinal,".txt",sep= ""), "w") 


	DTIND <- merge ( DTGOOntology, DTGOHierarchy , by.x= 2,by.y = 1)

	names(DTIND)<-c("X0",names(DTGOOntology))
	DTIND <- rbind.data.frame ( DTGOOntology, DTIND[,2:3] )

	OrdreLignes <- do.call("order", DTIND[,c( 1,2)])
	DTIND <- DTIND [ OrdreLignes,]

	DTIND <- unique.data.frame (DTIND) 

	dimens <- nrow(DTIND)
	DTIND <- cbind.data.frame( DTIND, row.names = c(1: dimens))

	write.table ( DTIND, FichFinal ,col.names =FALSE,row.names=FALSE,sep="\t")

	close( FichFinal )
	return (DTIND)

}




TroisiemeFiltrage <- function (DTGoTermsAndIds,DTLoc2GoFiltre1, FichOntology ,initialeOntology) {

	OntDTGoTermsAndIds <- DTGoTermsAndIds[ DTGoTermsAndIds[,3] == initialeOntology,]
	DTOntology <- DTLoc2GoFiltre1 [ DTLoc2GoFiltre1[,2] %in% OntDTGoTermsAndIds[,1],]

	DTOntology <- unique.data.frame (DTOntology) 

	OrdreLignes <- do.call("order", DTOntology[,c( 1,2)])  
	DTOntology <- DTOntology [ OrdreLignes,]

	dimens <- nrow(DTOntology)
	DTOntology <- cbind.data.frame( DTOntology, row.names = c(1: dimens))

	write.table ( DTOntology, FichOntology ,col.names =FALSE,row.names=FALSE,sep="\t")
	return (DTOntology)
}



GeneInfo <- function(){


	download.file ( paste("ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz",sep="") ,paste( getwd(),"/Annotations/Temp.gz",sep=""), mode = "wb" )

	LLTempgz <- gzfile (paste(getwd(),"/Annotations/Temp.gz",sep=""),"rb") 

	LLTemp <- file (paste(getwd(),"/Annotations/Temp.txt",sep=""), "w+")
	vecTemp <- readLines(LLTempgz)
	cat ( vecTemp , file = LLTemp, sep = "\n")

	rm(vecTemp)


	DTLLbrut <- read.table ( file = LLTemp,na.strings = "-",fill=TRUE,colClasses="character",sep= "\t",header = FALSE, quote = "",comment.char = "#")
	DTLL <- cbind(DTLLbrut[,1],DTLLbrut[,2],DTLLbrut[,7],DTLLbrut[,8],DTLLbrut[,9])

	close(LLTempgz)
	close (LLTemp)
	unlink (paste(getwd(),"/Annotations/Temp.gz",sep=""))
	unlink (paste(getwd(),"/Annotations/Temp.txt",sep=""))


	return (DTLL) 


}



GeneInfoSC <- function(){


	download.file ( paste("ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz",sep="") ,paste( getwd(),"/Annotations/Temp.gz",sep=""), mode = "wb" )

	LLTempgz <- gzfile (paste(getwd(),"/Annotations/Temp.gz",sep=""),"rb") 

	LLTemp <- file (paste(getwd(),"/Annotations/Temp.txt",sep=""), "w+")
	vecTemp <- readLines(LLTempgz)
	cat ( vecTemp , file = LLTemp, sep = "\n")

	rm(vecTemp)


	DTLLbrut <- read.table ( file = LLTemp,na.strings = "-",fill=TRUE ,colClasses="character",sep= "\t",header = FALSE, quote = "",comment.char = "#")
	DTLL <- cbind(DTLLbrut[,1],DTLLbrut[,2],DTLLbrut[,3],DTLLbrut[,4],DTLLbrut[,7],DTLLbrut[,8],DTLLbrut[,9])

	close(LLTempgz)
	close (LLTemp)
	unlink (paste(getwd(),"/Annotations/Temp.gz",sep=""))
	unlink (paste(getwd(),"/Annotations/Temp.txt",sep=""))


	return (DTLL) 


}





MakeCorrespondanceLLGO <- function (DTLoc2Go, DTGoTermsAndIds,DTGOHierarchy, espece = espece) {  

	if ( espece == "hs" ) {
		DIResp <- "HS"
		DTLoc2Go <- DTLoc2Go[DTLoc2Go[,1] == "9606",2:3]
	} else {
		if ( espece == "mm"){
			DIResp <- "MM"
			DTLoc2Go <- DTLoc2Go[DTLoc2Go[,1] == "10090",2:3]
		}else if (espece == "sc"){
			DIResp <- "SC"
			DTLoc2Go <- DTLoc2Go[DTLoc2Go[,1] == "4932",2:3]
		}
	}



	FichBioProcess <- file (paste(getwd(),"/Annotations/GO/",DIResp,"/DIR/biological_process.txt",sep= ""), "w") 
	FichCellComp <- file (paste(getwd(),"/Annotations/GO/",DIResp,"/DIR/cellular_component.txt",sep= ""), "w") 
	FichMolFunc <- file (paste(getwd(),"/Annotations/GO/",DIResp,"/DIR/molecular_function.txt",sep= ""), "w")



	DTLoc2GoFiltre1 <- DTLoc2Go[ DTLoc2Go[,2] %in% DTGoTermsAndIds[,1],1:2]





	if ( espece == "hs" ) {
		HS.GO.DIR.BP.file.annot <<- TroisiemeFiltrage (DTGoTermsAndIds,DTLoc2GoFiltre1, FichBioProcess ,"P")
		HS.GO.DIR.CC.file.annot <<- TroisiemeFiltrage (DTGoTermsAndIds,DTLoc2GoFiltre1, FichCellComp ,"C")
		HS.GO.DIR.MF.file.annot <<- TroisiemeFiltrage (DTGoTermsAndIds,DTLoc2GoFiltre1, FichMolFunc ,"F")
		HS.GO.IND.BP.file.annot <<- AnnotationIndirecte (HS.GO.DIR.BP.file.annot,DTGOHierarchy,DIResp,"biological_process")
		HS.GO.IND.CC.file.annot <<- AnnotationIndirecte (HS.GO.DIR.CC.file.annot,DTGOHierarchy,DIResp,"cellular_component")
		HS.GO.IND.MF.file.annot <<- AnnotationIndirecte (HS.GO.DIR.MF.file.annot,DTGOHierarchy,DIResp,"molecular_function")
	} else {
		if ( espece == "mm"){
			MM.GO.DIR.BP.file.annot <<- TroisiemeFiltrage (DTGoTermsAndIds,DTLoc2GoFiltre1, FichBioProcess ,"P")
			MM.GO.DIR.CC.file.annot <<- TroisiemeFiltrage (DTGoTermsAndIds,DTLoc2GoFiltre1, FichCellComp ,"C")
			MM.GO.DIR.MF.file.annot <<- TroisiemeFiltrage (DTGoTermsAndIds,DTLoc2GoFiltre1, FichMolFunc ,"F")
			MM.GO.IND.BP.file.annot <<- AnnotationIndirecte (MM.GO.DIR.BP.file.annot,DTGOHierarchy,DIResp,"biological_process")
			MM.GO.IND.CC.file.annot <<- AnnotationIndirecte (MM.GO.DIR.CC.file.annot,DTGOHierarchy,DIResp,"cellular_component")
			MM.GO.IND.MF.file.annot <<- AnnotationIndirecte (MM.GO.DIR.MF.file.annot,DTGOHierarchy,DIResp,"molecular_function")
		}else if(espece == "sc"){
			SC.GO.DIR.BP.file.annot <<- TroisiemeFiltrage (DTGoTermsAndIds,DTLoc2GoFiltre1, FichBioProcess ,"P")
			SC.GO.DIR.CC.file.annot <<- TroisiemeFiltrage (DTGoTermsAndIds,DTLoc2GoFiltre1, FichCellComp ,"C")
			SC.GO.DIR.MF.file.annot <<- TroisiemeFiltrage (DTGoTermsAndIds,DTLoc2GoFiltre1, FichMolFunc ,"F")
			SC.GO.IND.BP.file.annot <<- AnnotationIndirecte (SC.GO.DIR.BP.file.annot,DTGOHierarchy,DIResp,"biological_process")
			SC.GO.IND.CC.file.annot <<- AnnotationIndirecte (SC.GO.DIR.CC.file.annot,DTGOHierarchy,DIResp,"cellular_component")
			SC.GO.IND.MF.file.annot <<- AnnotationIndirecte (SC.GO.DIR.MF.file.annot,DTGOHierarchy,DIResp,"molecular_function")


		}
	}


	close(FichBioProcess)
	close(FichCellComp)
	close (FichMolFunc)


}



MakeLocusNames <- function( espece, DTLL ) {  

	if ( espece == "hs" ) {
		DTLLesp <- DTLL[DTLL[,1] == "9606",2:5]
		DIResp <- "HS"
	} else {
		if ( espece == "mm"){
			DTLLesp <- DTLL[DTLL[,1] == "10090",2:5]
			DIResp <- "MM"
		}else if(espece == "sc"){
			DTLLesp <- DTLL[DTLL[,1] == "4932",2:5]
			DIResp <- "SC"
		}
	}



	DTLLloc <- DTLLesp[, c(1,2,3)]
	DTLLloc <- unique.data.frame (DTLLloc) 
	dimens <- nrow(DTLLloc)
	DTLLloc <- cbind.data.frame( DTLLloc, row.names = c(1: dimens)) #
	DTLLloc <- DTLLloc[2:nrow(DTLLloc),]

	FichierLocation <- file (paste(getwd(),"/Annotations/LL/",DIResp,"/location.txt",sep= ""), "w") 

	write.table ( DTLLloc, FichierLocation ,col.names =FALSE,row.names=FALSE, sep="\t")


	DTLLesp <- DTLLesp[, c(1,4)]
	DTLLesp <- unique.data.frame (DTLLesp) 
	dimens <- nrow(DTLLesp)
	DTLLesp <- cbind.data.frame( DTLLesp, row.names = c(1: dimens)) 
	DTLLesp <- DTLLesp[2:nrow(DTLLesp),]




	FichierFinal <- file (paste(getwd(),"/Annotations/LL/",DIResp,"/locus_name.txt",sep= ""), "w") 

	write.table ( DTLLesp, FichierFinal ,col.names =FALSE,row.names=FALSE, sep="\t")


	if ( espece == "hs" ) {
		HS.locus.name <<- DTLLesp
		HS.locus.cyto <<- DTLLloc
	} else {
		if ( espece == "mm"){
			MM.locus.name <<- DTLLesp
			MM.locus.cyto <<- DTLLloc
		}else if(espece == "sc"){
			SC.locus.name <<- DTLLesp
			SC.locus.cyto <<- DTLLloc

		}
	}

	close ( FichierFinal )
	close(FichierLocation)

}




goAndLL <- function (DATE="") { 


	dir.create (paste(getwd(),"/Annotations/GO",sep=""))



	download.file ( "http://www.geneontology.org/doc/GO.terms_and_ids" ,paste(getwd(),"/Annotations/GO/GoTermsTmp.txt",sep=""), mode = "w" )

	FichGoTermsAndIds <- file (paste(getwd(),"/Annotations/GO/GoTermsTmp.txt",sep=""), "r") 

	DTGoTermsAndIds <- read.table ( file = FichGoTermsAndIds ,na.strings = "",fill=TRUE,colClasses="character",col.names =c(1:3),sep= "\t",header = FALSE, quote = "",comment.char = "")
	close(FichGoTermsAndIds)
	unlink(paste(getwd(),"/Annotations/GO/GoTermsTmp.txt",sep=""))





	download.file ( "ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz" ,paste(getwd(),"/Annotations/GO/gene2go.gz",sep=""), mode = "wb" )

	Gene2GOTempgz <- gzfile (paste(getwd(),"/Annotations/GO/gene2go.gz",sep=""),"rb") 

	Gene2GOTemp <- file (paste(getwd(),"/Annotations/GO/gene2go.txt",sep=""), "w+")
	vecTemp <- readLines(Gene2GOTempgz)
	cat ( vecTemp , file = Gene2GOTemp, sep = "\n")


	DTLoc2Go <- read.table ( file = Gene2GOTemp,na.strings = "",fill=TRUE ,colClasses="character",col.names =c(1:8),sep= "\t",header = FALSE, quote = "",comment.char = "#")
	close (Gene2GOTemp)
	unlink(paste(getwd(),"/Annotations/GO/gene2go.txt",sep=""))
	close (Gene2GOTempgz)
	unlink(paste(getwd(),"/Annotations/GO/gene2go.gz",sep=""))



	MakeGoTerms (DTGoTermsAndIds) 

	DTGoHierarchy <- MakeGoHierarchy (DATE=DATE)


	DTLL <- GeneInfo()


	dir.create (paste(getwd(),"/Annotations/GO/HS",sep=""))
	dir.create (paste(getwd(),"/Annotations/GO/HS/DIR",sep=""))
	dir.create (paste(getwd(),"/Annotations/GO/HS/IND",sep=""))
	MakeCorrespondanceLLGO (DTLoc2Go,DTGoTermsAndIds, DTGoHierarchy , "hs")


	dir.create (paste(getwd(),"/Annotations/GO/MM",sep=""))
	dir.create (paste(getwd(),"/Annotations/GO/MM/DIR",sep=""))
	dir.create (paste(getwd(),"/Annotations/GO/MM/IND",sep=""))
	MakeCorrespondanceLLGO (DTLoc2Go, DTGoTermsAndIds, DTGoHierarchy , "mm")


	dir.create (paste(getwd(),"/Annotations/GO/SC",sep=""))
	dir.create (paste(getwd(),"/Annotations/GO/SC/DIR",sep=""))
	dir.create (paste(getwd(),"/Annotations/GO/SC/IND",sep=""))
	MakeCorrespondanceLLGO (DTLoc2Go, DTGoTermsAndIds, DTGoHierarchy , "sc")




	dir.create(paste(getwd(),"/Annotations/LL",sep=""))


	dir.create (paste(getwd(),"/Annotations/LL/HS",sep=""))
	MakeLocusNames ( "hs", DTLL)


	dir.create (paste(getwd(),"/Annotations/LL/MM",sep=""))
	MakeLocusNames ( "mm", DTLL)


	dir.create (paste(getwd(),"/Annotations/LL/SC",sep=""))
	MakeLocusNames ( "sc", DTLL)
}


annotations <- function (date.annot="") {   


	dir.create (paste (getwd(),"/Annotations", sep= ""))
	Kegg ()
	goAndLL(DATE=date.annot)
	annot.date <- format(Sys.time(), "%Y %b %d")
	
# 	GO.terms.hierarchy <-GO.terms.hierarchy[(GO.terms.hierarchy[,1] %in% GO.terms.name[,1]) && (GO.terms.hierarchy[,2] %in% GO.terms.name[,1]),]
# 	
# 	HS.GO.DIR.BP.file.annot <- HS.GO.DIR.BP.file.annot[(HS.GO.DIR.BP.file.annot[,2] %in% GO.terms.name[,1]) & (HS.GO.DIR.BP.file.annot[,1] %in% HS.locus.name[,1]),]
# 	HS.GO.DIR.CC.file.annot <- HS.GO.DIR.CC.file.annot[(HS.GO.DIR.CC.file.annot[,2] %in% GO.terms.name[,1]) & (HS.GO.DIR.CC.file.annot[,1] %in% HS.locus.name[,1]),]
# 	HS.GO.DIR.MF.file.annot <- HS.GO.DIR.MF.file.annot[(HS.GO.DIR.MF.file.annot[,2] %in% GO.terms.name[,1]) & (HS.GO.DIR.MF.file.annot[,1] %in% HS.locus.name[,1]),]
# 	HS.GO.IND.BP.file.annot <- HS.GO.IND.BP.file.annot[(HS.GO.IND.BP.file.annot[,2] %in% GO.terms.name[,1]) & (HS.GO.IND.BP.file.annot[,1] %in% HS.locus.name[,1]),]
# 	HS.GO.IND.CC.file.annot <- HS.GO.IND.CC.file.annot[(HS.GO.IND.CC.file.annot[,2] %in% GO.terms.name[,1]) & (HS.GO.IND.CC.file.annot[,1] %in% HS.locus.name[,1]),]
# 	HS.GO.IND.MF.file.annot <- HS.GO.IND.MF.file.annot[(HS.GO.IND.MF.file.annot[,2] %in% GO.terms.name[,1]) & (HS.GO.IND.MF.file.annot[,1] %in% HS.locus.name[,1]),]
# 	
# 	MM.GO.DIR.BP.file.annot <- MM.GO.DIR.BP.file.annot[(MM.GO.DIR.BP.file.annot[,2] %in% GO.terms.name[,1]) & (MM.GO.DIR.BP.file.annot[,1] %in% MM.locus.name[,1]),]
# 	MM.GO.DIR.CC.file.annot <- MM.GO.DIR.CC.file.annot[(MM.GO.DIR.CC.file.annot[,2] %in% GO.terms.name[,1]) & (MM.GO.DIR.CC.file.annot[,1] %in% MM.locus.name[,1]),]
# 	MM.GO.DIR.MF.file.annot <- MM.GO.DIR.MF.file.annot[(MM.GO.DIR.MF.file.annot[,2] %in% GO.terms.name[,1]) & (MM.GO.DIR.MF.file.annot[,1] %in% MM.locus.name[,1]),]
# 	MM.GO.IND.BP.file.annot <- MM.GO.IND.BP.file.annot[(MM.GO.IND.BP.file.annot[,2] %in% GO.terms.name[,1]) & (MM.GO.IND.BP.file.annot[,1] %in% MM.locus.name[,1]),]
# 	MM.GO.IND.CC.file.annot <- MM.GO.IND.CC.file.annot[(MM.GO.IND.CC.file.annot[,2] %in% GO.terms.name[,1]) & (MM.GO.IND.CC.file.annot[,1] %in% MM.locus.name[,1]),]
# 	MM.GO.IND.MF.file.annot <- MM.GO.IND.MF.file.annot[(MM.GO.IND.MF.file.annot[,2] %in% GO.terms.name[,1]) & (MM.GO.IND.MF.file.annot[,1] %in% MM.locus.name[,1]),]
# 	
# 	SC.GO.DIR.BP.file.annot <- SC.GO.DIR.BP.file.annot[(SC.GO.DIR.BP.file.annot[,2] %in% GO.terms.name[,1]) & (SC.GO.DIR.BP.file.annot[,1] %in% SC.locus.name[,1]),]
# 	SC.GO.DIR.CC.file.annot <- SC.GO.DIR.CC.file.annot[(SC.GO.DIR.CC.file.annot[,2] %in% GO.terms.name[,1]) & (SC.GO.DIR.CC.file.annot[,1] %in% SC.locus.name[,1]),]
# 	SC.GO.DIR.MF.file.annot <- SC.GO.DIR.MF.file.annot[(SC.GO.DIR.MF.file.annot[,2] %in% GO.terms.name[,1]) & (SC.GO.DIR.MF.file.annot[,1] %in% SC.locus.name[,1]),]
# 	SC.GO.IND.BP.file.annot <- SC.GO.IND.BP.file.annot[(SC.GO.IND.BP.file.annot[,2] %in% GO.terms.name[,1]) & (SC.GO.IND.BP.file.annot[,1] %in% SC.locus.name[,1]),]
# 	SC.GO.IND.CC.file.annot <- SC.GO.IND.CC.file.annot[(SC.GO.IND.CC.file.annot[,2] %in% GO.terms.name[,1]) & (SC.GO.IND.CC.file.annot[,1] %in% SC.locus.name[,1]),]
# 	SC.GO.IND.MF.file.annot <- SC.GO.IND.MF.file.annot[(SC.GO.IND.MF.file.annot[,2] %in% GO.terms.name[,1]) & (SC.GO.IND.MF.file.annot[,1] %in% SC.locus.name[,1]),]
	
	save(GO.terms.hierarchy,GO.terms.name,HS.GO.DIR.BP.file.annot,HS.GO.DIR.CC.file.annot,HS.GO.DIR.MF.file.annot,
		HS.GO.IND.BP.file.annot,HS.GO.IND.CC.file.annot,HS.GO.IND.MF.file.annot,HS.KEGG.file.annot,HS.locus.name,HS.locus.cyto,KEGG.terms.name,
		MM.GO.DIR.BP.file.annot,MM.GO.DIR.CC.file.annot,MM.GO.DIR.MF.file.annot,MM.GO.IND.BP.file.annot,MM.GO.IND.CC.file.annot,
		MM.GO.IND.MF.file.annot,MM.KEGG.file.annot,MM.locus.name,MM.locus.cyto,SC.GO.DIR.BP.file.annot,SC.GO.DIR.CC.file.annot,SC.GO.DIR.MF.file.annot,
		SC.GO.IND.BP.file.annot,SC.GO.IND.CC.file.annot,SC.GO.IND.MF.file.annot,SC.KEGG.file.annot,SC.locus.name,SC.locus.cyto,
		annot.date,file="sysdata.rda",compress=TRUE) 
	q(save="no")

}


