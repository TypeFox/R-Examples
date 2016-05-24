#' File Fetcher
#'
#' Download Protein Alignment and Accessory Files
#' @param Loci HLA Loci to be fetched. Limited Loci available.
#' @note This function is for internal BIGDAWG use only.
GetFiles <- function(Loci) {
  #downloads *_prot.txt alignment files
  #downloads hla_nom_p.txt file
  #downloads hla.xml file
  
  # Get P-Groups Files
  download.file("http://hla.alleles.org/wmda/hla_nom_p.txt",destfile="hla_nom_p.txt")
  
  # Remove Number from Loci names
  Loci.rep <- unique(gsub("[[:digit:]]","",Loci))
  
  # Get Locus Based Alignments
  for(i in 1:length(Loci.rep)) {
    Locus <- Loci.rep[i]
    FileName <- paste(Locus,"_prot.txt",sep="")
    URL <- paste("ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/alignments/",FileName,sep="")
    download.file(URL,FileName)
    
    if (i==1) {
      tmp <- read.table(FileName,fill=T,sep="\t",stringsAsFactors=F)
      Release <- tmp[2,]
      write.table(Release,file="Release.txt",sep="\t",col.names=F,row.names=F,quote=F)
    }
    
  }
  
}

#' HLA P group File Formatter
#'
#' Format the hla_nom_p.txt read table object for a specific locus.
#' @param x P group object from read.table command.
#' @param Locus Locus to be filtered on.
#' @note This function is for internal BIGDAWG use only.
PgrpFormat <- function(x,Locus) {  
  #Format
  x[,1] <- gsub("\\*","",x[,1])
  x.sub <- x[which(x[,1]==Locus),]
  rownames(x.sub) <- NULL
  x.sub[,2] <- sapply(x.sub[,2],function(i) paste(paste(Locus,"*",unlist(strsplit(i,"/")),sep=""),collapse="/"))
  colnames(x.sub) <- c("Locus","Allele","P.Group")
  
  #Expand
  x.list <- list()
  for(i in 1:nrow(x.sub)) {
    if(grepl("/",x.sub[i,'Allele'],fixed=T)) {
      tmp <- unlist(strsplit(x.sub[i,'Allele'],"/"))
      tmp <- cbind(rep(Locus,length(tmp)),tmp,rep(x.sub[i,'P.Group'],length(tmp)))
      colnames(tmp) <- colnames(x.sub)
      x.list[[i]] <- tmp
    } else {
      x.list[[i]] <- x.sub[i,]
    }
  }
  x.list <- do.call(rbind,x.list)
  x.list <- x.list[order(x.list[,'Allele']),]
  rownames(x.list) <- NULL
  colnames(x.list) <- c("Locus","Allele","P.Group")
  return(x.list)
}

#' HLA P group Finder
#'
#' Identify P group for a given allele if exists.
#' @param x Allele of interest.
#' @param y Formatted P groups.
#' @note This function is for internal BIGDAWG use only.
PgrpExtract <- function(x,y) {
  getRow <- grep(x,y[,'Allele'],fixed=T)
  if(length(getRow)>=1) {
    if(length(getRow)>1) { getRow <- getRow[which(sapply(as.character(y[getRow,'Allele']),nchar)==nchar(x))] }
    return(as.character(y[getRow,'P.Group']))
  } else { return("") }
}

#' Protein Exon Alignment Formatter
#'
#' Dynamically creates an alignmnet of Allele exons for Analysis.
#' @param Locus Locus alignment to be formatted.
#' @param RefTab Reference exon protein information for alignment formatting.
#' @note This function is for internal BIGDAWG use only.
ExonPtnAlign.Create <- function(Locus,RefTab) {
  
  AlignMatrix <- NULL; rm(AlignMatrix)
  
  #Read in P-Groups 
  Pgrps <- read.table("hla_nom_p.txt",fill=T,header=F,sep=";",stringsAsFactors=F,strip.white=T,colClasses="character")
  Pgrps <- PgrpFormat(Pgrps,Locus)
  
  #Read in Alignment
  Locus.sub <- gsub("[[:digit:]]","",Locus)
  Name <- paste(Locus.sub,"_prot.txt",sep="")
  Align <- read.table(Name,fill=T,header=F,sep="\t",stringsAsFactors=F,strip.white=T,colClasses="character")
  
  #Trim
  Align <- as.matrix(Align[-c(1:5),]) #Remove Header
  Align <- as.matrix(Align[-nrow(Align),]) #Remove Footer
  Align <- as.matrix(Align[-grep("\\|",Align[,1]),1]) #Remove Pipes
  
  #Begin Formatting
  Align[,1] <- sapply(Align[,1],FUN=sub,pattern=" ",replacement="~")
  Align[,1] <- sapply(Align[,1],FUN=gsub,pattern=" ",replacement="")
  Align <- strsplit(Align[,1],"~")
  Align <- as.matrix(do.call(rbind,Align))
  
  #Adjust rows where Sequence column == Allele Name
  Align[which(Align[,1]==Align[,2]),2] <- ""
  
  #Find start of repeating blocks
  Start <- which(Align[,1]=="Prot")
  End <- c(Start[2:length(Start)] - 1,nrow(Align))
  
  #Rearrange Imported Alignment, build Alignment Block file
  Block <- Align[Start[1]:End[1],]
  i=2
  Block.size <- End[1]-Start[1]
  while (i <= length(Start) ) {
    if(End[i]-Start[i]==Block.size) { Block <- cbind(Block,Align[Start[i]:End[i],2]) }
    i <- i + 1
  }; rm(i)
  
  #Remove first line with AA and Position Designations
  Block <- Block[-1,]
  
  #Paste Sequence into Single Column
  Block <- cbind(Block[,1],apply(Block[,2:ncol(Block)],MARGIN=1,paste,collapse=""))
  
  #Split Allele name into separate Locus and Allele, Send Back to Align object
  AlignAlleles <- do.call(rbind,strsplit(Block[,1],"[*]"))
  AlignAlleles <- cbind(AlignAlleles,apply(AlignAlleles,MARGIN=c(1,2),FUN=GetField,Res=2)[,2])
  rownames(AlignAlleles) <- NULL
  
  Align <- cbind(AlignAlleles,Block)
  colnames(Align) <- c("Locus","Allele","Trimmed","FullName","Sequence")
  
  #Define Reference
  RefSeq <- Align[1,]
  
  #Ensure Locus Specific Rows - Add Reference
  Align <- Align[which(Align[,'Locus']==Locus),]
  Align <- rbind(RefSeq,Align)
  Align[1,1:4]  <- "RefSeq"
  rownames(Align) <- NULL
  
  #Get Refence Exon Map
  RefExon <- RefTab[which(RefTab[,'Locus']==Locus),'Reference.Peptide']
  RefStart <- as.numeric(RefTab[which(RefTab[,'Locus']==Locus),'Reference.Start'])
  RefAllele <- RefTab[which(RefTab[,'Locus']==Locus),'Reference.Allele']
  
  #Extract Exon Specific Location Based on RefExon
  Estart <- as.numeric(regexpr(RefExon,Align[1,'Sequence']))
  Eend <- (Estart + nchar(RefExon))-1
  Align[,'Sequence'] <- substr(Align[,'Sequence'],Estart,Eend)
  
  #Split Sub Alignment into composite elements
  Align.split <- sapply(Align[,'Sequence'],strsplit,split="*")
  Align.split <- lapply(Align.split,function(x) c(x,rep("-",nchar(RefExon)-length(x))))
  Align.split <- do.call(rbind,Align.split)
  rownames(Align.split) <- NULL
  
  AlignMatrix <- cbind(Align[,1:4],Align.split)
  colnames(AlignMatrix) <- c("Locus","Allele","Trimmed","FullName",paste("Position",seq(RefStart,ncol(AlignMatrix[,5:ncol(AlignMatrix)])+RefStart-1),sep="."))
  
  #Propagate Consensus Positions
  for(i in 5:ncol(AlignMatrix)) {
    x <- AlignMatrix[,i]
    x[which(x=="-")] <- x[1]
    AlignMatrix[,i] <- x
  }
  
  #Remove Reference
  AlignMatrix <- AlignMatrix[-1,]
  
  #Assign P groups
  AlignMatrix <- cbind(AlignMatrix,sapply(AlignMatrix[,'FullName'],PgrpExtract,y=Pgrps))
  colnames(AlignMatrix)[ncol(AlignMatrix)] <- "P group"
  
  #Add Absent Allele (Absence due to lack of allele and not lack of typing information)
  AlignMatrix <- rbind(c(Locus,"00:00:00:00","00:00",paste(Locus,"*00:00:00:00",sep=""),rep("^",ncol(AlignMatrix)-4)),
                       AlignMatrix)
  
  #Tally Unknowns as separate column 
  AlignMatrix <- cbind(AlignMatrix,
                       apply(AlignMatrix[,5:(ncol(AlignMatrix)-1)],MARGIN=1,FUN=function(x) length(which(unlist(grep("*",x,fixed=T))>0))) )
  colnames(AlignMatrix)[ncol(AlignMatrix)] <- "Unknowns"
  
  #Tally Null Positions as separate column
  AlignMatrix <- cbind(AlignMatrix,
                       apply(AlignMatrix[,5:(ncol(AlignMatrix)-2)],MARGIN=1,FUN=function(x) length(which(unlist(grep("-",x,fixed=T))>0))) )
  colnames(AlignMatrix)[ncol(AlignMatrix)] <- "NullPositions"
  
  #Tally InDels
  AlignMatrix <- cbind(AlignMatrix,
                       apply(AlignMatrix[,5:(ncol(AlignMatrix)-3)],MARGIN=1,FUN=function(x) length(which(unlist(grep(".",x,fixed=T))>0))) )
  colnames(AlignMatrix)[ncol(AlignMatrix)] <- "InDels"
  
  rownames(AlignMatrix) <- NULL
  FileName <- paste("ExonPtnAlign_",Locus,".obj",sep="")
  
  save(AlignMatrix,file=FileName)
  
}

#' Alignment Object Creator
#'
#' Synthesize Object for Exon Protein Alignments.
#' @param Loci Loci to be bundled.
#' @param Release IMGT/HLA database release version.
#' @param RefTab Data of reference exons used for protein alignment creation.
#' @note This function is for internal BIGDAWG use only.
AlignObj.Create <- function(Loci,Release,RefTab) {
  
  AlignMatrix <- NULL; rm(AlignMatrix)
  
  ExonPtnList <- list()
  for(i in 1:length(Loci)) {
    Locus <- Loci[i]
    FileName <- paste("ExonPtnAlign_",Locus,".obj",sep="")
    load(FileName) #Loads AlignMatrix
    ExonPtnList[[Locus]] <- AlignMatrix
  }
  ExonPtnList[['Release']] <- Release
  ExonPtnList[['RefExons']] <- RefTab
  save(ExonPtnList,file="ExonPtnAlign.obj")
  
}

#' Updated Alignment Object Creator
#'
#' Synthesize Object for Exon Protein Alignments.
#' @param Loci Loci to be bundled.
#' @param Release IMGT/HLA database release version.
#' @param RefTab Data of reference exons used for protein alignment creation.
#' @note This function is for internal BIGDAWG use only.
AlignObj.Update <- function(Loci,Release,RefTab) {
  
  AlignMatrix <- NULL; rm(AlignMatrix)
  
  UpdatePtnList <- list()
  for(i in 1:length(Loci)) {
    Locus <- Loci[i]
    FileName <- paste("ExonPtnAlign_",Locus,".obj",sep="")
    load(FileName) #Loads AlignMatrix
    UpdatePtnList[[Locus]] <- AlignMatrix
  }
  UpdatePtnList[['Release']] <- Release
  UpdatePtnList[['RefExons']] <- RefTab
  save(UpdatePtnList,file="UpdatePtnAlign.RData")
  
}
