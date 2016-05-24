#library(parallel)
# Load needed library
##library(XML)
##load("geneinfo.rda")
##globalVariables("geneinfo")
#
getTCGA <- function(disease="GBM", data.type="RNASeq2", type="", filter="Y", p=getOption("mc.cores", 2), clinical=FALSE, cvars="OS"){
  
  # Check if it's valid datatype
  data.good <-  c("RNASeq2", "RNASeq", "miRNASeq", "CNA_SNP", "CNV_SNP", "CNA_CGH", 
                  "Methylation", "Mutation", "mRNA_Array", "miRNA_Array")
  if(! (data.type %in% data.good)){
    message("Error: Not recognized datatype for Firehose\n")
    dat <- NULL
    return(dat)
  }
  
  # Parse links
  ldoc <- tryCatch({
    ldoc <- XML::htmlTreeParse("http://gdac.broadinstitute.org/runs/stddata__latest/", useInternalNodes = T)
  } ,error = function(e){
    ldoc = NULL
  })
  
  if(is.null(ldoc)){
    message("Error: Problem connect to Firehose. Please ensure Internet connection is working.\n")
    dat <- NULL
    return(dat)
  }
  
  # parse to get data set and URL for each datasets
  datasets <- XML::xpathSApply(ldoc, "//a[contains(@href, 'Standardized+Data+Run+Release+Notes')]", XML::xmlValue)
  
  if(!(disease %in% datasets)){
    message("Error: Not recognized Disease Abbreviation for Firehose\n")
    dat <- NULL
    return(dat)
  }
  
  if(Sys.getenv("TAR") == "" | Sys.getenv("R_GZIPCMD") ==""){
	message("Error: TAR is not installed in the system. Data unzip failed.\n")
    dat <- NULL
    return(dat)
  }
  
  dataset <- disease
  llinks = unlist(XML::xpathApply(ldoc, "//a[@href]", XML::xmlGetAttr, "href"))
  dlinks = llinks[grepl(paste("/data/", dataset, "/", sep=""), llinks)]
  ddoc = XML::htmlTreeParse(dlinks, useInternalNodes = T)
  
  # -----------------
  # Get data
  if(data.type=="RNASeq2"){
    dat <- RNASeqV2(ddoc=ddoc, dlinks=dlinks, dataset=dataset)
    
    if(is.null(dat)){
      return(dat)
    }
    if(!clinical){
      return(list(dat=dat, clinical=NULL, merged.dat=NULL))
    }
    
    if(clinical){
      gdats <- list()
      cli <- Clinical(ddoc=ddoc, dlinks=dlinks, dataset=dataset)
      mdat <- NULL
      if(!is.null(cli)){
        # combine clinical and data together
        mdat <- MatrixMerge(dat=dat, cli=cli, cvars=cvars)
      }
      return(list(dat=dat, clinical=cli, merged.dat=mdat))
    }
  }
  
  # ----------
  if(data.type=="RNASeq"){
    if(type==""){
      dat <- RNASeq(ddoc=ddoc, dlinks=dlinks, dataset=dataset)
    }
    else{
      dat <- RNASeq(ddoc=ddoc, dlinks=dlinks, dataset=dataset, type=type)
    }
    
    if(is.null(dat)){
      return(dat)
    }
    if(!clinical){
      return(list(dat=dat, clinical=NULL, merged.dat=NULL))
    }
    
    if(clinical){
      gdats <- list()
      cli <- Clinical(ddoc=ddoc, dlinks=dlinks, dataset=dataset)
      mdat <- NULL
      if(!is.null(cli)){
        # combine clinical and data together
        mdat <- MatrixMerge(dat=dat, cli=cli, cvars=cvars)
      }
      return(list(dat=dat, clinical=cli, merged.dat=mdat))
    }
  }
  # ----------
  if(data.type=="miRNASeq"){
    if(type==""){
      dat <- miRNASeq(ddoc=ddoc, dlinks=dlinks, dataset=dataset)
    }
    else{
      dat <- miRNASeq(ddoc=ddoc, dlinks=dlinks, dataset=dataset, type=type)
    }
    if(is.null(dat)){
      return(dat)
    }
    if(!clinical){
      return(list(dat=dat, clinical=NULL, merged.dat=NULL))
    }
    
    if(clinical){
      gdats <- list()
      cli <- Clinical(ddoc=ddoc, dlinks=dlinks, dataset=dataset)
      mdat <- NULL
      if(!is.null(cli)){
        # combine clinical and data together
        mdat <- MatrixMerge(dat=dat, cli=cli, cvars=cvars)
      }
      return(list(dat=dat, clinical=cli, merged.dat=mdat))
    }
  }
  # ----------
  if(data.type=="CNA_SNP"){
    dat <- CNA_SNP(ddoc=ddoc, dlinks=dlinks, dataset=dataset, filter=filter)
    
    if(is.null(dat)){
      return(dat)
    }
    if(!clinical){
      return(list(dat=dat, clinical=NULL, merged.dat=NULL))
    }
    
    if(clinical){
      gdats <- list()
      cli <- Clinical(ddoc=ddoc, dlinks=dlinks, dataset=dataset)
      mdat = NULL
      if(!is.null(cli)){
        # combine clinical and data together
        mdat <- MatrixMerge(dat=dat, cli=cli, cvars=cvars)
      }
      return(list(dat=dat, clinical=cli, merged.dat=mdat))
    }
  }
  # ----------
  if(data.type=="CNV_SNP"){
    dat <- CNV_SNP(ddoc=ddoc, dlinks=dlinks, dataset=dataset, filter=filter)
    
    if(is.null(dat)){
      return(dat)
    }
    if(!clinical){
      return(list(dat=dat, clinical=NULL, merged.dat=NULL))
    }
    
    if(clinical){
      gdats <- list()
      cli <- Clinical(ddoc=ddoc, dlinks=dlinks, dataset=dataset)
      mdat = NULL
      if(!is.null(cli)){
        # combine clinical and data together
        mdat <- MatrixMerge(dat=dat, cli=cli, cvars=cvars)
      }
      return(list(dat=dat, clinical=cli, merged.dat=mdat))
    }
  }
  # ----------
  if(data.type=="CNA_Seq"){
    dat <- CNA_Seq(ddoc=ddoc, dlinks=dlinks, dataset=dataset)
    
    if(is.null(dat)){
      return(dat)
    }
    if(!clinical){
      return(dat)
    }
    
    if(clinical){
      gdats <- list()
      cli <- Clinical(ddoc=ddoc, dlinks=dlinks, dataset=dataset)
      
      if(!is.null(cli)){
        # combine clinical and data together
        mdat <- MatrixMerge(dat=dat, cli=cli, cvars=cvars)
      }
      return(list(dat=dat, clinical=cli, merged.dat=mdat))
    }
  }
  # ----------
  if(data.type=="CNA_CGH"){
    if(type==""){
      dat <- CNA_CGH(ddoc=ddoc, dlinks=dlinks, dataset=dataset, filter=filter)
    }
    else{
      dat <- CNA_CGH(ddoc=ddoc, dlinks=dlinks, dataset=dataset, type=type, filter=filter)
    }
        
    if(is.null(dat)){
      return(dat)
    }
    if(!clinical){
      return(list(dat=dat, clinical=NULL, merged.dat=NULL))
    }
    
    if(clinical){
      gdats <- list()
      cli <- Clinical(ddoc=ddoc, dlinks=dlinks, dataset=dataset)
      mdat = NULL
      if(!is.null(cli)){
        # combine clinical and data together
        mdat <- MatrixMerge(dat=dat, cli=cli, cvars=cvars)
      }
      return(list(dat=dat, clinical=cli, merged.dat=mdat))
    }
  }
  # ----------
  if(data.type=="Methylation"){
    
    if(type==""){
      dat <- Methylation(ddoc=ddoc, dlinks=dlinks, dataset=dataset, p=p)
    }
    else{
      dat <- Methylation(ddoc=ddoc, dlinks=dlinks, dataset=dataset, p=p, type=type)
    }
    
    if(is.null(dat)){
      return(dat)
    }
    if(!clinical){
      return(list(dat=as.matrix(dat[, -c(1,2,3)]), cpgs=dat[,1:3], clinical=NULL, merged.dat=NULL))
    }
    
    if(clinical){
      gdats <- list()
      cli <- Clinical(ddoc=ddoc, dlinks=dlinks, dataset=dataset)
      mdat <- NULL
      if(!is.null(cli)){
        # combine clinical and data together
        mdat <- MatrixMerge(dat=dat, cli=cli, cvars=cvars)
      }
      return(list(dat=as.matrix(dat[, -c(1,2,3)]), cpgs=dat[,1:3], clinical=cli, merged.dat= mdat))
    }
  }
  # ----------
  if(data.type=="Mutation"){
    if(type==""){
      dat <- Mutation(ddoc=ddoc, dlinks=dlinks, dataset=dataset)
    }
    else{
      dat <- Mutation(ddoc=ddoc, dlinks=dlinks, dataset=dataset, type=type)
    }
    if(is.null(dat)){
      return(dat)
    }
    if(!clinical){
      return(list(dat=dat, clinical=NULL, merged.dat=NULL))
    }
    
    if(clinical){
      gdats <- list()
      cli <- Clinical(ddoc=ddoc, dlinks=dlinks, dataset=dataset)
      mdat <- NULL
      
      if(!is.null(cli)){
        # combine clinical and data together
        mdat <- MatrixMerge(dat=dat, cli=cli, cvars=cvars)
      }
      return(list(dat=dat, clinical=cli, merged.dat=mdat))
    }
  }
  # ----------
  if(data.type=="mRNA_Array"){
    if(type==""){
      dat <- mRNA_Array(ddoc=ddoc, dlinks=dlinks, dataset=dataset)
    }
    else{
      dat <- mRNA_Array(ddoc=ddoc, dlinks=dlinks, dataset=dataset, type=type)
    }
    
    if(is.null(dat)){
      return(dat)
    }
    if(!clinical){
      return(list(dat=dat, clinical=NULL, merged.dat=NULL))
    }
    
    if(clinical){
      gdats <- list()
      cli <- Clinical(ddoc=ddoc, dlinks=dlinks, dataset=dataset)
      mdat <- NULL
      if(!is.null(cli)){
        # combine clinical and data together
        mdat <- MatrixMerge(dat=dat, cli=cli, cvars=cvars)
      }
      return(list(dat=dat, clinical=cli, merged.dat=mdat))
    }
  }
  # ----------
  if(data.type=="miRNA_Array"){
    dat <- miRNA_Array(ddoc=ddoc, dlinks=dlinks, dataset=dataset)
    
    if(is.null(dat)){
      return(dat)
    }
    if(!clinical){
      return(list(dat=as.matrix(dat), clinical=NULL, merged.dat=NULL))
    }
    
    if(clinical){
      gdats <- list()
      cli <- Clinical(ddoc=ddoc, dlinks=dlinks, dataset=dataset)
      mdat <- NULL
      if(is.null(cli)){
        message("Error: clinical data not available.")
      }
      if(!is.null(cli)){
        # combine clinical and data together
        mdat <- MatrixMerge(dat=dat, cli=cli, cvars=cvars)
      }
      return(list(dat=as.matrix(dat), clinical=cli, merged.dat=as.matrix(mdat)))
    }
  }
}
#
#

RNASeqV2 <- function(ddoc, dlinks, dataset){
  # get RNAseq_Gene_v2
  keyWord = paste("","Level_3__RSEM_genes_normalized__data.Level_3",sep="")
  keyWord = paste("//a[contains(@href, '",keyWord,"')]",sep="")
  plinks = XML::xpathSApply(ddoc, keyWord, XML::xmlGetAttr, 'href')
  plinks = plinks[grepl(paste("*", dataset, ".Merge_rnaseqv2__illuminahiseq*._rnaseqv2__.*.tar[.]gz$", sep=""),plinks)]
  
  if(length(plinks) == 0){
    message("Error: No data available for download. Please ensure the data is available from TCGA. \n")
    dat <- NULL
    return(dat)
  }
  
  # get data processing time-stamp
  timestamp <- unlist(strsplit(dlinks, "/"))
  timestamp <- timestamp[length(timestamp)]
  
  gdats <- list()
  for(i in 1:length(plinks)){
    # Get the right download links
    download_link = paste(dlinks,trim(plinks[i]),sep="/")
    
    message("RNASeqV2 data will be imported! This may take some time!")
    
    # Download data in tar format and then extract the right text file
    utils::download.file(url=download_link,destfile=paste(dataset,"-RNAseq2GeneNorm.tar.gz",sep=""),method="auto",quiet = TRUE, mode = "w")
    fileList <- utils::untar(paste(dataset,"-RNAseq2GeneNorm.tar.gz",sep=""),list=TRUE)
    grepSearch = paste("*.",dataset,"[.]rnaseqv2__.*.__Level_3__RSEM_genes_normalized__data.data.txt$",sep="")
    fileList = fileList[grepl(grepSearch,fileList)]
    utils::untar(paste(dataset,"-RNAseq2GeneNorm.tar.gz",sep=""),files=fileList)
    
    # Rename and clean up
    fname = paste(dataset,"_", timestamp, "-RNAseq2GeneNorm.txt",sep="")
    file.rename(from=fileList,to=fname)
    file.remove(paste(dataset,"-RNAseq2GeneNorm.tar.gz",sep=""))
    delFodler <- paste(getwd(),"/",strsplit(fileList,"/")[[1]][1],sep="")
    ###message(delFodler)
    unlink(delFodler, recursive = TRUE)
    
    #Get selected type only
    tmpCols = utils::read.delim(fname,nrows=1,colClasses="character")
    tmpdat = utils::read.delim(fname, skip=1, sep="\t", stringsAsFactors=F)
    
    colOrder <- 1:ncol(tmpCols)
    colOrder <- colOrder[tmpCols[1,] == "normalized_count"]
    
    #	badgene <- grepl("\\?\\|",tmpdat[,1])
    gnames <- sapply(tmpdat[,1], function(s) unlist(strsplit(s, "\\|"))[1])
    badg = which(gnames == "?" | duplicated(gnames))
    
    gdat <- tmpdat[-badg, colOrder]
    colnames(gdat) <- colnames(tmpCols)[colOrder]
    colnames(gdat) <- gsub("\\.", "-", colnames(gdat))
    rownames(gdat) <- gnames[-badg]
    
    message(paste(nrow(gdat) ,"genes have been imported!"))
    gdats[[i]] <- as.matrix(gdat)
    file.remove(fname)
  }
  
  if(length(gdats) == 1){
    gdats <- gdats[[1]]
  }
  
  return(gdats)
}
#
# type could be "count" and "RPKM"
RNASeq <- function(ddoc, dlinks, dataset, type="count"){
  
  if(! (type %in% c("count", "RPKM"))){
    message("Error: Invalid type.")
    gdat <- NULL
    return(gdat)
  }
  
  if(type=="count"){
    type = "raw_counts"  
  }
  
  keyWord = paste("","Level_3__gene_expression__data.Level_3",sep="")
  keyWord = paste("//a[contains(@href, '",keyWord,"')]",sep="")
  plinks = XML::xpathSApply(ddoc, keyWord, XML::xmlGetAttr, 'href')
  plinks = plinks[grepl("*.Merge_rnaseq__.*._rnaseq__.*.tar[.]gz$",plinks)]
  
  if(length(plinks) == 0){
    message("Error: No data available for download. Please ensure the data is available from TCGA. \n")
    dat <- NULL
    return(dat)
  }
  
  if(length(plinks) > 1){
    plinks = plinks[grepl("*_illuminahiseq_*", plinks)]
  }
  
  # get data processing time-stamp
  timestamp <- unlist(strsplit(dlinks, "/"))
  timestamp <- timestamp[length(timestamp)]
  
  gdats <- list()
  for(i in 1:length(plinks)){
    # Get the right download links
    download_link = paste(dlinks,trim(plinks[i]),sep="/")
    
    message("RNAseq data will be imported! This may take some time!")
    
    # Download data in tar format and then extract the right text file
    utils::download.file(url=download_link,destfile=paste(dataset,"-RNAseqGene.tar.gz",sep=""),method="auto",quiet = TRUE, mode = "w")
    fileList <- utils::untar(paste(dataset,"-RNAseqGene.tar.gz",sep=""),list=TRUE)
    grepSearch = paste("*.",dataset,"[.]rnaseq__.*.__Level_3__gene_expression__data.data.txt$",sep="")
    fileList = fileList[grepl(grepSearch,fileList)]
    utils::untar(paste(dataset,"-RNAseqGene.tar.gz",sep=""),files=fileList)
    
    # Rename and clean up
    fname = paste(dataset,"_", timestamp, "-RNAseqGene.txt",sep="")
    
    file.rename(from=fileList,to=fname)
    file.remove(paste(dataset,"-RNAseqGene.tar.gz",sep=""))
    delFodler <- paste(getwd(),"/",strsplit(fileList,"/")[[1]][1],sep="")
    unlink(delFodler, recursive = TRUE)
    
    #Get selected type only
    tmpCols = utils::read.delim(fname,nrows=1,colClasses="character")
    tmpdat = utils::read.delim(fname, skip=1, sep="\t", stringsAsFactors=F)
    
    colOrder <- 1:ncol(tmpCols)
    colOrder <- colOrder[tmpCols[1,] == type]
    
    gnames <- sapply(tmpdat[,1], function(s) unlist(strsplit(s, "\\|"))[1])
    badg = which(gnames == "?" | duplicated(gnames))
    
    gdat <- tmpdat[-badg, colOrder]
    colnames(gdat) <- colnames(tmpCols)[colOrder]
    colnames(gdat) <- gsub("\\.", "-", colnames(gdat))
    rownames(gdat) <- gnames[-badg]
    
    message(paste(nrow(gdat) ,"genes have been imported!"))
    gdats[[i]] <- as.matrix(gdat)
    
    file.remove(fname)
  }
  
  if(length(gdats) == 1){
    gdats <- gdats[[1]]
  }
  
  return(gdats)
}
#
# type could be "count" and "rpmmm"
miRNASeq <- function(ddoc, dlinks, dataset, type="count"){
  
  if(! (type %in% c("count", "rpmmm"))){
    message("Error: Invalid type.")
    gdat <- NULL
    return(gdat)
  }
  
  if(type == "count"){
    type = "read_count"
  }
  
  if(type == "rpmmm"){
    type = "reads_per_million_miRNA_mapped"
  }
  
  keyWord = paste("","Level_3__miR_gene_expression__data.Level_3",sep="")
  keyWord = paste("//a[contains(@href, '",keyWord,"')]",sep="")
  plinks = XML::xpathSApply(ddoc, keyWord, XML::xmlGetAttr, 'href')
  plinks = plinks[grepl(paste("*.",dataset,"[.]Merge_mirnaseq__.*.hiseq_mirnaseq__.*.tar[.]gz$",sep=""),plinks)]
  
  if(length(plinks) == 0){
    message("Error: No data available for download. Please ensure the data is available from TCGA. \n")
    dat <- NULL
    return(dat)
  }
  
  # get data processing time-stamp
  timestamp <- unlist(strsplit(dlinks, "/"))
  timestamp <- timestamp[length(timestamp)]
  
  gdats <- list()

  for(i in 1:length(plinks)){
    # Get the right download links
    download_link = paste(dlinks,trim(plinks[i]),sep="/")
    
    message("miRNAseq data will be imported! This may take some time!")
    
    # Download data in tar format and then extract the right text file
    utils::download.file(url=download_link,destfile=paste(dataset,"-miRNAseqGene.tar.gz",sep=""),method="auto",quiet = TRUE, mode = "w")
    fileList <- utils::untar(paste(dataset,"-miRNAseqGene.tar.gz",sep=""),list=TRUE)
    grepSearch = paste("*.",dataset,"[.]mirnaseq__.*.__Level_3__miR_gene_expression__data.data.txt$",sep="")
    fileList = fileList[grepl(grepSearch,fileList)]
    utils::untar(paste(dataset,"-miRNAseqGene.tar.gz",sep=""),files=fileList)
    
    # Rename and clean up
    fname = paste(dataset,"_", timestamp, "-miRNAseqGene.txt",sep="")
    file.rename(from=fileList,to=fname)
    file.remove(paste(dataset,"-miRNAseqGene.tar.gz",sep=""))
    delFodler <- paste(getwd(),"/",strsplit(fileList,"/")[[1]][1],sep="")
    ##message(delFodler)
    unlink(delFodler, recursive = TRUE)
    
    #Get selected type only
    
    tmpCols = utils::read.delim(fname,nrows=1,colClasses="character")
    tmpdat = utils::read.delim(fname, skip=1, sep="\t", stringsAsFactors=F)
    
    colOrder <- 1:ncol(tmpCols)
    colOrder <- colOrder[tmpCols[1,] == type]
    
    gnames <- tmpdat[,1]
    badg <- which(duplicated(gnames))
    
    if(length(badg) > 0){
      gdat <- tmpdat[-badg, colOrder]
      colnames(gdat) <- colnames(tmpCols)[colOrder]
      colnames(gdat) <- gsub("\\.", "-", colnames(gdat))
      rownames(gdat) <- gnames[-badg]
    }
    if(length(badg) == 0) {
      gdat <- tmpdat[, colOrder]
      colnames(gdat) <- colnames(tmpCols)[colOrder]
      colnames(gdat) <- gsub("\\.", "-", colnames(gdat))
      rownames(gdat) <- gnames
    }
    message(paste(nrow(gdat) ,"genes have been imported!"))
    
    gdats[[i]] <- as.matrix(gdat)
    
    file.remove(fname)
  }
  
  if(length(gdats) == 1){
    gdats <- gdats[[1]]
  }
  
  return(gdats)
}
#
CNA_SNP <- function(ddoc, dlinks, dataset,filter="Y"){
  keyWord = paste("","Level_3__segmented_scna_hg19__seg.Level_3",sep="")
  keyWord = paste("//a[contains(@href, '",keyWord,"')]",sep="")
  plinks = XML::xpathSApply(ddoc, keyWord, XML::xmlGetAttr, 'href')
  plinks = plinks[grepl(paste("*.",dataset,"[.]Merge_snp__.*.__Level_3__segmented_scna_hg19__seg.Level_3.*.tar[.]gz$",sep=""),plinks)]
  
  if(length(plinks) == 0){
    message("Error: No data available for download. Please ensure the data is available from TCGA. \n")
    dat <- NULL
    return(dat)
  }
  
  # get data processing time-stamp
  timestamp <- unlist(strsplit(dlinks, "/"))
  timestamp <- timestamp[length(timestamp)]
  
  gdats <- list()
  for(i in 1:length(plinks)){
    # Get the right download links
    download_link = paste(dlinks,trim(plinks[i]),sep="/")
    
    message("CNA_SNP data will be imported! This may take some time!")
    
    # Download data in tar format and then extract the right text file
    utils::download.file(url=download_link,destfile=paste(dataset,"-CNASNPHg19.tar.gz",sep=""),method="auto",quiet = TRUE, mode = "w")
    fileList <- utils::untar(paste(dataset,"-CNASNPHg19.tar.gz",sep=""),list=TRUE)
    grepSearch = paste("*.",dataset,"[.]snp__.*.__Level_3__segmented_scna_hg19__seg.seg.txt$",sep="")
    fileList = fileList[grepl(grepSearch,fileList)]
    utils::untar(paste(dataset,"-CNASNPHg19.tar.gz",sep=""),files=fileList)
    
    # Rename and clean up
    fname = paste(dataset,"_", timestamp, "-CNASNPHg19.txt",sep="")
    file.rename(from=fileList,to=fname)
    file.remove(paste(dataset,"-CNASNPHg19.tar.gz",sep=""))
    delFodler <- paste(getwd(),"/",strsplit(fileList,"/")[[1]][1],sep="")
    ##message(delFodler)
    unlink(delFodler, recursive = TRUE)
    
    gdat <- utils::read.delim(fname, sep="\t", stringsAsFactors=F)
    
    gdats[[i]] <- gdat
    
    file.remove(fname)
  }
  
  # make the matrix and segments
  gsegs <- gdats
  gdats <- list()
  for(i in 1:length(gsegs)){
    tmat <- gsegs[[i]]
    colnames(tmat) <- c("ID", "chrom", "loc.start", "loc.end", "num.mark", "seg.mean")

    # --
    GeneMatrix <- GetSegments.ByGene(tmat, geneInfof="", filter=filter)
    # Check and rename genes sharing same name in both X and Y chromosomes
    gene.symbols <- GeneMatrix[,5]
    dup.genes <- gene.symbols[duplicated(gene.symbols)]
    rchr <- ifelse(gene.symbols %in% dup.genes, paste("_", GeneMatrix[,1], sep=""), "")
    gene.symbols <- paste(gene.symbols, rchr, sep="")
    
    gdat <- GeneMatrix[,6:ncol(GeneMatrix)]
    rownames(gdat) <- gene.symbols
    # --
    
    gdats[[i]] <- gdat
  }
  
  if(length(gdats) == 1){
    gdats <- gdats[[1]]
    gsegs <- gsegs[[1]]
  }
  
  return(gdats)
}
#
CNV_SNP <- function(ddoc, dlinks, dataset, filter="Y"){
  keyWord = paste("","Level_3__segmented_scna_minus_germline_cnv_hg19__seg.Level_3",sep="")
  keyWord = paste("//a[contains(@href, '",keyWord,"')]",sep="")
  plinks = XML::xpathSApply(ddoc, keyWord, XML::xmlGetAttr, 'href')
  plinks = plinks[grepl(paste("*.",dataset,"[.]Merge_snp__.*.Level_3__segmented_scna_minus_germline_cnv_hg19__seg.Level_3.*.tar[.]gz$",sep=""),plinks)]
  
  if(length(plinks) == 0){
    message("Error: No data available for download. Please ensure the data is available from TCGA. \n")
    dat <- NULL
    return(dat)
  }
  
  # get data processing time-stamp
  timestamp <- unlist(strsplit(dlinks, "/"))
  timestamp <- timestamp[length(timestamp)]
  
  gdats <- list()

  for(i in 1:length(plinks)){
    # Get the right download links
    download_link = paste(dlinks,trim(plinks[i]),sep="/")
    
    message("CNV_SNP data will be imported! This may take some time!")
    
    # Download data in tar format and then extract the right text file
    utils::download.file(url=download_link,destfile=paste(dataset,"-CNVSNPHg19.tar.gz",sep=""),method="auto",quiet = TRUE, mode = "w")
    fileList <- utils::untar(paste(dataset,"-CNVSNPHg19.tar.gz",sep=""),list=TRUE)
    grepSearch = paste("*.",dataset,"[.]snp__.*.Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt$",sep="")
    fileList = fileList[grepl(grepSearch,fileList)]
    utils::untar(paste(dataset,"-CNVSNPHg19.tar.gz",sep=""),files=fileList)
    
    # Rename and clean up
    fname = paste(dataset,"_", timestamp, "-CNVSNPHg19.txt",sep="")
    file.rename(from=fileList,to=fname)
    file.remove(paste(dataset,"-CNVSNPHg19.tar.gz",sep=""))
    delFodler <- paste(getwd(),"/",strsplit(fileList,"/")[[1]][1],sep="")
    ##message(delFodler)
    unlink(delFodler, recursive = TRUE)
    
    gdat <- utils::read.delim(fname, sep="\t", stringsAsFactors=F)
    
    gdats[[i]] <- gdat
    
    file.remove(fname)
  }
  
  # make the matrix and segments
  gsegs <- gdats
  gdats <- list()
  for(i in 1:length(gsegs)){
    tmat <- gsegs[[i]]
    colnames(tmat) <- c("ID", "chrom", "loc.start", "loc.end", "num.mark", "seg.mean")
    # --
    GeneMatrix <- GetSegments.ByGene(tmat, geneInfof="", filter=filter)
    # Check and rename genes sharing same name in both X and Y chromosomes
    gene.symbols <- GeneMatrix[,5]
    dup.genes <- gene.symbols[duplicated(gene.symbols)]
    rchr <- ifelse(gene.symbols %in% dup.genes, paste("_", GeneMatrix[,1], sep=""), "")
    gene.symbols <- paste(gene.symbols, rchr, sep="")
    
    gdat <- GeneMatrix[,6:ncol(GeneMatrix)]
    rownames(gdat) <- gene.symbols
    # --
    
    gdats[[i]] <- gdat
  }
  
  
  if(length(gdats) == 1){
    gdats <- gdats[[1]]
    gsegs <- gsegs[[1]]
  }

  return(gdats)
}
#
CNA_Seq <- function(ddoc, dlinks, dataset){
  keyWord = paste("","__Level_3__segmentation__seg.Level_3",sep="")
  keyWord = paste("//a[contains(@href, '",keyWord,"')]",sep="")
  plinks = XML::xpathSApply(ddoc, keyWord, XML::xmlGetAttr, 'href')
  plinks = plinks[grepl(paste("*.",dataset,"[.]Merge_cna__.*.dnaseq.*.__Level_3__segmentation__seg.Level_3.*.tar[.]gz$",sep=""),plinks)]
  
  if(length(plinks) == 0){
    message("Error: No data available for download. Please ensure the data is available from TCGA. \n")
    dat <- NULL
    return(dat)
  }
  
  # get data processing time-stamp
  timestamp <- unlist(strsplit(dlinks, "/"))
  timestamp <- timestamp[length(timestamp)]
  
  gdats <- list()

  for(i in 1:length(plinks)){
    # Get the right download links
    download_link = paste(dlinks,trim(plinks[i]),sep="/")
    
    message("CNA_Seq data will be imported! This may take some time!")
    
    # Download data in tar format and then extract the right text file
    utils::download.file(url=download_link,destfile=paste(dataset,"-CNAseq.tar.gz",sep=""),method="auto",quiet = TRUE, mode = "w")
    fileList <- utils::untar(paste(dataset,"-CNAseq.tar.gz",sep=""),list=TRUE)
    grepSearch = paste("*.",dataset,"[.]cna__.*.__Level_3__segmentation__seg.seg.txt$",sep="")
    fileList = fileList[grepl(grepSearch,fileList)]
    utils::untar(paste(dataset,"-CNAseq.tar.gz",sep=""),files=fileList)
    
    # Rename and clean up
    fname = paste(dataset,"_", timestamp, "-CNAseq.txt",sep="")
    file.rename(from=fileList,to=fname)
    file.remove(paste(dataset,"-CNAseq.tar.gz",sep=""))
    delFodler <- paste(getwd(),"/",strsplit(fileList,"/")[[1]][1],sep="")
    ##message(delFodler)
    unlink(delFodler, recursive = TRUE)
    
    gdat <- utils::read.delim(fname, sep="\t", stringsAsFactors=F)
    
    gdats[[i]] <- gdat
    
    file.remove(fname)
  }
  
  if(length(gdats) == 1){
    gdats <- gdats[[1]]
  }
  
  return(gdats)
}
#
CNA_CGH <- function(ddoc, dlinks, dataset, type="415K", filter="Y"){
  
  if(! (type %in% c("415K", "244A"))){
    message("Error: Invalid type for CGH data.")
    gdat <- NULL
    return(gdat)
  }
  
  keyWord = paste("","__Level_3__segmentation__seg.Level_3",sep="")
  keyWord = paste("//a[contains(@href, '",keyWord,"')]",sep="")
  plinks = XML::xpathSApply(ddoc, keyWord, XML::xmlGetAttr, 'href')
  
  if(type=="415K"){
    plinks = plinks[grepl(paste("*.",dataset,"[.]Merge_cna__.*.cgh_415k.*.__Level_3__segmentation__seg.Level_3.*.tar[.]gz$",sep=""),plinks)]
  }

  if(type=="244A"){
    plinks = plinks[grepl(paste("*.",dataset,"[.]Merge_cna__.*.cgh_244a.*.__Level_3__segmentation__seg.Level_3.*.tar[.]gz$",sep=""),plinks)]  
  }
    
  if(length(plinks) == 0){
    message("Error: No data available for download. Please make sure the data is available from TCGA. \n")
    dat <- NULL
    return(dat)
  }
  
  # get data processing time-stamp
  timestamp <- unlist(strsplit(dlinks, "/"))
  timestamp <- timestamp[length(timestamp)]
  
  gdats <- list()
  for(i in 1:length(plinks)){
    # Get the right download links
    download_link = paste(dlinks,trim(plinks[i]),sep="/")
    
    if(i==1){
      message("CNA_CGH data will be imported! This may take some time!")
    }
    
    # Download data in tar format and then extract the right text file
    utils::download.file(url=download_link,destfile=paste(dataset,"-CNACGH.tar.gz",sep=""),method="auto",quiet = TRUE, mode = "w")
    fileList <- utils::untar(paste(dataset,"-CNACGH.tar.gz",sep=""),list=TRUE)
    grepSearch = paste("*.",dataset,"[.]cna__.*.__Level_3__segmentation__seg.seg.txt$",sep="")
    fileList = fileList[grepl(grepSearch,fileList)]
    utils::untar(paste(dataset,"-CNACGH.tar.gz",sep=""),files=fileList)
    
    # Rename and clean up
    fname = paste(dataset,"_", timestamp, "-CNACGH.txt",sep="")
    file.rename(from=fileList,to=fname)
    file.remove(paste(dataset,"-CNACGH.tar.gz",sep=""))
    delFodler <- paste(getwd(),"/",strsplit(fileList,"/")[[1]][1],sep="")
    ##message(delFodler)
    unlink(delFodler, recursive = TRUE)
    
    temp <- utils::read.delim(fname, sep="\t", stringsAsFactors=F)
    gdats[[i]] <- temp
    
    file.remove(fname)
  }
  
  # make the matrix and segments
  gsegs <- gdats
  gdats <- list()
  for(i in 1:length(gsegs)){
    tmat <- gsegs[[i]]
    colnames(tmat) <- c("ID", "chrom", "loc.start", "loc.end", "num.mark", "seg.mean")
    
    # --
    GeneMatrix <- GetSegments.ByGene(tmat, geneInfof="", filter=filter)
    # Check and rename genes sharing same name in both X and Y chromosomes
    gene.symbols <- GeneMatrix[,5]
    dup.genes <- gene.symbols[duplicated(gene.symbols)]
    rchr <- ifelse(gene.symbols %in% dup.genes, paste("_", GeneMatrix[,1], sep=""), "")
    gene.symbols <- paste(gene.symbols, rchr, sep="")
    
    gdat <- GeneMatrix[,6:ncol(GeneMatrix)]
    rownames(gdat) <- gene.symbols
    
    # --  
  
    gdats[[i]] <- gdat
  }
  
  
  if(length(gdats) == 1){
    gdats <- gdats[[1]]
    gsegs <- gsegs[[1]]
  }
  
  return(gdats)
}
#
# Create a matrix of genes x samples of CN-change
GetSegments.ByGene <- function(dat, geneInfof="", filter="Y"){
  geneinfo=NULL
  utils::data("geneinfo", envir = environment())
  if(sum(!(filter == "" | is.null(filter)))>0){
    dat = dat[!(dat[,"chrom"] %in% filter), ]  
  }

  seg = CNTools::CNSeg(dat)
  mat= CNTools::getRS(seg, by="gene", imput=FALSE, XY=TRUE, geneMap=geneinfo, what="median")
  matrs = mat@rs
  
  # Remove out the genes if filtered out above, if necessary
  if(sum(!(filter == "" | is.null(filter)))>0){
    matrs = matrs[!(matrs[,"chrom"] %in% filter), ]
  }
  return(matrs)
}
#
Methylation <- function(ddoc, dlinks, dataset, p=getOption("mc.cores", 2L), type="27K"){
  
  if(! (type %in% c("27K", "450K"))){
    message("Error: Invalid type for methylation data.")
    gdat <- NULL
    return(gdat)
  }
  
  keyWord = paste("","__Level_3__within_bioassay_data_set_function__data.Level_3",sep="")
  keyWord = paste("//a[contains(@href, '",keyWord,"')]",sep="")
  plinks = XML::xpathSApply(ddoc, keyWord, XML::xmlGetAttr, 'href')
  
  if(type=="27K"){
    plinks = plinks[grepl(paste("*.",dataset,"[.]Merge_methylation__.*.methylation27.*.__Level_3__within_bioassay_data_set_function__data.Level_3.*.tar[.]gz$",sep=""),plinks)]
  }
  
  if(type=="450K"){
    plinks = plinks[grepl(paste("*.",dataset,"[.]Merge_methylation__.*.methylation450.*.__Level_3__within_bioassay_data_set_function__data.Level_3.*.tar[.]gz$",sep=""),plinks)]
  }
  
  if(length(plinks) == 0){
    message("Error: No data available for download. Please ensure the data is available from TCGA. \n")
    dat <- NULL
    return(dat)
  }
  
  # get data processing time-stamp
  timestamp <- unlist(strsplit(dlinks, "/"))
  timestamp <- timestamp[length(timestamp)]
  
  gdats <- list()
  for(i in 1:length(plinks)){
    # Get the right download links
    download_link = paste(dlinks,trim(plinks[i]),sep="/")
    
    if(i==1){
      message("Methylation data will be imported! This may take some time!")
    }
    
    # Download data in tar format and then extract the right text file
    utils::download.file(url=download_link,destfile=paste(dataset,"-Methylation.tar.gz",sep=""),method="auto",quiet = TRUE, mode = "w")
    fileList <- utils::untar(paste(dataset,"-Methylation.tar.gz",sep=""),list=TRUE)
    grepSearch = paste("*.",dataset,"[.]methylation__.*.__Level_3__within_bioassay_data_set_function__data.data.txt$",sep="")
    fileList = fileList[grepl(grepSearch,fileList)]
    utils::untar(paste(dataset,"-Methylation.tar.gz",sep=""),files=fileList)
    
    # Rename and clean up
    fname = paste(dataset,"_", timestamp, "-Methylation.txt",sep="")
    file.rename(from=fileList,to=fname)
    file.remove(paste(dataset,"-Methylation.tar.gz",sep=""))
    delFodler <- paste(getwd(),"/",strsplit(fileList,"/")[[1]][1],sep="")
    ## message(delFodler)
    unlink(delFodler, recursive = TRUE)
    
    #Get selected type only
    tmpCols = utils::read.delim(fname,nrows=1,colClasses="character")
    colOrder <- 1:ncol(tmpCols)
    colOrder <- colOrder[tmpCols[1,] == "Beta_value"]
    colOrder <- c(1,3,4,5, colOrder)
    
    readf = 1
    readl = 2
    block <- 5000*p
    #		block <- 5000
    gdat <- c()
    while(readf){
      
      sblock = block/p
      skips = c(readl, readl+(1:(p-1)*sblock))
      
      wrapper <- function(j) {
        tt = NULL
        tryCatch({ tt = utils::read.delim(fname, skip=j, nrows=sblock, stringsAsFactors=F, header=F)
        }
        ,error = function(e){tt = NULL}
        )
        return(tt)
      }
      temps <- parallel::mclapply(skips, wrapper, mc.cores=p)
      
      temp <- c()
      for(j in 1:p){
        if(!is.null(temps[[j]])){
          temp <- rbind(temp, temps[[j]])
        }
      }
      
      if(nrow(temp) > 0){
        gdat <- rbind(gdat, temp[,colOrder])
        readl = readl + nrow(temp)
      }
      if(nrow(temp) == 0 | nrow(temp) < block){
        readf = 0
      }
    }
    colnames(gdat) <- gsub("\\.", "-", c(tmpCols[colOrder[1:4]], names(tmpCols)[colOrder[5:length(colOrder)]]))
    rownames(gdat) <- gdat[,1]
    gdat <- gdat[,-1]
 
    # keep only cpg probes
    gdat <- gdat[grep("^cg", rownames(gdat)), ]
    
    message(paste(nrow(gdat) ,"CPG probes have been imported!"))
    
    gdats[[i]] <- gdat
    
    file.remove(fname)
  }
  
  if(length(gdats) ==1){
    gdats <- gdats[[1]]
  }
  
  return(gdats)
}
#
#
# Mutations type allows "all" or "somatic"
Mutation <- function(ddoc, dlinks, dataset, type="somatic"){
  
  if(!(type %in% c("somatic", "all"))){
    message("Error: Invalid type.")  
    gdat <- NULL
    return(gdat)
  }
  
  keyWord = paste("","Mutation_Packager_Calls",sep="")
  keyWord = paste("//a[contains(@href, '",keyWord,"')]",sep="")
  plinks = XML::xpathSApply(ddoc, keyWord, XML::xmlGetAttr, 'href')
  plinks = plinks[grepl(paste("*.",dataset,"[.]Mutation_Packager_Calls[.]Level_3[.].*.tar[.]gz$",sep=""),plinks)]
  
  if(length(plinks) == 0){
    message("Error: No data available for download. Please ensure the data is available from TCGA. \n")
    dat <- NULL
    return(dat)
  }
  
  # get data processing time-stamp
  timestamp <- unlist(strsplit(dlinks, "/"))
  timestamp <- timestamp[length(timestamp)]
  
  gdats <- list()
  mafs <- list()
  for(i in 1:length(plinks)){
    # Get the right download links
    download_link = paste(dlinks,trim(plinks[i]),sep="/")
    
    if(i==1){
      message("Mutation data will be imported! This may take some time!")
    }
    
    # Download data in tar format and then extract the right text file
    utils::download.file(url=download_link,destfile=paste(dataset,"-Mutation.tar.gz",sep=""),method="auto",quiet = TRUE, mode = "w")
    fileList <- utils::untar(paste(dataset,"-Mutation.tar.gz",sep=""),list=TRUE)
    grepSearch = "*.maf.txt$"
    fileList = fileList[grepl(grepSearch,fileList)] 
    utils::untar(paste(dataset,"-Mutation.tar.gz",sep=""),files=fileList)
    
    ###
    gdat <- c()
    
    for(files in fileList){
      temp <- utils::read.delim(files, header=T, colClasses="character")
      gdat <- rbind(gdat, temp)
    }
    
    delFodler <- paste(getwd(),"/",strsplit(fileList[1],"/")[[1]][1],sep="")
    unlink(delFodler, recursive = TRUE)
    file.remove(paste(dataset,"-Mutation.tar.gz",sep=""))
    
    mafs[[i]] <- gdat
    
    # Make it into nice matrix form
    if(type=="somatic"){
	    gdat <- gdat[gdat[,"Mutation_Status"] == "Somatic" & gdat[,"Variant_Classification"] != "Silent", ]
    }
    
    gdats[[i]] <- Proc.Mutation(gdat)
  }
  
  if(length(gdats) ==1){
    gdats <- gdats[[1]]
    mafs <- mafs[[1]]
  }
  return(gdats=gdats)
}

#
Proc.Mutation <- function(dat){  
  dat <- unique(dat)
  
  # Get patients
  dat[, "Tumor_Sample_Barcode"] <- substr(dat[,"Tumor_Sample_Barcode"], 1, 12)
  dim(dat)
  dat <- unique(dat)
  dim(dat)
  
  genes <- unique(dat[,"Hugo_Symbol"])
  samples <- unique(dat[, "Tumor_Sample_Barcode"])
  
  dat.mat <- t(sapply(samples, function(s) as.numeric(genes %in% (dat[dat[,"Tumor_Sample_Barcode"] == s, ][,"Hugo_Symbol"]))))
  rownames(dat.mat) <- samples
  colnames(dat.mat) <- genes
  
  return(t(dat.mat))
}
#
mRNA_Array <- function(ddoc, dlinks, dataset, type="G450"){
  
  if(!(type %in% c("G450", "U133", "Huex"))){
    message("Error: Invalid miRNA-array platform.")
    gdat <- NULL
    return(gdat)
  }
  
  if(type=="G450"){
    keyWord = paste("","Merge_transcriptome__agilentg4502a_07",sep="")
    keyWord = paste("//a[contains(@href, '",keyWord,"')]",sep="")
    plinks = XML::xpathSApply(ddoc, keyWord, XML::xmlGetAttr, 'href')
    plinks = plinks[grepl(paste("*.",dataset,"[.]Merge_transcriptome__agilentg4502a_.*.__Level_3__unc_lowess_normalization_gene_level__data.Level_3.*.tar[.]gz$",sep=""),plinks)]  
  }
  
  if(type=="U133"){
    keyWord = paste("","Merge_transcriptome__ht_hg_u133a",sep="")
    keyWord = paste("//a[contains(@href, '",keyWord,"')]",sep="")
    plinks = XML::xpathSApply(ddoc, keyWord, XML::xmlGetAttr, 'href')
    plinks = plinks[grepl(paste("*.",dataset,"[.]Merge_transcriptome__ht_hg_u133a__.*.__Level_3__gene_rma__data.Level_3.*.tar[.]gz$",sep=""),plinks)]
  }
  
  if(type=="Huex"){
    keyWord = paste("","Merge_exon__huex_1_0_st_v2",sep="")
    keyWord = paste("//a[contains(@href, '",keyWord,"')]",sep="")
    plinks = XML::xpathSApply(ddoc, keyWord, XML::xmlGetAttr, 'href')
    plinks = plinks[grepl(paste("*.",dataset,"[.]Merge_exon__huex_1_0_st_v2__.*.__Level_3__quantile_normalization_gene__data.Level_3.*.tar[.]gz$",sep=""),plinks)]  
  }
  
  
  if(length(plinks) == 0){
    message("Error: No data available for download. Please ensure this data is available from TCGA. \n")
    dat <- NULL
    return(dat)
  }
  
  # get data processing time-stamp
  timestamp <- unlist(strsplit(dlinks, "/"))
  timestamp <- timestamp[length(timestamp)]
  
  gdats <- list()
  for(i in 1:length(plinks)){
    download_link = paste(dlinks,trim(plinks[i]),sep="/")
    
    message("mRNA Array data will be imported! This may take some time!")
    
    # Download data in tar format and then extract the right text file
    utils::download.file(url=download_link,destfile=paste(dataset,"-mRNAArray.tar.gz",sep=""),method="auto",quiet = TRUE, mode = "w")
    fileList <- utils::untar(paste(dataset,"-mRNAArray.tar.gz",sep=""),list=TRUE)
    grepSearch = paste("*.",dataset,".*__data.data.txt$",sep="")
    fileList = fileList[grepl(grepSearch,fileList)]
    
    utils::untar(paste(dataset,"-mRNAArray.tar.gz",sep=""),files=fileList)
    
    # Rename and clean up
    suff = ifelse(i==1, "", paste("_", i, sep=""))
    fname = paste(dataset,"_", timestamp, "-mRNAArray", suff, ".txt",sep="")
    
    file.rename(from=fileList,to=fname)
    file.remove(paste(dataset,"-mRNAArray.tar.gz",sep=""))
    delFodler <- paste(getwd(),"/",strsplit(fileList,"/")[[1]][1],sep="")
    ##message(delFodler)
    unlink(delFodler, recursive = TRUE)
    
    #Get selected type only
    tmpCols = utils::read.delim(fname,header=F, nrows=1,colClasses="character")
    tmpdat = utils::read.delim(fname, skip=1, sep="\t", stringsAsFactors=F)
    
    gdat <- tmpdat[,-1]
    gdat <- apply(gdat, 2, function(x) suppressWarnings(as.numeric(x)))
    colnames(gdat) <- tmpCols[1,-1]
    colnames(gdat) <- gsub("\\.", "-", colnames(gdat))
    rownames(gdat) <- tmpdat[,1]
    
    message(paste(nrow(gdat) ,"genes have been imported!"))
    gdats[[i]] <- gdat
    names(gdats)[i] <- trim(plinks[i])
    file.remove(fname)
  }
  
  if(length(gdats) == 1){
    gdats <- gdats[[1]]
  }
  else{
    temp <-gdats[[1]]
    for(i in 2:length(gdats)){
      temp <- cbind(temp, gdats[[i]])
    }
    gdats <- temp
  }
  
  return(gdats)
}
#
#
miRNA_Array <- function(ddoc, dlinks, dataset){
  keyWord = paste("","h_mirna_8x15k",sep="")
  keyWord = paste("//a[contains(@href, '",keyWord,"')]",sep="")
  plinks = XML::xpathSApply(ddoc, keyWord, XML::xmlGetAttr, 'href')
  plinks = plinks[grepl("*.Merge_mirna__h_mirna_8x15k.*.data.Level_3.*.tar[.]gz$",plinks)]
  
  if(length(plinks) == 0){
    message("Error: No data available for download. Please ensure the data is available from TCGA. \n")
    dat <- NULL
    return(dat)
  }
  
  # get data processing time-stamp
  timestamp <- unlist(strsplit(dlinks, "/"))
  timestamp <- timestamp[length(timestamp)]
  
  gdats <- list()
  for(i in 1:length(plinks)){
    download_link = paste(dlinks,trim(plinks[i]),sep="/")
    
    message("miRNA Array data will be imported! This may take some time.....")
    
    # Download data in tar format and then extract the right text file
    utils::download.file(url=download_link,destfile=paste(dataset,"-miRNAArray.tar.gz",sep=""),method="auto",quiet = TRUE, mode = "w")
    fileList <- utils::untar(paste(dataset,"-miRNAArray.tar.gz",sep=""),list=TRUE)
    grepSearch = paste("*.",dataset,".*__data.data.txt$",sep="")
    fileList = fileList[grepl(grepSearch,fileList)]
    
    utils::untar(paste(dataset,"-miRNAArray.tar.gz",sep=""),files=fileList)
    
    # Rename and clean up
    suff = ifelse(i==1, "", paste("_", i, sep=""))
    fname = paste(dataset,"_", timestamp, "-miRNAArray", suff, ".txt",sep="")
    
    file.rename(from=fileList,to=fname)
    file.remove(paste(dataset,"-miRNAArray.tar.gz",sep=""))
    delFodler <- paste(getwd(),"/",strsplit(fileList,"/")[[1]][1],sep="")
    #message(delFodler)
    unlink(delFodler, recursive = TRUE)
    
    #Get selected type only
    tmpCols = utils::read.delim(fname,header=F, nrows=1,colClasses="character")
    tmpdat = utils::read.delim(fname, skip=1, sep="\t", stringsAsFactors=F)
    
    gdat <- tmpdat[,-1]
    colnames(gdat) <- tmpCols[1,-1]
    colnames(gdat) <- gsub("\\.", "-", colnames(gdat))
    rownames(gdat) <- tmpdat[,1]
    
    # Here check and remove control probes
    tname <- rownames(gdat)
    controls <- grep("(^dmr)|(^NC)|(Corner)|(Control)|(^hur)|(^mr)", tname) 
    if(length(controls)>0){
      gdat <- gdat[-controls, ]
    }
    # 
    
    message(paste(nrow(gdat) ,"miRNAs have been imported!"))
    gdats[[i]] <- gdat
    names(gdats)[i] <- trim(plinks[i])
    
    file.remove(fname)
  }
  
  if(length(gdats) == 1){
    gdats <- gdats[[1]]
  }
  return(gdats)
}
#
trim <- function (x) gsub("^\\s+|\\s+$", "", x)
#
Clinical <- function(ddoc, dlinks, dataset){
  keyWord = paste("",".Clinical_Pick_Tier1.Level_4",sep="")
  keyWord = paste("//a[contains(@href, '",keyWord,"')]",sep="")
  plinks = XML::xpathSApply(ddoc, keyWord, XML::xmlGetAttr, 'href')
  plinks = plinks[grepl("*.tar[.]gz$",plinks)]
  
  
  if(length(plinks) == 0){
    message("Error: No Clinical data available for download. Please ensure the data is available from TCGA. \n")
    dat <- NULL
    return(dat)
  }
  
  # get data processing time-stamp
  timestamp <- unlist(strsplit(dlinks, "/"))
  timestamp <- timestamp[length(timestamp)]
  
  gdats <- list()
  for(i in 1:length(plinks)){
    download_link = paste(dlinks,trim(plinks[i]),sep="/")
    
    message("Clinical data will be imported. \n")
    
    # Download data in tar format and then extract the right text file
    utils::download.file(url=download_link,destfile=paste(dataset,"-Clinical.tar.gz",sep=""),method="auto",quiet = TRUE, mode = "w")
    fileList <- utils::untar(paste(dataset,"-Clinical.tar.gz",sep=""),list=TRUE)
    fileList = fileList[grepl("*.clin.merged.picked.txt$",fileList)]
    utils::untar(paste(dataset,"-Clinical.tar.gz",sep=""),files=fileList)
    
    # Rename and clean up
    suff = ifelse(i==1, "", paste("_", i, sep=""))
    fname = paste(dataset,"_", timestamp, "-clinical", suff, ".txt",sep="")
    
    file.rename(from=fileList,to=fname)
    file.remove(paste(dataset,"-Clinical.tar.gz",sep=""))
    delFodler <- paste(getwd(),"/",strsplit(fileList,"/")[[1]][1],sep="")
    #message(delFodler)
    unlink(delFodler, recursive = TRUE)
    
    # Read in clinical data
    cli.raw <- utils::read.delim(fname, stringsAsFactors=F, sep="\t")
    cli <- t(cli.raw)
    colnames(cli) <- cli[1,]
    cli <- cli[-1,]
    rownames(cli) <- gsub("\\.", "-", toupper(rownames(cli)))
    colnames(cli) <- gsub("_", "", colnames(cli))
    gdats[[i]] <- cli
    names(gdats)[i] <- trim(plinks[i])
    file.remove(fname)
  }
  
  if(length(gdats) == 1){
    gdats <- gdats[[1]]
  }
  return(gdats)
}
#
#

SurvivalData <- function(dat=dat, cli=cli){
  ctbl <- cli[,c("vitalstatus", "daystodeath", "daystolastfollowup")]
  
  if(class(dat) == "data.frame" | class(dat) == "matrix"){
    tl <- sapply(colnames(dat), function(s) nchar(s))
    if(sum(tl==12) == ncol(dat)){
      tum.pat <- dat
    }
    else{
      tum <- sapply(colnames(dat), function(s) unlist(strsplit(s, "-"))[4]) 
      tum.pat <- dat[, grep("^01", tum)]  
      colnames(tum.pat) <- substr(colnames(tum.pat), 1, 12)
    }
    
    temp1 <- OverallSurvival(ctbl)
    temp2 <- t(tum.pat)
    
    temp1 <- data.frame(rownames(temp1), temp1, stringsAsFactors=F)
    temp2 <- data.frame(rownames(temp2), temp2, stringsAsFactors=F)
    colnames(temp1)[1] <- "bcr"
    colnames(temp2)[1] <- "bcr"
    
    dat.surv <- merge(temp1, temp2, by.x="bcr", by.y="bcr", all=FALSE)
    
    return(dat.surv)
  }
  if(class(dat)=="list"){
    dat.survs <- list()
    for(i in 1:length(dat)){
      temp <- dat[[i]]
      
      tl <- sapply(colnames(temp), function(s) nchar(s))
      if(sum(tl==12) == ncol(temp)){
        tum.pat <- temp
      }
      else{        
        tum <- sapply(colnames(temp), function(s) unlist(strsplit(s, "-"))[4]) 
        tum.pat <- temp[, grep("^01", tum)]  
        colnames(tum.pat) <- substr(colnames(tum.pat), 1, 12)
      }
      
      temp1 <- OverallSurvival(ctbl)
      temp2 <- t(tum.pat)
      
      temp1 <- data.frame(rownames(temp1), temp1, stringsAsFactors=F)
      temp2 <- data.frame(rownames(temp2), temp2, stringsAsFactors=F)
      colnames(temp1)[1] <- "bcr"
      colnames(temp2)[1] <- "bcr"
      
      dat.surv <- merge(temp1, temp2, by.x="bcr", by.y="bcr", all=FALSE)
      
      dat.survs[[i]] <- dat.surv
    }
    return(dat.survs)
  }
}
#
OverallSurvival <- function(ctbl){
  status <- as.numeric(unlist(ctbl[,"vitalstatus"]))
  dtd <- as.numeric(unlist(ctbl[,"daystodeath"]))
  dtlf <- as.numeric(unlist(ctbl[,"daystolastfollowup"]))
  OS <- ifelse(status==1, dtd, dtlf)
  
  overall <- data.frame(status, OS)
  rownames(overall) <- rownames(ctbl)
  return(overall)
}

#
MatrixMerge <- function(dat=dat, cli=cli, cvars="OS"){
  
  # Get clinical
  if(is.null(cli)){
    message("Error: no clinical data available for this disease type.")
    return(NULL)    
  }
  if(!(toupper(cvars) %in% toupper(c("OS", colnames(cli))))){
    message("Error: invalid clinical covariates.")
    return(NULL)  
  }
  
  # Now everything is good
  if(toupper(cvars) == "OS"){
    ctbl <- cli[,c("vitalstatus", "daystodeath", "daystolastfollowup")]
  }
  else{
    colnames(cli) <- toupper(colnames(cli))
    ctbl <- cli[,toupper(cvars), drop=FALSE]
  }
  
  if(class(dat) == "data.frame" | class(dat) == "matrix"){
    if(sum(colnames(dat)[1:3] == c("Gene_Symbol", "Chromosome", "Genomic_Coordinate")) == 3){
      # For methylation
      dat <- dat[, -c(1:3)]
    }
    
    tl <- sapply(colnames(dat), function(s) nchar(s))
    if(sum(tl==12) == ncol(dat)){
      tum.pat <- dat
    }
    else{
      tum <- sapply(colnames(dat), function(s) unlist(strsplit(s, "-"))[4]) 
      tum.pat <- dat[, grep("^01", tum)]  
      colnames(tum.pat) <- substr(colnames(tum.pat), 1, 12)
    }
    
    if(toupper(cvars) == "OS"){
      temp1 <- OverallSurvival(ctbl)
    }
    else{
      temp1 <- ctbl
    }
    temp2 <- t(tum.pat)
    
    temp1 <- data.frame(rownames(temp1), temp1, stringsAsFactors=F)
    temp2 <- data.frame(rownames(temp2), temp2, stringsAsFactors=F)
    colnames(temp1)[1] <- "bcr"
    colnames(temp2)[1] <- "bcr"
    
    dat.cli <- merge(temp1, temp2, by.x="bcr", by.y="bcr", all=FALSE)
    
    return(dat.cli)
  }
  if(class(dat)=="list"){
    dat.clis <- list()
    for(i in 1:length(dat)){
      temp <- dat[[i]]
      
      if(sum(colnames(temp)[1:3] == c("Gene_Symbol", "Chromosome", "Genomic_Coordinate")) == 3){
        # For methylation
        temp <- temp[, -c(1:3)]
      }
      
      tl <- sapply(colnames(temp), function(s) nchar(s))
      if(sum(tl==12) == ncol(temp)){
        tum.pat <- temp
      }
      else{        
        tum <- sapply(colnames(temp), function(s) unlist(strsplit(s, "-"))[4]) 
        tum.pat <- temp[, grep("^01", tum)]  
        colnames(tum.pat) <- substr(colnames(tum.pat), 1, 12)
      }
      
      if(toupper(cvars) == "OS"){
        temp1 <- OverallSurvival(ctbl)
      }
      else{
        temp1 <- ctbl
      }
      temp2 <- t(tum.pat)
      
      temp1 <- data.frame(rownames(temp1), temp1, stringsAsFactors=F)
      temp2 <- data.frame(rownames(temp2), temp2, stringsAsFactors=F)
      colnames(temp1)[1] <- "bcr"
      colnames(temp2)[1] <- "bcr"
      
      dat.cli <- merge(temp1, temp2, by.x="bcr", by.y="bcr", all=FALSE)
      
      dat.clis[[i]] <- dat.cli
    }
    return(dat.clis)
  }
}
#
#
SampleSplit <- function(dat){
  temp <- dat
  
  if(is.null(temp)){
    message("Empy data object")
    return(NULL)
  }
  
  if(class(temp) == "data.frame" | class(temp) == "matrix"){
    if(sum(colnames(temp)[1:3] == c("Gene_Symbol", "Chromosome", "Genomic_Coordinate")) == 3){
      # Process the methylation data
      annot <- temp[,1:3]
      temp <- temp[,-c(1:3)]
      temp.type <- sapply(colnames(temp), function(s) unlist(strsplit(s, "-"))[4])
      primary.tumor <-  cbind(annot, temp[, grep("^01", temp.type)])
      recurrent.tumor <-  cbind(annot, temp[, grep("^02", temp.type)])
      normal <-  cbind(annot, temp[, c(grep("^10", temp.type), grep("^11", temp.type), grep("^12", temp.type))])
    }
    else{
      temp.type <- sapply(colnames(temp), function(s) unlist(strsplit(s, "-"))[4])
      primary.tumor <-  temp[, grep("^01", temp.type)]
      recurrent.tumor <-  temp[, grep("^02", temp.type)]
      normal <-  temp[, c(grep("^10", temp.type), grep("^11", temp.type), grep("^12", temp.type))]  
    }
    return(list(primary.tumor=primary.tumor, recurrent.tumor=recurrent.tumor, normal=normal))
  }
}
#
#
TumorNormalMatch <- function(dat){
  temp <- dat
  
  if(is.null(temp)){
    message("Empy data object")
    return(NULL)
  }
  
  if(class(temp) == "data.frame" | class(temp) == "matrix"){
    temp.type <- sapply(colnames(temp), function(s) unlist(strsplit(s, "-"))[4])
    primary.tumor <-  temp[, grep("^01", temp.type)]
    normal <-  temp[, c(grep("^10", temp.type), grep("^11", temp.type), grep("^12", temp.type))]

    tum.bcr <- substr(colnames(primary.tumor), 1, 12)
    norm.bcr <- substr(colnames(normal), 1, 12)
    
    matching <- intersect(tum.bcr, norm.bcr)
    
    if(length(matching) == 0){
      message("No matching samples on tumor and normal.")
      return(NULL)
    }
    
    if(length(matching) > 0){
      colnames(primary.tumor) <- substr(colnames(primary.tumor), 1, 12)
      colnames(normal) <- substr(colnames(normal), 1, 12)
      
      primary.tumor.match <- primary.tumor[, matching]
      normal.match <- normal[, matching]
      return(list(primary.tumor=primary.tumor.match, normal=normal.match))
    }
  }
}

# Merge multiple object
OMICSBind <- function(dat1, dat2){
  
  if(is.null(dat1) || is.null(dat2)){
    message("Empy data object")
    return(NULL)
  }
  
  if((class(dat1) == "data.frame" || class(dat1) == "matrix") & (class(dat2) == "data.frame" || class(dat2) == "matrix")){
    t1 <- sapply(colnames(dat1), function(s) nchar(s))
    t2 <- sapply(colnames(dat2), function(s) nchar(s))
    
    dat1.tumor <- dat1
    dat2.tumor <- dat2
      
    if(sum(t1>12) == ncol(dat1)){
      dat1.type <- sapply(colnames(dat1), function(s) unlist(strsplit(s, "-"))[4])
      dat1.tumor <- dat1[, grep("^01", dat1.type)]
      colnames(dat1.tumor) <- substr(colnames(dat1.tumor), 1, 12)
    }
    
    if(sum(t2>12) == ncol(dat2)){
      dat2.type <- sapply(colnames(dat2), function(s) unlist(strsplit(s, "-"))[4])
      dat2.tumor <- dat2[, grep("^01", dat2.type)]
      colnames(dat2.tumor) <- substr(colnames(dat2.tumor), 1, 12)
    }
    
    matching <- intersect(colnames(dat1.tumor), colnames(dat2.tumor))
    
    if(length(matching) == 0){
      message("No matched samples")
      return(NULL)
    }
    
    dat1.good <- dat1.tumor[,matching]
    dat2.good <- dat2.tumor[,matching]
    
    rownames(dat1.good) <- paste("d1.", rownames(dat1.good), sep="")
    rownames(dat2.good) <- paste("d2.", rownames(dat2.good), sep="")
    
    mdata <- t(rbind(dat1.good, dat2.good))
    return(list(merged.data=t(mdata), X=dat1.good, Y=dat2.good))
  }
}

