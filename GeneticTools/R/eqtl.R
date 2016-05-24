eQTL <- function(gex, geno, xAnnot=NULL, xSamples=NULL, genoSamples=NULL, windowSize=0.5, method="LM", mc=1, sig=NULL, which=NULL,nper=2000 , verbose=TRUE){

  # Read in the genotype data if name is given, otherwise assume already imported SNPS are given as input
    if(is.character(geno)==TRUE)
    {
      # Check if there is a ped/map pair or vcf given as input
      # CASE: VCF
      if(toupper(substr(geno,nchar(geno)-2, nchar(geno)))=="VCF"){
        stop("The vcf format is not yet supported. Please use the ped/map format.")
        # CASE: PED/MAP
      } else {
        if(verbose==TRUE) cat("Start reading the genotype information at",date(),"(ped/map format assumed)\n")
        genotData <- read.pedfile(file=paste(geno,".ped",sep=""),snps=paste(geno,".map",sep=""))
      }
    } else {
      genotData <- geno
    }

  # Input checks
    if(is.vector(gex) & is.null(xAnnot))
    {
      warning("No annotations given, we will test all given SNPs against all given expressions!\n")
      xAnnot <- data.frame(gene="",Chr=as.character(genotData$map[1,1]),Start=genotData$map[1,4],End=genotData$map[1,4])
      windowSize=NULL
    }

  # If no separate genOSamples object is given, we take those from the snpStats object:
    if(is.null(genoSamples)) genoSamples <-  as.character(genotData$fam$member)
  
  # If the annotations are given as data frame, we will transform them into a list
    if(is.data.frame(xAnnot)){
      if(is.factor(xAnnot[,1])) xAnnot[,1] <- as.character(xAnnot[,1])
      if(is.factor(xAnnot[,2])) xAnnot[,2] <- as.character(xAnnot[,2])
      if(verbose==TRUE) cat("We will transform the gene annotations into a list at",date(),"!\n")
      xAnnot <- makeAnnotList(xAnnot)
      if(verbose==TRUE) cat("We finished to transform the gene annotations at",date(),"\n")
    }

  # Take only those genes from the annotation list that were requested
    if(!is.null(which)) xAnnot <- xAnnot[is.element(names(xAnnot),which)]
  
  # Input checks
    single <- FALSE
    th <- windowSize
    method <- match.arg(method,c("LM","directional"))

  # If only one gene is given it could be a vector, transform it here to a column matrix
    if(is.vector(gex)){
      single <- TRUE
      tempNames <- names(gex)
      if(!is.null(xSamples)) tempNames <- xSamples
      gex <- matrix(gex,ncol=1)
      rownames(gex) <- tempNames
      colnames(gex) <- names(xAnnot)
      gexColNames <- colnames(gex)
    }
 
  # In case that the row names have been changed, bring them into an order
    rownames(genotData$map) <- 1:nrow(genotData$map)

  # Sample statistics
    overlap <- is.element(rownames(gex),genoSamples)
    olPerc <- round(sum(overlap)/nrow(gex)*100,1)
    if(sum(overlap)==0) stop("No matching expression / genotype sample names!\n")
    if(verbose==TRUE) cat("We have for",olPerc,"% of the samples in the expression data the genotype information. \n")

  # Location statistics
    overlap <- is.element(colnames(gex),names(xAnnot))
    olPerc <- round(sum(overlap)/ncol(gex)*100,1)
    if(sum(overlap)==0) stop("No matching expression probe names / probe name annotations!\n")
    if(verbose==TRUE) cat("We have for",olPerc,"% of the expression data the annotations. \n")

  # Probe statistics
    matchingProbes <- colnames(gex)[is.element(colnames(gex),names(xAnnot))]
    if(verbose==TRUE) cat("We will investigate for",length(matchingProbes),"genes possible eQTLs! \n")
    result <- list()

  # Reducing the expression data to those rows, where we have also genotype information available
    gex <- gex[is.element(rownames(gex),genoSamples),]
    if(single==TRUE){
      gex <- t(t(gex))
      colnames(gex) <- gexColNames
    }

  # Now go through all possible probes
  eqtl <- list()
  for(probeRun in 1:length(matchingProbes))
  {
    # Do that for each possible location of the probe (might not be unique...)
    tempAnnot <- xAnnot[[which((names(xAnnot)==matchingProbes[probeRun])==TRUE)]]
    eqtlTemp <- list()
    for(tempRun in 1:nrow(tempAnnot))
    {
      # Temporary values
      SNPloc <- getSNPlocations(genotInfo=genotData$map,annot=tempAnnot[tempRun,],th=th)
      SNPmatrix <- genotData$genotypes[,SNPloc$SNPcol]
      genoGroups <- as(SNPmatrix,"numeric")
      genoGroups <- rearrange(genoGroups,rownames(gex),genoSamples)
      
      # Run eQTL only if there are SNPs inside the window
      if(dim(SNPloc$SNPloc)[1]>0){
      # Type of eQTL
      if(method=="LM"){
	if(is.null(sig))
        {
	  eqtlTemp[[tempRun]] <- list(ProbeLoc=rep(tempRun,ncol(genoGroups)),TestedSNP=SNPloc[[1]],p.values=eqtlLM(genoGroups,gex[,probeRun]))
	} else {
	  p.values <- eqtlLM(genoGroups,gex[,probeRun])
	  pPos <- p.values<sig
	  eqtlTemp[[tempRun]] <- cbind(SNPloc[[1]][pPos,c(1,2,4)],p.values[pPos])
	}
      } else if(method=="directional"){
        if(is.null(sig))
        { 
	  eqtlTemp[[tempRun]] <- list(ProbeLoc=rep(tempRun,ncol(genoGroups)),TestedSNP=SNPloc[[1]],p.values=eqtlDir.P(genoGroups,gex[,probeRun],mc=mc,nper=nper))
	} else {
	  p.values <- eqtlDir.P(genoGroups,gex[,probeRun],mc=mc,nper=nper)
	  pPos <- p.values<sig
	  eqtlTemp[[tempRun]] <- cbind(SNPloc[[1]][pPos,c(1,2,4)],p.values[pPos])
	}
      } 
     } else {
        warning("There were no variants within the window")
        eqtlTemp[[tempRun]] <- list(ProbeLoc=-1,TestedSNP=-1,p.values=-1)
     }
    }
    
    # Join the output
    if(is.null(sig))
    {
      eqtl[[probeRun]] <- joinEQTL(eqtlTemp)
      eqtl[[probeRun]]$GeneInfo <- tempAnnot
    } else {
 #     if(is.null(windowSize))
 #     {
	bedTemp <- joinEQTLsig(eqtlTemp)
	bedTemp <- cbind(bedTemp,rep(matchingProbes[probeRun],nrow(bedTemp)))
	colnames(bedTemp) <- c("chr", "SNP", "Location", "p.value", "Assoc.Gene")
	eqtl[[probeRun]] <- bedTemp
 #     } else {
#	eqtl[[probeRun]] <- joinEQTLsig(eqtlTemp)
 #     }
    }
    if(verbose==TRUE) cat ("We investigated",probeRun,"gene eQTLs at",date(),"\n")
  }

  # Return the result
  if(is.null(sig))
  {
    names(eqtl) <- matchingProbes
    result <- list(bed=NULL,eqtl=eqtl,gex=gex, geno=geno, xAnnot=xAnnot, xSamples=xSamples, genoSamples=genoSamples, windowSize=windowSize, method=method, mc=mc, which=which, type="full")
  } else {
    result <- list(bed=joinEQTLsig(eqtl),gex=gex, geno=geno, xAnnot=xAnnot, xSamples=xSamples, genoSamples=genoSamples, windowSize=windowSize, method=method, mc=mc, which=which, type="sig")
  }
  class(result) <- "eqtl"
  result
}



