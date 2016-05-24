# This script performs a QTL study for QTL studies 

QTL <- function(pheno, geno, phenoSamples=NULL, genoSamples=NULL, method="LM", mc=1, sig=NULL, nper=2000 , verbose=TRUE){

  method <- match.arg(method,c("LM","directional"))

  # Read in the genotype data if name is given, otherwise assume already imported SNPS are given as input
    if(is.character(geno)==TRUE)
    {
      if(verbose==TRUE) cat("Start reading the genotype information at",date(),"\n")
      genotData <- read.pedfile(file=paste(geno,".ped",sep=""),snps=paste(geno,".map",sep=""))
    } else {
      genotData <- geno
    }

  # The input for pheno has to be numeric, as I assume it later, when I loop over the phenotypes 
    if(!is.numeric(pheno)) stop("There are non-numeric entries in the input 'pheno'. Please fix that before continuing by dropping those columns or apply a as.numeric() to them.")

  # If no separate genOSamples object is given, we take those from the snpStats object:
    if(is.null(genoSamples)) genoSamples <-  as.character(genotData$fam$member)

  # If only one phenotype is given, it could be a vector, transform it here to a column matrix
    single <- FALSE
    if(is.vector(pheno)){
      tempNames <- names(pheno)
      single <- TRUE
      if(!is.null(phenoSamples)) tempNames <- phenoSamples
      pheno <- matrix(pheno,ncol=1)
      rownames(pheno) <- tempNames
      colnames(pheno) <- "Phenotype"
      phenoColNames <- colnames(pheno)
    }

  # Sample statistics
    overlap <- is.element(rownames(pheno),genoSamples)
    olPerc <- round(sum(overlap)/nrow(pheno)*100,1)
    if(sum(overlap)==0) stop("No matching phenotype / genotype sample names!\n")
    if(verbose==TRUE) cat("We have for",olPerc,"% of the samples in the phenotype data also the genotype information. \n")

  # Reducing the expression data to those rows, where we have also genotype information available
    pheno <- pheno[is.element(rownames(pheno),genoSamples),]
    if(single==TRUE){
      pheno <- t(t(pheno))
      colnames(pheno) <- phenoColNames
    }

    
}