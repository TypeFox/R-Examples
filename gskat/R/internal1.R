#######################################################################
#These R functions were taken from the package by Gao et al.
#            http://simplem.sourceforge.net/    
#            We have put it here instead of direting to the website or installing the 
#            package to make it easier to use.
########################################################################





#============================================================================
# simpleM_Ex.R 

#============================================================================
# License:  GPL version 2 or newer. 
# NO warranty. 

#============================================================================
# citation: 
#
# Gao X, Starmer J and Martin ER (2008) A Multiple Testing Correction Method for
# Genetic Association Studies Using Correlated Single Nucleotide Polymorphisms. 
# Genetic Epidemiology 32:361-369
#
# Gao X, Becker LC, Becker DM, Starmer J, Province MA (2009) Avoiding the high 
# Bonferroni penalty in genome-wide association studies. Genetic Epidemiology 
# (Epub ahead of print) 

#============================================================================
# readme: 
# example SNP file format:
# row => SNPs
# column => Unrelated individuals 

# The data file should contain only POLYMORPHIC SNPs. 

# Missing values should be imputed. 
# There should be NO missing values in the SNP data file.
# SNPs are coded as 0, 1 and 2 for the number of reference alleles. 
# SNPs are separated by one-character spaces. 

# You may need to change file path (search for "fn_In" variable) 
# depending on where your snp file is stored at.

#============================================================================
# Meff through the PCA approach 
# use a part of the eigen values according to how much percent they contribute
# to the total variation 
Meff_PCA <- function(eigenValues, percentCut){
	totalEigenValues <- sum(eigenValues)
	myCut <- percentCut*totalEigenValues
	num_Eigens <- length(eigenValues)
	myEigenSum <- 0
	index_Eigen <- 0
	
	for(i in 1:num_Eigens){
		if(myEigenSum <= myCut){
			myEigenSum <- myEigenSum + eigenValues[i]
			index_Eigen <- i
		}
		else{
			break
		}
	}	
	return(index_Eigen)
}

#============================================================================
# infer the cutoff => Meff
inferCutoff <- function(dt_My){
	CLD <- cor(dt_My)
	eigen_My <- eigen(CLD)
		
	# PCA approach
	eigenValues_dt <- abs(eigen_My$values)
	Meff_PCA_gao <- Meff_PCA(eigenValues_dt, PCA_cutoff)
	return(Meff_PCA_gao)
}

#============================================================================
PCA_cutoff <- 0.995

#============================================================================
# fix length, simpleM

simpleM<-function (Z) {

	#fn_In <- "D:/simpleM_Ex/snpSample.txt"				# <---- change path here!!!
	#mySNP_nonmissing <- read.table(fn_In, colClasses="integer")
	mySNP_nonmissing<- as.data.frame(t(Z))

	numLoci <- length(mySNP_nonmissing[, 1])

	simpleMeff <- NULL
	fixLength <- 133 
	i <- 1
	myStart <- 1
	myStop <- 1
	while(myStop < numLoci){
		myDiff <- numLoci - myStop 
		if(myDiff <= fixLength) break
		
		myStop <- myStart + i*fixLength - 1
		snpInBlk <- t(mySNP_nonmissing[myStart:myStop, ])
		MeffBlk <- inferCutoff(snpInBlk)
		simpleMeff <- c(simpleMeff, MeffBlk)
		myStart <- myStop+1
	}
	snpInBlk <- t(mySNP_nonmissing[myStart:numLoci, ])
	MeffBlk <- inferCutoff(snpInBlk)
	simpleMeff <- c(simpleMeff, MeffBlk)

	#cat("Total number of SNPs is: ", numLoci, "\n")
	#cat("Inferred Meff is: ", sum(simpleMeff), "\n")

  return(simpleMeff)
}

#============================================================================
# end 
