rowMAFs <- function(x, check=TRUE){
	if(!is.matrix(x))
		stop("x must be a matrix.")
	if(check && any(!x %in% c(1:3, NA), na.rm=TRUE))
		stop("The genotypes must be coded by 1 (for the homozygous reference),\n",
			"2 (heterozygous), and 3 (homozygous variant).")
	maf <- rowMeans(x==3, na.rm=TRUE) + 0.5 * rowMeans(x==2, na.rm=TRUE)
	if(any(maf>0.5))
		warning("At least one MAF is larger than 0.5.", call.=FALSE)	
	maf
}

