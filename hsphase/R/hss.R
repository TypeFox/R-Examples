readGenotype <- function(genotypePath, separatorGenotype = " ", check = TRUE)
{
	debut <- proc.time()
	print("Reading the genotype file ...")
	nline <- readLines(genotypePath, 2)
	if (check) 
		if (length(unlist(strsplit(nline[1], split = separatorGenotype))) != (length(unlist(strsplit(nline[2], split = separatorGenotype)))) - 1)
		{
			print(paste("Line 1:", length(unlist(strsplit(nline[1], split = separatorGenotype)))))
			print(paste("Line 2:", length(unlist(strsplit(nline[2], split = separatorGenotype)))))
			print("The first line must have one element less than line two, possible reason: ")
			stop("The animal IDs must not have a column name (The first line must only contain SNPs' names)")
		}
	ncols <- length(unlist(strsplit(nline[2], split = separatorGenotype)))
	if (check) 
		if (ncols < 101)
		{
			stop("The genotype file must have more than 100 columns.")
		}
	dataTemp <- scan(genotypePath, what = c("character", rep("numeric", ncols - 1)), skip = 1, , sep = separatorGenotype)
	index <- seq(from = 1, to = length(dataTemp), by = ncols)
	indexneg <- setdiff(1:length(dataTemp), index)
	if (check) 
		if (length(indexneg)%%(ncols - 1) != 0) 
			stop("All rows must have the same number of elements (The separator may set incorrectly)")
	genotype <- matrix(as.integer(dataTemp[indexneg]), ncol = (ncols - 1), byrow = T)
	rownames(genotype) <- dataTemp[index]
	colnames(genotype) <- unlist(strsplit(nline[1], split = separatorGenotype))
	if (check) 
		if (any(is.na(genotype))) 
			stop("Genotype must contain only integer data (After excluding SNP namse and IDs")
	if (check) 
		if (length(genotype[genotype != 9 & genotype != 0 & genotype != 2 & genotype != 1]) > 0) 
			stop("Genotype must contain only 0, 1, 2 or 9")
	if (check) 
		if (nrow(genotype) < 4)
		{
			stop("The genotype file must have more than 4 rows")
		}
	print(proc.time() - debut)
	print("Finished!")
	genotype
}

hss <- function(pedigree, genotype, check = TRUE)
{
	debut <- proc.time()
	
	if (check) 
	{
		if(length(rownames(genotype)[duplicated(rownames(genotype))])>0)
		{
			stop("There are duplicated ID in the input file")
		}
		if (any(is.na(genotype))) 
			stop("Genotype must contain only integer data (After excluding SNP namse and IDs")
		if (length(genotype[genotype != 9 & genotype != 0 & genotype != 2 & genotype != 1]) > 0) 
			stop("Genotype must contain only 0, 1, 2 or 9")
		if (nrow(genotype) < 4)
		{
			stop("The genotype file must have more than 4 rows")
		}
		if (ncol(pedigree) < 2) 
			stop("Pedigree must have at least 2 columns")
	}
	pedigree <- as.matrix(pedigree[, 1:2])
	pedigree <- pedigree[!(pedigree[, 2] == ""), ]
	pedigree <- pedigree[pedigree[, 2] != 0, ]
	pedigree <- pedigree[pedigree[, 1] %in% rownames(genotype), ]
	## sires <- aggregate(pedigree[, 1], by = list(pedigree[, 2]), function(x) x, simplify = F)
	sires <- tapply(pedigree[, 1], pedigree[, 2], function(x) x, simplify = F) # Faster
	nhalfsibs <- as.integer(lapply(sires, length))
	sires <- sires[which(nhalfsibs > 3)]
	gc()
	print(paste("There are", nrow(sires), "half-sib groups (Based on the pedigree and genotype)"))
	halfsib <- list(nrow(sires))
	genotype <- as.matrix(genotype)
	if (check) 
		if (length(genotype[genotype != 9 & genotype != 0 & genotype != 2 & genotype != 1]) > 0) 
			stop("GenotypeMatrix must contain only 0, 1, 2 or 9")
	ncols <- numeric(nrow(genotype))
	print("Half-sib groups separator")
	for (i in 1:nrow(sires))
	{
		print(paste(i, nrow(sires)))
		halfsib[[i]] <- as.matrix(genotype[which(rownames(genotype) %in% as.character(as.vector(sires[[i]]))), ])
		names(halfsib)[i] <- names(sires)[i]
		if (ncol(halfsib[[i]]) == 1) 
			halfsib[[i]] <- t(halfsib[[i]])
		halfsib[[i]] <- na.omit(halfsib[[i]])
		ncols[i] <- ncol(halfsib[[i]])
	}
	halfsib[lapply(halfsib, nrow) < 4 | lapply(halfsib, ncol) < 100] <- NULL
	print(paste("There are", length(halfsib), "half-sib groups"))
	print(proc.time() - debut)
	print("Finished!")
	gc()
	halfsib
} 
