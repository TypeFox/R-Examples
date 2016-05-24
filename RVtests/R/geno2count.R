geno2count <- function(genotype)
{
if (!is.matrix(genotype)) {
	rownm<- rownames(genotype)
	genotype<- as.matrix(genotype)
	rownames(genotype)<- rownm
	}
indid<- NULL
snpid<- NULL
count<- NULL
for (j in 1:ncol(genotype)){
	indx<- (genotype[,j] > 0)
	if (sum(indx) > 0) {
		gt<- genotype[indx,j, drop = FALSE]
		ord<- order(rownames(gt))
		indid<- c(indid, rownames(gt)[ord])
		count<- c(count, gt[ord])
		snpid<- c(snpid, rep(colnames(genotype)[j], sum(indx)))
		}
	}
countg<- data.frame(indid, snpid, count, stringsAsFactors = FALSE)
rownames(countg)<- 1:nrow(countg)
return(countg)
}

