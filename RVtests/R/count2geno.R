count2geno <-
function(cgeno, indid)
{
#cgeno names=c("indid", "snpid", "count")
if (missing(indid)) indid<- unique(cgeno[,1])
else stopifnot(all(unique(cgeno[,1]) %in% indid))

snpid<- unique(cgeno[,2])
genotype<- matrix(0, nrow = length(indid), ncol = length(snpid))
dimnames(genotype)<- list(indid, snpid)
for (i in 1:nrow(cgeno))
	genotype[as.character(cgeno[i,1]), as.character(cgeno[i,2])]<- cgeno[i,3]

return(genotype)
}

