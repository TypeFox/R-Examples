"save.gwaa.data" <-
function(data,phenofile = "pheno.dat", genofile = "geno.raw",human=FALSE) {
	cat("writing phenotypes...\n")
	write.table(data@phdata,file=phenofile,sep=" ",na="NA",row.names=FALSE,col.names=TRUE)
	save.snp.data(data = data@gtdata, genofile =  genofile, human = human)
}

