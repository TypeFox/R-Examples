prepcombat <-
function (common.gene,geno.files,surv.data, batchID,x,y)
{
	m = NULL
	for (i in union(x,y)){
		m = rbind(m, get(geno.files[i])[,common.gene])
		if (i == x[length(x)]) train.ind = 1:nrow(m)
	}
	phyno = comb.surv.censor(geno.files,union(x,y),surv.data)

	lst = excl.missing(m,phyno)
	writeSamples(lst$mat, batchID, "sampleMerge")

       writeGeno(lst$mat, "genoMerge")
       mat = compute.combat ("genoMerge", "sampleMerge")
        return(list(mat = t(mat), phyno = lst$phyno, train.ind = train.ind))
}

