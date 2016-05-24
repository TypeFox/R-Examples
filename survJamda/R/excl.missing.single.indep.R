excl.missing.single.indep <-
function(geno.files,ind,surv.data, common.gene)
{

	mat = get(geno.files[ind])[,common.gene]
	phyno = comb.surv.censor(geno.files,c(ind),surv.data)

	mat = mat[!is.na(phyno$surv),]

	phyno$censor = phyno$censor[!is.na(phyno$surv)]
	phyno$surv = phyno$surv[!is.na(phyno$surv)]

	return(list(mat = mat, phyno = phyno))
}

