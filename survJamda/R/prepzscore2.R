prepzscore2 <-
function (common.gene, geno.files,surv.data,x,y)
{
	fLst = vector("list", length(x))
	m = NULL
	surv = censor = NULL

	for (i in x){
		fLst[[i]] = get(geno.files[i])[,common.gene]
		m = rbind(m, fLst[[i]])
		res = comb.surv.censor(geno.files,c(i), surv.data)
		surv = c(surv, res$surv)
		censor = c(censor, res$censor)
	
	}
	phyno = list(surv = surv, censor = censor)

	lst = excl.missing(m, phyno)
	m = znorm(lst$m)

        return(list(mat = m, phyno = lst$phyno))
}

