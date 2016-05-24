pool.zscores <-
function(common.gene, s, geno.files, surv.data)
{
        zstat = list(length = length(s))
        
        for (i in s){
  		lst.pheno = comb.surv.censor(geno.files,c(i), surv.data)
                surv = lst.pheno$surv
		censor = lst.pheno$censor

                curr_file= get(geno.files[i])[!is.na(surv),]
                censor= censor[!is.na(surv)]
                surv= surv[!is.na(surv)]
                zstat[[i]] = featureselection.meta(curr_file[,common.gene],surv,
censor) }
        return (zstat)
}

