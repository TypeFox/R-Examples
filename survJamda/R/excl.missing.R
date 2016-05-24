excl.missing <-
function(mat, phyno)
{
	mat = mat[!is.na(phyno$surv),]
	phyno$censor = phyno$censor[!is.na(phyno$surv)]
	phyno$surv = phyno$surv[!is.na(phyno$surv)]

	return(list(mat = mat, phyno = phyno))
}

