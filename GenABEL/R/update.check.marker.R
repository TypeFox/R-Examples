"update.check.marker"  <-
function(preo,curo) {
	out <- list()
	out$Xerrtab <- preo$Xerrtab
# snp chars
	out$nofreq <- unique(c(preo$nofreq,curo$nofreq))
	out$nocall <- unique(c(preo$nocall,curo$nocall))
	out$nohwe <- unique(c(preo$nohwe,curo$nohwe))
	out$Xmrkfail <- unique(c(preo$Xmrkfail,curo$Xmrkfail))
	out$Xmrkfail <- unique(c(preo$Xmrkfail,curo$Xmrkfail))
# rely that redundancy is checked only in the last run
	out$details.redundancy <- curo$details.redundancy
	out$redundant <- curo$redundant
# get all failes
	snpfail <- unique(c(out$nofreq,out$nocall,out$nohwe,out$Xmrkfail,out$redundant))
# personal
	out$hetfail <- unique(c(preo$hetfail,curo$hetfail))
	out$idnocall <- unique(c(preo$idnocall,curo$idnocall))
	out$ibsfail <- unique(c(preo$ibsfail,curo$ibsfail))
#	out$Xidfail <- unique(c(preo$Xidfail,curo$Xidfail))
	out$isfemale <- unique(c(preo$isfemale,curo$isfemale))
	out$ismale <- unique(c(preo$ismale,curo$ismale))
	out$otherSexErr <- unique(c(preo$otherSexErr,curo$otherSexErr))
	out$isXXY <- unique(c(preo$isXXY,curo$isXXY))
#	idfail <- unique(c(out$hetfail,out$idnocall,out$ibsfail,out$Xidfail,out$isXXY,out$isfemale,out$ismale))
	idfail <- unique(c(out$hetfail,out$idnocall,out$ibsfail,out$isXXY,out$isfemale,out$ismale,out$otherSexErr))
# generate snpok and idok
	if (length(snpfail)>0) 
		out$snpok <- preo$snpok[!(preo$snpok %in% snpfail)]
	else 
		out$snpok <- preo$snpok
	if (length(idfail)>0) 
		out$idok <- preo$idok[!(preo$idok %in% idfail)]
	else 
		out$idok <- preo$idok
# return object
	out
}
