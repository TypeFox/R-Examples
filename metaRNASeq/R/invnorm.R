invnorm <-
function(indpval, nrep, BHth = 0.05) 
{
    listres = vector("list", 4)
    qnormpval= do.call(cbind,lapply(indpval, FUN=function(x) qnorm(1-x)))
    nrepcorr=t(apply(qnormpval,1,FUN=function(x){ 
    nrep2=nrep
    nrep2[which(is.na(x))]=0
    nrep2}))
    nreptot=apply(nrepcorr,1,sum)
    weight=sqrt(nrepcorr/nreptot)
    wqnormp=weight*qnormpval
    statc=apply(wqnormp,1, FUN=function(x) sum(x,na.rm=TRUE))
	## Added by Andrea for genes filtered in all samples
	## (otherwise returns a value of 0)
	nan.index <- which(apply(wqnormp, 1, function(x) sum(is.nan(x))) == ncol(wqnormp))
	statc[nan.index] <- NA
    rpvalc = 1 - pnorm(statc)
    res = which(p.adjust(rpvalc, method = "BH") <= BHth)
    listres[[1]] = res
    listres[[2]] = statc
    listres[[3]] = rpvalc
    listres[[4]] = p.adjust(rpvalc, method = "BH")
    names(listres) = c("DEindices", "TestStatistic", "rawpval", "adjpval")
    return(listres)
}
