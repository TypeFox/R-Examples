
nullSim<- function(y, x, gdat, prdat, ped, gmap, hap,
	method=c("permutation","gene dropping"), vc=NULL, intcovar=NULL,
	test = c("None","F","Chisq"), minorGenoFreq=0.05, rmv=TRUE,
	recode.pedigree=FALSE, gr=2, ntimes=10){
	matr<- NULL
	method<- match.arg(method)
	test<- match.arg(test)
	if(method=="gene dropping"){
		pedR<- ped
		if(recode.pedigree) pedR<- pedRecode(pedR)
		if(missing(prdat)){
			ids<- rownames(gdat)
			if(any(!is.element(ids, pedR$old)))
				stop("Not all sample IDs in both 'prdat' and 'ped'?")
			pos<- gmap
			for(n in 1:ntimes){
				gdatTmp<- genoSim(pedR, gmap=gmap, ids=ids, hap=hap,
					method="Haldane", recode.pedigree=FALSE)
				llkTmp<- scanOne(y=y, x=x, gdat=gdatTmp, vc=vc, intcovar=intcovar,
					test=test, minorGenoFreq=minorGenoFreq, rmv=rmv)

				if(minorGenoFreq <= 0 && rmv) matr<- rbind(matr, llkTmp$p)
				else if(test != "None") stop("test should be 'None'")
				else matr<- rbind(matr, max(llkTmp$p))
			}
		}else{
			ids<- dimnames(prdat$pr)[[1]]
			if(any(!is.element(ids, pedR$old)))
				stop("Not all sample IDs in both 'prdat' and 'ped'?")
			pos<- data.frame(snp=prdat$snp, chr=prdat$chr, dist=prdat$dist)
			for(n in 1:ntimes){
				gdatTmp<- genoSim(pedR, gmap=gmap, ids=ids, hap=hap,
					method="Haldane", recode.pedigree=FALSE)
				prd<- genoProb(gdatTmp, gmap, gr=gr, pos=pos, method="Haldane", verbose = FALSE)
				llkTmp<- scanOne(y=y, x=x, prdat=prd, vc=vc, intcovar=intcovar,
					test=test, minorGenoFreq=minorGenoFreq, rmv=rmv)

				if(minorGenoFreq <= 0 && rmv) matr<- rbind(matr, llkTmp$p)
				else if(test != "None") stop("test should be 'None'")
				else matr<- rbind(matr, max(llkTmp$p))
			}
		}
	}else if(method=="permutation"){
		if(missing(prdat)){
			if(missing(gdat)) stop("either 'gdat' or 'prdat' should be provided")
			for(n in 1:ntimes){
				idx<- sample(1:nrow(gdat),replace=FALSE)
				llkTmp<- scanOne(y=y, x=x, gdat=gdat[idx,], vc=vc, intcovar=intcovar,
					test=test, minorGenoFreq=minorGenoFreq, rmv=rmv)

				matr<- rbind(matr, llkTmp$p)
			}
		}else{
			prd<- prdat
			for(n in 1:ntimes){
				idx<- sample(1:dim(prdat$pr)[1], replace=FALSE)
				prd$pr<- prdat$pr[idx,,]
				llkTmp<- scanOne(y=y, x=x, prdat=prd, vc=vc, intcovar=intcovar,
					test=test, minorGenoFreq=minorGenoFreq, rmv=rmv)

				matr<- rbind(matr, llkTmp$p)
			}
		}
	}else stop("permutation or gene dropping.")

	matr
}


