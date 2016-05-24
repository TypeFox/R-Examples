mcmc.otu <-
function(fixed=NULL,random=NULL,data,y.scale="proportion",globalMainEffects="remove",vprior="uninf",...) {
	if (globalMainEffects=="keep") { 
		data=data[data$otu!="summ",] 
		data$otu=factor(data$otu,levels=unique(data$otu))
	}
	ngenes=length(levels(data[,"otu"]))
	if(vprior=="uninf") { 
		vstr.g1=list(V=1, nu=0)
		vstr.gs=list(V=diag(ngenes), nu=0)
	} else {
		if (vprior=="iw") {
			vstr.g1=list(V=1, nu=0.002, alpha.mu=0, alpha.V=1000)
			vstr.gs=list(V=diag(ngenes), nu=ngenes-1+0.02)
		} else {
			if (vprior=="iw01") {
				vstr.g1=list(V=1, nu=0.002, alpha.mu=0, alpha.V=1000)
				vstr.gs=list(V=diag(ngenes)*0.1, nu=ngenes-1+0.02)
			} else { 
				stop("vprior is not recognized, should be uninf, iw, or iw01")
				}
		}
	}

# assembling the fixed effects formula:
	if (is.null(fixed)) {
		if (globalMainEffects=="keep") {
			ff="count~0+otu"
		} else {
			if (y.scale=="proportion") { 
				ff="count~otu"			
			} else {
				ff="count~0+otu"
			}
		}	
	} else {
		ff=gsub('\\s?\\+\\s?',"+otu:",x=fixed,perl=TRUE)
		if (globalMainEffects=="keep") {
			ff=paste("count~0+otu+otu:",ff,sep="")
		} else {
			if (y.scale=="proportion") { 
				ff=paste("count~otu+",fixed,"+otu:",ff,sep="")	
			} else {
				ff=paste("count~0+otu+",fixed,"+otu:",ff,sep="")					
			}		
		}
	}	

# assembling the random effects formula:
	rr="~sample"
	for (r in random){
			rr=paste(rr,"+idh(otu):",r,sep="")
	}

	Gstruct=list(G1=vstr.g1)
	va.r=diag(ngenes)
	va.r.cov=diag(ngenes)*0
	Rr=list(V=diag(ngenes), nu=ngenes-1+0.02) 

	if(length(random)>0){
		for (ri in 1:length(random)) {
			Gstruct[[paste("G",1+ri,sep="")]]=vstr.gs
		}
	}
	prior=list(R=Rr,G=Gstruct)
	print(list("PRIOR"=prior,"FIXED"=ff,"RANDOM"=rr))
	mc=MCMCglmm(formula(ff),random=formula(rr),rcov=~idh(otu):units,data=data,prior=prior,family="poisson",...)
	return(mc)
}
