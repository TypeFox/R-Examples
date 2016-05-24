mcmc.qpcr.classic <-
function(fixed=NULL,globalFixed=NULL,random=NULL,globalRandom=NULL,data,controls,genebysample=TRUE,geneSpecRes=FALSE,center=TRUE,...) {
	ngenes=length(levels(data[,"gene"]))

#checking if there are technical reps
	g1=levels(data[,"gene"])[1]
	ss=data[data[,"gene"]==g1,]
	if (length(levels(ss$sample))<length(ss[,1])) reps=1 else reps=0
	
# assembling the fixed effects formula:
	if (!is.null(globalFixed)) {
		gfs=paste(globalFixed,collapse="+")
		gfs=paste(gfs,"+",sep="")
	} else { gfs="" }
	if (is.null(fixed)) {
		ff=paste("count~0+",gfs,"gene",sep="")
	} else {
		ff=gsub('\\s?\\+\\s?',"+gene:",x=fixed,perl=TRUE)
		ff=paste("count~0+gene+",gfs,"+gene:",ff,sep="")
	}	

# assembling the random effects formula:
#	globalRandom="tank"

	rr="~"
	if (!is.null(globalRandom)) {
		grs=paste(globalRandom,collapse="+")
#		grs=paste(grs,"+",sep="")
	} else { grs="" }
	rr=paste(rr,grs,sep="")
	if(genebysample)	random=append(random,"sample")
	for (r in random){
			if (r==random[1] && is.null(globalRandom)) rr=paste(rr,"idh(gene):",r,sep="") else rr=paste(rr,"+idh(gene):",r,sep="") 
	}
# moving controls to the end of the factor list (for variance fixing later)
	controls=controls[length(controls):1]
	relev=c()
	for (g in levels(data[,"gene"])) {
		if (g %in% controls) {
			next
		}
		relev=append(relev,g)
	}
	relev=append(relev,controls)
	data[,"gene"]=factor(data[,"gene"],levels=relev)

	ddn=normalize.qpcr(data,controls,center)
	print(list("FIXED"=ff,"RANDOM"=rr))
	if (geneSpecRes){
		if (rr=="~") {
			mc=MCMCglmm(formula(ff),rcov=~idh(gene):units,data=ddn,...)
		} else { 
			mc=MCMCglmm(formula(ff),random=formula(rr),rcov=~idh(gene):units,data=ddn,...) 
		}
	} else {
		if (rr=="~") {
			mc=MCMCglmm(formula(ff),data=ddn,...)
		} else { 
			mc=MCMCglmm(formula(ff),random=formula(rr),data=ddn,...) 
		}
	}	
	return(mc)
}
