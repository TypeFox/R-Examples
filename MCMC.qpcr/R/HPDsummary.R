HPDsummary <-
function(model,data,xgroup=NULL,genes=NA,relative=FALSE,
log.base=2,summ.plot=TRUE,ptype="z",...) {
	
#model=msr3;data=qs;xgroup=NULL;genes=NA;relative=FALSE;log.base=2;summ.plot=TRUE;ptype="z"

	mm=model;base=log.base
	gene.results=list()
	trms=attr(terms(mm$Fixed$formula),"term.labels")[grep("gene:",attr(terms(mm$Fixed$formula),"term.labels"))]
	trms=sub("gene:","",trms)
	sols=colnames(mm$Sol)
	if (is.na(genes[1])) { 
		genes=sub("gene","",sols[grep("gene.*$",sols)])
		genes=unique(sub(":.*","",genes))	
	}
	facts =list()
	for (t in trms) {
		if (!grepl(":",t)) { facts =append(facts,list(levels(data[,t])))}
	}
	names(facts)=trms[1:length(facts)]
	nfactors=length(facts)
	if (nfactors>2) { 
		stop("not implemented for more than 2 crossed factors")
	}
	gsols=c();fac1=c();fac2=c();samps=c();skips=c()
	
	globalInt=rep(0,length(mm$Sol[,1]))	
	if(sols[1]=="(Intercept)") { 
		globalInt=rep(mean(mm$Sol[,"(Intercept)"]),length(mm$Sol[,1]))
	} 
	
	interaction=0
	if (nfactors==2) {
		sol=paste("gene",genes[2],":",names(facts)[1],facts[[1]][2],":",names(facts)[2],facts[[2]][2],sep="")
		if (sol %in% sols) {
			interaction=1
		}
	}
	
	for (gene in genes) {
		sol=paste("gene",gene,sep="")
		if (!(sol %in% sols)) colnames(mm$Sol)[1]=sol
	}
	sols=colnames(mm$Sol)	

	for (gene in genes) {
		fac1g=c();fac2g=c();sampsg=c();skip=FALSE
		for (lev1 in 1:length(facts[[1]])) {
			if (nfactors==2) {
				for (lev2 in 1:length(facts[[2]])) {
					if (lev1==1 & lev2==1) { 
						sol=paste("gene",gene,sep="")
						if(sum(grep(sol,sols))==0) {
							skip=TRUE
							next
						}
						samp=(globalInt+mm$Sol[,sol])*as.numeric(!relative)
						int0=samp
					} else {
						if (lev2==1) { 
							sol=paste("gene",gene,":",names(facts)[1],facts[[1]][lev1],sep="")
							if(sum(grep(sol,sols))==0) {
								skip=TRUE
								next
							}
							
							samp=int0+mm$Sol[,sol]
							int2=mm$Sol[,sol]
						} else {
							if (lev1==1) { 
								sol=paste("gene",gene,":",names(facts)[2],facts[[2]][lev2],sep="") 
								if(sum(grep(sol,sols))==0) {
									skip=TRUE
									next
								}
								samp=int0+mm$Sol[,sol]
								int1=mm$Sol[,sol]
							} else {
								sol=paste("gene",gene,":",
								names(facts)[1],facts[[1]][lev1],
								":",names(facts)[2],facts[[2]][lev2],sep="") 
								if(sum(grep(sol,sols))==0 & interaction==1) {
									skip=TRUE
									next
								}
								if (interaction==1) {
									samp=int0+int1+int2+mm$Sol[,sol]
								} else { 
									samp=int0+int1+int2
								}
							}
						}
					}
		#			print(paste(lev1,lev2,sol))
					if (skip) { 
#						genes=genes[!(genes %in% gene)]						
						next 
					}
					gsols=append(gsols,gene)
					fac1g=append(fac1g,facts[[1]][lev1])
					fac2g=append(fac2g,facts[[2]][lev2])
					sampsg=cbind(sampsg,samp)
				}
			} else {
				if (lev1==1) { 
						sol=paste("gene",gene,sep="") 
						if(sum(grep(sol,sols))==0) {
							skip=TRUE
							next
						}
						samp=(globalInt+mm$Sol[,sol])*as.numeric(!relative)
						int0=samp
				} else {
					sol=paste("gene",gene,":",
					names(facts)[1],facts[[1]][lev1],sep="")
					if(sum(grep(sol,sols))==0) {
						skip=TRUE
						next
					}
					samp=int0+mm$Sol[,sol]
				}
				if (skip) { 
#					genes=genes[!(genes %in% gene)]
					next 
				}
				gsols=append(gsols,gene)
				fac1g=append(fac1g,facts[[1]][lev1])
				sampsg=cbind(sampsg,samp)
			}
		}
		if (skip) { 
			skips=append(skips,gene)
			next 
		}
		fac1=append(fac1,fac1g)
		samps=cbind(samps,sampsg)
		if (nfactors==2) { fac2=append(fac2,fac2g) }
		sampsg=data.frame(sampsg)
		if (nfactors==2) {
			gres=matrix(nrow=length(fac1g),ncol=length(fac1g),
			dimnames=list("pvalue"=paste(names(facts)[1],fac1g,":",
			names(facts)[2],fac2g,sep=""),
			"difference"=paste(names(facts)[1],fac1g,":",
			names(facts)[2],fac2g,sep="")))
		} else { 
			gres=matrix(nrow=length(fac1g),ncol=length(fac1g),
			dimnames=list("pvalue"=fac1g,"difference"=fac1g))
		}
		for (i in 1:(length(fac1g)-1)) {
			for (j in (i+1):length(fac1g)) {
				diff=sampsg[,j]-sampsg[,i]
				gres[j,i]=mcmc.pval(diff,ptype=ptype)
				gres[i,j]=mean(diff)/log(base)
			}	
		}		
		gene.results=append(gene.results,list(gres))
		names(gene.results)[length(gene.results)]=gene		
	}

	genes=genes[!(genes %in% skips)]
	for (ge in skips) { 
		gsols=gsols[-c(grep(ge,gsols))]		
	}
	big.summary=data.frame(cbind("gene"=gsols,"f1"=fac1))
	names(big.summary)[2]=names(facts[1])
	if(nfactors==2) {
		big.summary$f2=fac2
		names(big.summary)[3]=names(facts[2])
	}
	samps=apply(samps,2, function(x){ return(x/log(base)) })
	mns=apply(samps,2,mean)
	sds=apply(samps,2,sd)
	lower=apply(samps,2,function(x) { return(quantile(x,0.05)) })
	upper=apply(samps,2,function(x) { return(quantile(x,0.95)) })
	big.summary=cbind(big.summary,"mean"=mns,
	"sd"=sds,"lower"=lower,"upper"=upper)
	ploo=NULL

	if(relative) {
		if (nfactors==2) { 
			remov=which(
				big.summary[,2]==facts[[1]][1] & big.summary[,3]==facts[[2]][1]
			)
			big.summary=big.summary[-remov,] 
		} else { 	
			big.summary=big.summary[big.summary[,2]!=facts[[1]][1],] 
			}
	}
	if (is.null(xgroup)) {
		xgroup=names(facts[1])
		if(nfactors==2) { facet=names(facts[2]) }
	} else {
		if (xgroup==names(facts[1])) {
			if(nfactors==2) { facet=names(facts[2]) }
		} else {
			facet=names(facts[1])
		}
	}

	if(summ.plot) { 
		if (!relative) { 
			if (nfactors==2) {
				ploo=summaryPlot(big.summary,xgroup=xgroup,
				facet=facet,type="line",log.base=log.base,...) 
			} else {
				ploo=summaryPlot(big.summary,xgroup=names(facts[1]),type="line",log.base=log.base,...) 
			}
		} else {
			if (nfactors==2) {
				ploo=summaryPlot(big.summary,xgroup=xgroup,
				facet=facet,type="bar",log.base=log.base,...) 
			} else {
				ploo=summaryPlot(big.summary,xgroup=names(facts[1]),type="bar",log.base=log.base,...) 
			}
		}
	}

	return(list("summary"=big.summary,"geneWise"=gene.results,"ggPlot"=ploo))
}
