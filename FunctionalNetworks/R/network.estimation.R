network.estimation <-
function(Nsim,Burn.in,data.generation.obj,bic.generation.obj)
{
	if(Nsim<=Burn.in)
		stop("Number of simulations most strictly be bigger than the burn in")
	nombres.data=c("gene.data","set.data","affy.loc","gene2set.mat","set2gene.mat","Set.obj","Set.src","G.obj","G.src")	
	nombres=names(data.generation.obj)
	aux=intersect(nombres.data,nombres)
	if(length(aux)!=length(nombres.data))
		stop("data.generation.obj does not appear to be a data.generation object. Please check")
	nombres.bic=c("BIC.Gene.0.Pred","BIC.Gene.1.Pred","BIC.Set.0.Pred","BIC.Set.1.Pred")
	nombres=names(bic.generation.obj)
	aux=intersect(nombres.bic,nombres)
	if(length(aux)!=length(nombres.bic))
		stop("bic.generation.obj does not appear to be a bic.generation object. Please check")	
	
Red.C=function(fcnets,nrfcnet,gnets,nrgnet,procagens,nrpag,ncpag,gensaproc,nrgap,ncgap,data,nrdata,ncdata,datafc,nrdatafc,ncdatafc,dondeaffy,tamaffy,day,tamday,sim,burnin,matResG,matResFC,bicresgen,bicres,ssrg,ssrfc,bicgen0,bicgen1,bicfc0,bicfc1)
{
	output=.C("Red",as.integer(fcnets),as.integer(nrfcnet),as.integer(gnets),as.integer(nrgnet),as.integer(procagens),as.integer(nrpag),as.integer(ncpag),as.integer(gensaproc),as.integer(nrgap),as.integer(ncgap),as.double(data),as.integer(nrdata),as.integer(ncdata),as.double(datafc),as.integer(nrdatafc),as.integer(ncdatafc),as.integer(dondeaffy),as.integer(tamaffy),as.integer(day),as.integer(tamday),as.integer(sim),as.integer(burnin),as.integer(matResG),as.integer(matResFC),as.double(bicresgen),as.double(bicres),as.double(ssrg),as.double(ssrfc),as.double(bicgen0),as.double(bicgen1),as.double(bicfc0),as.double(bicfc1))
	    resultado1=matrix(output[[23]],nrgnet,nrgnet)
	    resultado2=matrix(output[[24]],nrfcnet,nrfcnet)
		
		resultado3=matrix(output[[25]],sim,1)  #bic gene
		resultado4=matrix(output[[26]],sim,1)  #bic fc

		resultado5=matrix(output[[27]],sim,1)   # ssr gen
		resultado6=matrix(output[[28]],sim,1)    # ssr fc
		
	    resultado=list(resultado1,resultado2,resultado3,resultado4,resultado5,resultado6)
	    nombres=c("Gene.network","Set.network","BIC.gene","BIC.set","RSS.gene","RSS.set")
	    names(resultado)=nombres
		return(resultado)
}
		bic.res=rep(0,Nsim)
	bic.res.gen=rep(0,Nsim)
	ssr.gen=rep(0,Nsim)
	ssr.fc=rep(0,Nsim)
	mat.res=matrix(0,length(bic.generation.obj$BIC.Set.0.Pred),length(bic.generation.obj$BIC.Set.0.Pred))
	mat.res.gen=matrix(0,length(bic.generation.obj$BIC.Gene.0.Pred),length(bic.generation.obj$BIC.Gene.0.Pred))
	pheno=1:dim(data.generation.obj$set.data)[2]
	
cosa=Red.C(data.generation.obj$Set.src,length(data.generation.obj$Set.src),data.generation.obj$G.src,length(data.generation.obj$G.src),data.generation.obj$set2gene.mat,dim(data.generation.obj$set2gene.mat)[1],dim(data.generation.obj$set2gene.mat)[2],data.generation.obj$gene2set.mat,dim(data.generation.obj$gene2set.mat)[1],dim(data.generation.obj$gene2set.mat)[2],data.generation.obj$gene.data,dim(data.generation.obj$gene.data)[1],dim(data.generation.obj$gene.data)[2],data.generation.obj$set.data,dim(data.generation.obj$set.data)[1],dim(data.generation.obj$set.data)[2],data.generation.obj$affy.loc,length(data.generation.obj$affy.loc),pheno,length(pheno),Nsim,Burn.in,mat.res.gen,mat.res,bic.res.gen,bic.res,ssr.gen,ssr.fc,bic.generation.obj$BIC.Gene.0.Pred,bic.generation.obj$BIC.Gene.1.Pred,bic.generation.obj$BIC.Set.0.Pred,bic.generation.obj$BIC.Set.1.Pred)
	
nombres.proc=names(data.generation.obj$set2gene.mat)
nombres.gens=names(data.generation.obj$gene2set.mat)	

results=list(cosa,nombres.proc,nombres.gens)
nombres=c("Algorithm.results","Set.names","Gene.names")	
names(results)=nombres


return(results)
}
