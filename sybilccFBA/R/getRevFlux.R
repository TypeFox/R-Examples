getRevFlux<-function(model,modirrev,fdirrev){
# for reversible rxns get (_f-_b) 
fluxes=NULL;
 for (r in(react_id(model))){
	if(! react_rev(model)[react_id(model)==r]){
	 fluxes=rbind(fluxes,cbind(rxn=r,fwd=fdirrev[which(react_id(modirrev)==r)],bwd=0))
	}else{
		fluxes=rbind(fluxes,cbind(rxn=r,fwd=fdirrev[which(react_id(modirrev)==paste(r,'_f',sep=""))]
	               ,bwd=fdirrev[which(react_id(modirrev)==paste(r,'_b',sep=""))]))
	
	} 
}
return(fluxes)
}
