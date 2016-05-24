fdrTable <-
function(real.p,sim.p) {
	real.p=real.p[!is.na(real.p)]
	sim.p=sim.p[!is.na(sim.p)]
	ps=sort(unique(real.p))
	if (ps[1]==0) { minp=ps[2]} else { minp=ps[1] }
	if (minp<1e-10) { minp=1e-10 }
	sim.p=-log(sim.p+minp,10)
	real.p=-log(real.p+minp,10)
	fdr=c();dr=c()
	for(p in seq(max(real.p),0.1,-0.1)) {
		dr=append(dr,table(real.p>=p)[2])
		fdr=append(fdr,table(sim.p>=p)[2])
	}
	fdr[is.na(fdr)]=0
	FDR=round(fdr/dr,2)
	tabs=cbind("logp"=seq(max(real.p),0.1,-0.1),dr,fdr,FDR)
	row.names(tabs)=NULL
	fdrs=data.frame(tabs)
	return(fdrs)
}
