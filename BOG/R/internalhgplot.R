internalhgplot <-
function(stat)
{
	h=stat$plot$h
	h_ref=stat$plot$h_ref
	key=LETTERS
	e.c=values(h_ref,USE.NAMES=FALSE)
	o.c=values(h,USE.NAMES=FALSE)
	e.total=sum(e.c)
	o.total=sum(o.c)
	e.f=e.c/e.total*o.total
	o.f=o.c
	margin.y=1.5*max(e.f,o.f)
	actual_cog=stat$plot$actual_cog
	
	sorted.cog=stat$hyper[,1][stat$hyper[,3]<0.1]
	sorted.pval=stat$hyper[,3][stat$hyper[,3]<0.1]  #Adjusted P-value
	actual_cog=stat$plot$actual_cog[stat$hyper[,3]<0.1]
	n=length(actual_cog)
	ecount=rep(0,n)
	ocount=rep(0,n)
	etotal=rep(0,n)
	ototal=rep(0,n)
	
	plotmat=matrix(nrow=2,ncol=n,0)
	alternative=stat$alternative

	for(i in 1:n){
		o=values(h,key=sorted.cog[i],USE.NAMES=FALSE)
		e=values(h_ref,key=sorted.cog[i],USE.NAMES=FALSE)
		plotmat[1,i]= o 
		plotmat[2,i]= round(e/e.total*o.total,digits=0)
		ocount[i]=o
		ecount[i]=round(e/e.total*o.total,digits=0)
		ototal[i]=o.total
		etotal[i]=e.total
	}	
	label=vector()
	for(i in 1:length(sorted.cog)){
		v1=sorted.pval[i]
		if(v1<0.00001){
			v=paste(sorted.cog[i]," (adj p < 0.00001)\n","o=",ocount[i]," e=",ecount[i],sep="")
		}else{
			v=paste(sorted.cog[i]," (adj p =",round(sorted.pval[i],digits=5),")\n","o=",ocount[i]," e=",ecount[i],sep="")
		}
		label=append(label,v)
	}
	
	list(plotmat=plotmat,label=label,margin.y=margin.y,e.total=e.total,o.total=o.total)
}
