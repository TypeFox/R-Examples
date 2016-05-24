HPDplotBygeneBygroup = function(model,gene,group1,group2,group3=NULL,interval="ci",colors=c("coral","cyan3","grey50"),symbols=c(19,17,15),jitter=0.16,yscale="log2",...) {
	
#model=naive2;gene="mir9";group1=a1221;group2=dmso;group3=eb;interval="ci";colors=c("coral","cyan3","grey50");symbols=c(19,17,15);jitter=0.16;yscale="log2"	
	
	a=HPDplotBygene(model,gene,conditions=group1,interval=interval,plot=F,yscale=yscale,...)
	b=HPDplotBygene(model,gene,conditions=group2,interval=interval,plot=F,yscale=yscale,...)
	lims=c(a$min,a$max,b$min,b$max)
	jit=c(0-jitter/2,0+jitter/2)
	if (is.null(group3)==FALSE) {
		c=HPDplotBygene(model,gene,conditions=group3,interval=interval,plot=F,...)
		lims=append(lims,c(c$min,c$max))
		jit=c(0-jitter,0,0+jitter)
	}
marg=0.25
if (yscale=="proportion") {
#	lims=sin(lims)^2
	marg=0.02
}

g1= HPDplotBygene(model,gene,conditions=group1,interval=interval,jitter=jit[1],pch=symbols[1],col=colors[1],ylimits=c(min(lims)-marg,max(lims)+marg),main=gene,yscale=yscale,...)
g2= HPDplotBygene(model,gene,conditions=group2,interval=interval,jitter=jit[2],pch=symbols[2],col=colors[2],newplot=F,yscale=yscale,...)
if (is.null(group3)==FALSE) {
	g3= HPDplotBygene(model,gene,conditions=group3,interval=interval,jitter=jit[3],pch=symbols[3],col=colors[3],newplot=F,yscale=yscale,...)
	return(list(g1,g2,g3))
} else { return(list(g1,g2)) }
}