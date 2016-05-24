## DEFINE FUNCTIONS
# jump-diffusion simulator

# plotting
ex.traitgram=function(phy, hist, alpha=0, cex.node=2, scl=2, ...){
	mm=max(abs(alpha-unlist(hist$phenotype)))
	root=Ntip(phy)+1
	plot(x=NULL, y=NULL, xlim=range(pretty(c(0,max(hist$time)))), ylim=range(pretty(c(-mm+alpha, mm+alpha))), bty="n", xlab="time", ylab="phenotypic value")
	for(i in 1:nrow(hist)) {
		start=ifelse(hist$ancestor[i]==root, alpha, hist$phenotype[which(hist$descendant==hist$ancestor[i])])
		stime=ifelse(hist$ancestor[i]==root, 0, hist$time[which(hist$descendant==hist$ancestor[i])])
		
		end=hist$phenotype[i]
		etime=hist$time[i]
		lines(c(stime,etime),c(start,end),col=.transparency("gray25",0.75),...)
	}
	points(hist$time,hist$phenotype,bg=.transparency("white",0.75),pch=21,cex=ifelse(hist$descendant<=Ntip(phy), cex.node, cex.node/scl))
	points(0,alpha,bg=.transparency("white",0.75),pch=21,cex=cex.node/scl)
}

# jump simulation
ex.jumpsimulator=function(phy, alpha=0, sigmasq.brown=0.01, sigmasq.jump=1, jumps){
    
    .jumpsim=function(phy, alpha=0, sigmasq.brown=0.01, sigma.jump=0.2, lambda.jump=0.2){
        phy=reorder(phy)
        cedges=cumsum(phy$edge.length)
        tmax=sum(phy$edge.length)
        nn=phy$edge[,2]
        cs=c()
        tt=c()
        jumps=0
        if(lambda.jump>0){
            while(1){
                dt=rexp(1,lambda.jump)
                tt=c(tt,dt)
                cs=cumsum(tt)
                if(cs[length(cs)]>tmax){
                    cs=cs[-length(cs)]
                    jumps=length(cs)
                    break()
                }
            }
        }
        
        jump.edges=rep(0, length(nn))
        if(jumps>0) {
            for(j in 1:length(cs)){
                tmp=min(which(cs[j]<cedges))
                jump.edges[tmp]=jump.edges[tmp]+1
            }
        }
        
        hist=as.data.frame(matrix(cbind(phy$edge, phy$edge.length, NA, jump.edges, NA, NA, NA), ncol=8))
        names(hist)=c("ancestor","descendant","edge","phenotype","jumps","time", "effect_Brown", "effect_Jump")
        
        root=Ntip(phy)+1
        for(i in 1:nrow(hist)){
            start=ifelse(hist$ancestor[i]==root, alpha, hist$phenotype[which(hist$descendant==hist$ancestor[i])])
            stime=ifelse(hist$ancestor[i]==root, 0, hist$time[which(hist$descendant==hist$ancestor[i])])
            t=hist$edge[i]
            hist$phenotype[i]=start+(Beffect<-rnorm(1, mean=0, sd=sqrt(sigmasq.brown*t)))
            hist$effect_Brown[i]=abs(Beffect^2)/hist$edge[i]
            hist$time[i]=stime+hist$edge[i]
            if((jmp<-hist$jumps[i])>0){
                hist$phenotype[i]=hist$phenotype[i]+(Jeffect<-sum(rnorm(jmp, mean=0, sd=sigma.jump)))
            } else {
                Jeffect=0
            }
            hist$effect_Jump[i]=abs(Jeffect^2)/hist$edge[i]
        }
        
        return(hist)
    }

	lambda=jumps/sum(phy$edge.length)
	hist=.jumpsim(phy, alpha=alpha, sigmasq.brown=sigmasq.brown, sigma.jump=sqrt(sigmasq.jump), lambda.jump=lambda)
	dat=hist$phenotype[hist$descendant<=Ntip(phy)]
	names(dat)=phy$tip.label[hist$descendant[hist$descendant<=Ntip(phy)]]
	edges=hist$descendant[hist$jumps>0]
    hist$scl=hist$effect_Jump
    
    hist$cex=(hist$scl-min(hist$scl))/(max(hist$scl)-min(hist$scl))
    hist$cex=4*asin(sqrt(hist$cex))
    
	return(list(hist=hist, dat=dat, phy=phy))
}

# local rate simulation
ex.ratesimulator=function(phy, scl=64, min=4, ...){
	while(1) {
		# find an internal edge
		anc=.get.desc.of.node(Ntip(phy)+1,phy)
		branches=phy$edge[,2]
		ss=sapply(anc, function(x) length(.get.descendants.of.node(x, phy, tips=FALSE)))
		if(all(ss<min)) stop("'min' is too large")
		branches=branches[branches>Ntip(phy) & branches!=anc]
		branch=branches[sample(1:length(branches),1)]
		desc=.get.descendants.of.node(branch,phy)
		if(length(desc)>=min) break()
	}
	rphy=phy
	rphy$edge.length[match(desc,phy$edge[,2])]=phy$edge.length[match(desc,phy$edge[,2])]*scl
	e=numeric(nrow(phy$edge))
	e[match(c(branch,desc),phy$edge[,2])]=1
	cols=c("red","gray")
	dev.new()
	plot(phy,edge.col=ifelse(e==1,cols[1],cols[2]), edge.width=2, ...)
	
	mtext("expected pattern of rates")
    
	rphy
}