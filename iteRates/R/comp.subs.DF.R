comp.subs.DF <-
function(tree,thr=6,srt="drop",min.val=0.01,mod.id=c(1,0,0,0)){

	if(length(mod.id)!=4){
		stop("\nmod.id error - should have four elements\n")
				} 
				
		mod.id<-as.logical(mod.id)		
	min.branch <- max(branching.times(tree))*min.val
	N=length(tree$edge.l)
	tree$node.label <- (length(tree$tip.label)+1):(2*length(tree$tip.label)-1)
	x <- tree$edge.length
	delta <- tree$edge[,2]>length(tree$tip.label)
	subs=subtrees2.6(tree)
	tests=matrix(NA,nrow=N/2,ncol=14)
	cat("\n1 of",N/2,"subtrees2.6\n*")
		for(i in 2:(N/2)){
			if(i%%20==0){cat("\n",i,"of",N/2,"subtrees2.6\n")}
			if(i%%20!=0){cat("*")}
				
	sZ<-length(subs[[i]]$edge.length)
		if(sZ>thr){
			if(sZ<(N-thr)){
			labs <- class.edge.all(tree,i,subroot=srt)
			lab <- class.edge.all(tree,i,subroot="drop")
		y1 <- na.omit(x[labs])
		y2 <- na.omit(x[labs==FALSE])
		y0 <- c(y1,y2)
		b1 <- na.omit(delta[labs])
		b2 <- na.omit(delta[labs==FALSE])
		b0 <- c(b1,b2)
		test1 <- model.sel.distbn(b1,y1,min.branch,mod.id)
		test2 <- model.sel.distbn(b2,y2,min.branch,mod.id)
		test0 <- model.sel.distbn(b0,y0,min.branch,mod.id)
		chi.df <- 1+test1$n+test2$n-test0$n 
		tests[i,]=cbind(test0$P1,test0$P2,test1$P1,test1$P2,test2$P1,test2$P2,test0$LL,test1$LL+test2$LL,test0$Mod, test1$Mod,test2$Mod,tree$edge[is.na(lab),1],tree$edge[is.na(lab),2],pchisq(-2*(test0$LL-test1$LL-test2$LL),chi.df,lower.tail=FALSE))
		}}		
		}
		cat("\n")
		colnames(tests) <- c("Par1.tot","Par2.tot","Par1.tr1","Par2.tr1","Par1.tr2","Par2.tr2","llk.1r","llk.2r","mod.1r.tot","mod.2r.tr1","mod.2r.tr2","node1","node2","p.val")
		testsout <- as.data.frame(tests)
		testsout
	}

