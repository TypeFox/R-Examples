FP.comp.subs <-
function(tree.size,na.present,sims=100,missing=0,alpha=0.05,verbose=FALSE,...){
	loop<-0
	MCres<-NA
	MCp<-NA
	while(loop<sims){
		cat("*")
	
		x<-birthdeath.tree.geiger1.3.1(1,0,taxa.stop=tree.size+missing+1)
		x<-drop.tip2.6(x,tip=as.character(tree.size+missing+1))
		if(missing>0)x<-drop.tip2.6(x,sample(as.character(1:tree.size),size=missing))
		isNA<-tree.na.Count(x)
			if(isNA==na.present){
				loop<-loop+1
				if(loop %% 10 == 0)cat("\nSimulating",loop,"of",sims)
				res<-Qcomp.subs(x)
				MCres[loop]<-sum(na.omit(res$p.val)<alpha)
				}		
		}
		if (verbose==TRUE)cat("\n\nTree size:",tree.size,"\nMissing taxa:",missing,"\nSimulated trees:",sims,"\nSignificance threshold (alpha=",alpha,"):",quantile(MCres,0.95),"\n\n\n")
		reS<-list(tree.size=tree.size,missing=missing,sims=sims,FPRthres=quantile(MCres,0.95))
	return(reS)
	}

