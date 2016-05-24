GAmaxcomp <-
function(xmat,ng,npop,ngen){
	n=dim(xmat)[1];m=dim(xmat)[2];nl=max(xmat)
#  compute match matrix
	matches=matrix(nrow=n,ncol=n)
	for(i in 1:n){
		for(j in 1:n){
			matches[i,j]=sum(xmat[i,]==xmat[j,])
		}
	}
	keymat=matrix(nrow=m,ncol=ng)
	dvec=1:n
### initialize population for the genetic algorithm
	np1=min(n/2-10,npop/2)
	pans=1:nl;obs4=1:nl
	key=matrix(nrow=m,ncol=2)
	popmat=matrix(1,nrow=npop,ncol=n)
	gord=1:ng;gnum=1:ng;newpop=1:n
## random initialization
	for(ipop in 1:npop){
		popmat[ipop,]=trunc(runif(n)*ng+1)
		for(i in 1:ng){gnum[i]=sum(popmat[ipop,]==i)}
		gord=order(gnum)
		if(!all(gord==ng:1)){
			for(i in 1:ng){newpop[popmat[ipop,]==gord[i]]=ng+1-i}
			popmat[ipop,]=newpop
		}
	}
	fitness=1:npop
## compute fitnesses for the initial population
	for(ipop in 1:npop){
		fitness[ipop]=fitfun(matches,nl,popmat[ipop,],ng)
	}
##   Now loop through generations					
	nrep=0
	obs=1:npop
	q1=trunc(npop*.2)
#	q3=trunc(npop*.95)
	nbig=1000;nm=trunc(npop/10);check=TRUE;nm2=nm
	score=1:ng;og=1:ng
	ord=order(-fitness)
	popmat=popmat[ord,]
	fitness=fitness[ord]
	lastbest=popmat[n,]
	while(nrep<ngen&check){
		nrep=nrep+1
		if(nrep==5){nm2=trunc(nm2/2)+1}
		if(nrep==10){nm2=trunc(nm2/2)+1}
		if(nrep==20){nm2=1}
		print(cbind(fitness[q1],fitness[1]))
		print(popmat[1,])
		qc=round(q1*1.25);fitcut=fitness[qc]
#  mutate!  randomly throughout population  (need to recompute fitness)
		imut=trunc(runif(nm)*(npop-1)+2)  ## don't mutate the best fit
		igene=trunc(runif(nm)*n+1)
		for(im in 1:nm){
			popmat[imut[im],igene[im]]=popmat[imut[im],igene[im]]+1
			if(popmat[imut[im],igene[im]]>ng){popmat[imut[im],igene[im]]=1}
			ipop=imut[im]
			for(i in 1:ng){gnum[i]=sum(popmat[ipop,]==i)}
			gord=order(gnum)
			if(!all(gord==ng:1)){
				for(i in 1:ng){newpop[popmat[ipop,]==gord[i]]=ng+1-i}
				popmat[ipop,]=newpop
			}
			fitness[ipop]=fitfun(matches,nl,popmat[ipop,],ng)
		}
		ord=order(-fitness)
		popmat=popmat[ord,]
		fitness=fitness[ord]
#  mutate!  non-randomly using mike's switching  (need to recompute fitness)
		imut=trunc(runif(nm2)*(npop-1)+2)  ## don't mutate the best fit
		for(im in 1:nm2){
#  get the answer key for the chosen groups
			try=popmat[imut[im],]
			if(sum(try==ng)>3){
			for(ig in 1:ng){
				anssub=mlesol1(matches[try==ig,try==ig],nl)
				dvec[try==ig]=anssub$par
				pans=1:nl;obs4=1:nl
				for(i in 1:m){	
					for(g in 1:ng){
						for(l in 1:nl){pans[l]=sum(dvec[xmat[try==g,i]==l])}
						keymat[i,g]=min(obs4[pans==max(pans)])
					}
				}
			}
## for each informant, swap groups if answers closer to another key	
			for(ii in 1:n){
				for(ig in 1:ng){score[ig]=sum(xmat[ii,]==keymat[,ig])}
				if(score[try[ii]]<max(score)){try[ii]=min(og[score==max(score)])}								
			}	
			for(i in 1:ng){gnum[i]=sum(try==i)}
			gord=order(gnum)
			if(!all(gord==ng:1)){
				for(i in 1:ng){newpop[try==gord[i]]=ng+1-i}
				try=newpop
			}
			newfit=fitfun(matches,nl,try,ng)
			if(newfit>fitness[imut[im]]){
				popmat[imut[im],]=try
				fitness[imut[im]]=newfit
			}
			}
		}
# reproduction : replace lower 3/4 with offspring -- combine elite with other  
		ord=order(-fitness)
		popmat=popmat[ord,]
		fitness=fitness[ord]
		for(ipop in q1:npop){
			mom=trunc(runif(1)*npop+1)
			dad=trunc(runif(1)*q1+1)
			digits=runif(n)>.6   ### more genes from (elite) father!
			popmat[ipop,digits]=popmat[mom,digits]
			popmat[ipop,!digits]=popmat[dad,!digits]
			for(i in 1:ng){gnum[i]=sum(popmat[ipop,]==i)}
			gord=order(gnum)
			if(!all(gord==ng:1)){
				for(i in 1:ng){newpop[popmat[ipop,]==gord[i]]=ng+1-i}
				popmat[ipop,]=newpop
			}
# find fitness of baby
			fitness[ipop]=fitfun(matches,nl,popmat[ipop,],ng)
		}
		ord=order(-fitness)
		popmat=popmat[ord,]
		fitness=fitness[ord]
# try targeted mutation
		if(!all(popmat[1,]==lastbest)){
			lastbest=popmat[1,]
			replace=0
			for(i in 1:min(n,npop/ng-1)){
				try=lastbest
				for(j in 1:(ng-1)){
					try[i]=try[i]+1;if(try[i]>ng){try[i]=1}
					fittry=fitfun(matches,nl,try,ng)
					if(fittry>fitness[1]){
						popmat[npop-replace,]=try;fitness[npop-replace]=fittry
						replace=replace+1
					}
				}
			}
			if(replace>0){
				ord=order(-fitness)
				popmat=popmat[ord,]
				fitness=fitness[ord]
			}
		}
		if((fitness[1]-fitness[q1])/fitness[1]<1e-5){check=FALSE}
	}	
	keymat=matrix(nrow=m,ncol=ng)
	pans=1:nl;obs4=1:nl
	bestsplit=popmat[1,]
	compvec=1:n
	for(g in 1:ng){
		if(sum(bestsplit==g)>1){
			mle=mlesol1(matches[bestsplit==g,bestsplit==g],nl)
			compvec[bestsplit==g]=mle$par
		}else if(sum(bestsplit==g)==1){
			compvec[bestsplit==g]=.99995
		}
	}
	for(i in 1:m){	
		for(g in 1:ng){
			for(l in 1:nl){pans[l]=sum(compvec[xmat[bestsplit==g,i]==l])}
			keymat[i,g]=min(obs4[pans==max(pans)])
		}
	}
	ans=new.env()
	ans$compsum=fitness[1]
	ans$bgrp=bestsplit
	ans$numgen=nrep
	ans$keymat=keymat
	ans$comp=compvec
	ans	
}
