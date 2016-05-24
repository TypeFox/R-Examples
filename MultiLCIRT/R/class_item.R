class_item <- function(S,yv,k,link=1,disc=0,difl=0,fort=FALSE,disp=FALSE,tol=10^-10){
	
	# print warnign
	cat("\n * Warning: this function can take a long execution time *\n\n")
	# preliminaries
	r = dim(S)[2]
	multi0 = matrix(1:r,r,1)
	label0 = -(1:r)
	out0 = est_multi_poly(S,yv,k,link=link,disc=disc,difl=difl,multi=multi0,fort=fort,disp=disp,tol=tol)
	lk0 = out0$lk; np0 = out0$np
	lk = NULL; np = NULL; merge = NULL
	# try to join two clusters
	rmulti0 = nrow(multi0); g = 0
	while(rmulti0>2){
		g = g+1
		bestlk = -Inf
		for(j1 in 1:(rmulti0-1)) for(j2 in (j1+1):rmulti0){
			AA = matrix(multi0[-c(j1,j2),],rmulti0-2,ncol(multi0))
			multi1 = rbind(cbind(AA,matrix(0,rmulti0-2,ncol(multi0))),
	    		           c(multi0[j1,],multi0[j2,]))
			for(h in 1:nrow(multi1)){
				ve = multi1[h,]
				ve = sort(ve[ve>0])
				multi1[h,] = c(ve,rep(0,ncol(multi1)-length(ve)))
			}
			ind = which(colSums(multi1)==0)
			if(length(ind)>0) multi1 = multi1[,-ind]
			label1 = c(label0[-c(j1,j2)],g) 
			out1 = est_multi_poly(S,yv,k,link=link,disc=disc,difl=difl,multi=multi1,fort=fort,disp=disp,tol=tol)
			if(out1$lk>bestlk){
				bestj1 = j1; bestj2 = j2
				bestlk = out1$lk; bestnp = out1$np; bestmulti1 = multi1; bestlabel1 = label1
			}
		}
		lk = c(lk,bestlk); np = c(np,bestnp); multi0 = bestmulti1;
		merge = rbind(merge,label0[c(bestj1,bestj2)]) 
		label0 = bestlabel1
		rmulti0 = nrow(multi0)
	}
	g = g+1
	multi1 = 1:r
	out1 = est_multi_poly(S,yv,k,link=link,disc=disc,difl=difl,multi=multi1,fort=fort,disp=disp,tol=tol)
	lk = c(lk,out1$lk); np = c(np,out1$np); multi0 = multi1
	merge = rbind(merge,c(label0[1],label0[2])) 
	
	#create gourps
	groups = vector("list",r-1)
	groups[[1]] = -merge[1,]
	for(g in 2:(r-1)){
		if(all(merge[g,]<0)){
			a = -merge[g,1]; b = -merge[g,2]
		}else if(all(merge[g,]>0)){
			a = groups[[merge[g,1]]]; b = groups[[merge[g,2]]]
		}else{
			a = -min(merge[g,]); b = groups[[max(merge[g,])]]
		}
		if(min(a)<min(b)) groups[[g]] = c(a,b) else groups[[g]] = c(b,a)
	}

	# output
	height=-2*(lk-lk0)
	dend = list(merge=merge,order=groups[[r-1]],height=-2*(lk-lk0))
	class(dend) = "hclust"
	out = list(merge=merge,height=height,lk=lk,np=np,lk0=lk0,np0=np0,
	           groups=groups,dend=dend,call = match.call())
	class(out) = "class_item"
	out
    
}
