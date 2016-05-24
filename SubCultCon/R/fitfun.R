fitfun <-
function(matches,nl,grp,ng){
	match1=matches[grp==1,grp==1]
	ans1=mlesol1(match1,nl)
	fitness=sum(ans1$par)
	for(ig in 2:ng){
		if(sum(grp==ig)>1){
			match2=matches[grp==ig,grp==ig]
			ans2=mlesol1(match2,nl)
			fitness=fitness+sum(ans2$par)
		}else if(sum(grp==ig)==1){
			fitness=fitness+1
		}
	}
	fitness
}
