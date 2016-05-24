fitfun2 <-
function(matches,nl,grp){
	match1=matches[grp==1,grp==1]
	ans1=mlesol1(match1,nl)
	if(sum(grp==2)>1){
		match2=matches[grp==2,grp==2]
		ans2=mlesol1(match2,nl)
		fitness=sum(ans1$par)+sum(ans2$par)
	}else if(sum(grp==2)==1){
		fitness=sum(ans1$par)+1
	}else{fitness=sum(ans1$par)}
	fitness
}
