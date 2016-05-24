startv <-function(dat,m,method='both'){
if(method=='survival') return(log(sum(dat$cdeath)/sum(dat$tdeath))*rep(1,m))
else { 
	na=is.na(dat$tprog1);
	tpr=ifelse(na,dat$tprog0,(dat$tprog0+dat$tprog1)/2)
	dpr=as.numeric(!na)
	lam1=log(sum(dpr)/sum(tpr))
	if (method=='progression') return ( lam1*rep(1,m))
	else {

		td=ifelse(na,0,dat$tdeath-tpr)
		lam2=log(sum(dat$cdeath)/sum(td))
		return(c(lam2*rep(1,m*(m+1)/2),lam1*rep(1,m)))
		}
	}
}
