HelpShiftsPPn <- function(t,timevec,lambdavec,muvec,rho=1,survival=1,root=1){
	out<- 0+root*(-log(length(t)))
	for (i in 2:length(t)){
		out<-out+log(Fderifuncshift(t[i],timevec,lambdavec,muvec,rho))+log(Ffuncshift(t[1],timevec,lambdavec,muvec,rho))-log(Ffuncshift(t[1],timevec,lambdavec,muvec,rho)-1)
	}
	out
}
