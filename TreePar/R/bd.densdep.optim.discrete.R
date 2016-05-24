bd.densdep.optim.discrete<-function(x,maxN,minN,muset,model=-1,rho){
minlik<-10^100
for (k in minN:maxN){
	init<-c(2,1)
	if (muset<0) {init<-c(2)}
	res<-subplex(init,LikDD,model=model,root=1,x=x,Ndec=k,minN=minN,muset=muset,sampling=rho)
	print(k)
	print(res)
	if (res$value<minlik){
		minlik<-res$value
		mini<-res
		miniN<-k
		}
	#print(miniN)
	print(mini)
}
mini$par<-c(mini$par,miniN)
res<-mini
res
}
