bd.densdep.optim<-function(x,minK=0,maxK=0,discrete=TRUE,continuous=FALSE,lambdainit=2,muinit=1,Kinit=0,Yule=FALSE,muset=0,rho=1,model=-1){
	survival <- 0
	if (survival == 1){"Conditioning on survival is not yet implemented. Please wait for the next version. Method proceeds with not conditioning on survival."}
	if (Yule==TRUE){muset<- -100000}
	x<-sort(x)
	if (minK==0) {minK<-length(x)+1}
	if (maxK==0) {maxK<-round(minK*1.5)}
	resdiscrete<-0
	rescont<-0
	if (discrete ==TRUE){
	resdiscrete<-bd.densdep.optim.discrete(x,maxK,minK,muset,model=model,rho=rho)}
	if (continuous==TRUE){
		if (Kinit==0) {Kinit<-(length(x)+2)}
	init<-c(lambdainit,muinit,Kinit)
	if (muset<0) {init<- c(lambdainit,Kinit)}	
	rescont<-subplex(init,LikDD,root=1,model=model,x=x,sampling=rho,minN=minK,muset=muset,control=list(reltol=10^(-10)))}
res<-list(resdiscrete,rescont)
res
}
