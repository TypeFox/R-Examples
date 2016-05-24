#	library(foreign)
#	library(msm)
#	library(MASS)
#	library(VGAM)
#	library(distr)
#	library(truncnorm)
#	library(Rcpp)
#	library(RcppArmadillo)
#	library(inline)
#	library(RcppEigen)
	expit<-function(x) exp(x)/(1+exp(x))

#############################
#############################
#############################

sfa<-function(M, missing.mat=NULL, gibbs=100, burnin=100, max.optim=50, thin=1, save.curr="UDV_curr", save.each=FALSE, thin.save=25, maxdim=NULL){

cat("Step 1 of 4: Formatting Data","\n")

	votes.mat<-M
	if(length(missing.mat)==0) missing.mat<-votes.mat*0
	empir.cdf=NULL
#source('FunctionsInternal_Count_EM.R')
	test.gamma.pois<-test.gamma.pois_EM
	make.Z<-make.Z_EM
	update_UDV<-update_UDV_EM

#votes.mat=dtm2;missing.mat=dtm2*0;max.loop=500; burnin=50; thin=1;empir.cdf=NULL;cutoff.seq=seq.vec;save.curr="UDV_curr";save.each=FALSE;thin.save=25;max.optim=300

		dubcent<-function(X){
		colmean.vec<-colMeans(X,na.rm=T)
		rowmean.vec<-rowMeans(X,na.rm=T)
		
		for(i in 1:dim(X)[1]) for (j in 1:dim(X)[2]) 
		X[i,j]<-X[i,j]-rowmean.vec[i]-colmean.vec[j]
		
		X-mean(X)
		}

#if(cutoff.seq==0){
dtm2<-votes.mat
unique.dtm<-sort(unique(as.vector(dtm2)))

seq.vec<-c(-1:(max(dtm2)+2))
for(i in 2:length(seq.vec)) ifelse(seq.vec[i]%in%unique.dtm, seq.vec[i]<-seq.vec[i-1]+1,seq.vec[i]<-seq.vec[i-1] )

cutoff.seq<-seq.vec
rm(seq.vec)
#}




	n<-dim(votes.mat)[1]
	k<-dim(votes.mat)[2]
	##For only fitting first few dimensions; not in this version
	if(length(maxdim)==0) maxdim<-min(n,k) else maxdim<-min(n,k,maxdim)
	
	
	
	votes.mat<-as.matrix(votes.mat)
	
	row.type<-apply(votes.mat,1,FUN=function(x)  (max(x)<=1))
	row.type[row.type==TRUE]<-"bin"
	row.type[row.type==FALSE]<-"count"
	row.type.ord<-apply(votes.mat,1,FUN=function(x)  (min(x)<0))
	row.type[row.type.ord==TRUE]<-"ord"

tau.ord<-1
if(min(votes.mat)<0) {
	ord.all<-(votes.mat[row.type=="ord",])[missing.mat[row.type=="ord"]==0]
	tau.ord<-qnorm(mean(ord.all>-1))-qnorm(mean(ord.all==0))
}

	cat('Distribution of observed data\n')
	print(table(row.type))
	sigma<-1

##Initialize D, lambda, Theta, Z
	votes.mat.init<-t(apply(votes.mat,1,FUN=function(x) x/sd(x)))
	votes.mat.init[!is.finite(votes.mat.init)]<-0

	mu.r<-rowMeans(sigma*votes.mat.init)
	mu.c<-colMeans(sigma*votes.mat.init)
	mu.grand<- mean(sigma*votes.mat.init)
	
	ones.r<-rep(1,k)
	ones.c<-rep(1,n)
	
	
	svd.mat<- sigma*votes.mat.init - ones.c%*%t(mu.c)-mu.r%*%t(ones.r)+mu.grand
	svd.mat<- svd.mat-ones.c%*%t(colMeans(svd.mat))-rowMeans(svd.mat)%*%t(ones.r)+mean(svd.mat)
	
	
	Dtau<-rep(1, dim(votes.mat)[1])
	Theta.last<-Z.next<-Z.last<-0*votes.mat
	Theta.last<-Z.last<-U.last<-0*votes.mat

##Containers
	U.run<-U.run.2<-matrix(NA, nrow=dim(votes.mat)[1], ncol=gibbs+burnin)
	V.run<-V.run.2<-matrix(NA, nrow=dim(votes.mat)[2], ncol=gibbs+burnin)
	D.run<-D.mean.run<-matrix(NA, nrow=maxdim, ncol=gibbs+burnin)
	accept.run<-NULL
	lambda.run<-NULL
	dev.run<-NULL
	V.average<-U.average<-0
	Z.run<-Z.mode.run<-0
	if(sum(row.type=="bin")>2){
	Theta.last[row.type=="bin",][votes.mat[row.type=="bin",]==1]<-dnorm(0)/pnorm(.5)
	Theta.last[row.type=="bin",][votes.mat[row.type=="bin",]==0]<- -dnorm(0)/pnorm(.5)
	Theta.last[row.type=="bin",][missing.mat[row.type=="bin",]==1]<- 0}
		dubcent<-function(X){
		colmean.vec<-colMeans(X,na.rm=T)
		rowmean.vec<-rowMeans(X,na.rm=T)
		
		for(i in 1:dim(X)[1]) for (j in 1:dim(X)[2]) 
		X[i,j]<-X[i,j]-rowmean.vec[i]-colmean.vec[j]
		
		X-mean(X)
		}

	if(sum(row.type=="count")>2){
	Theta.last[row.type=="count",]<- (log(1*(votes.mat[row.type=="count",]==0)+votes.mat[row.type=="count",]))
	#Theta.last[row.type=="count",]<-(Theta.last[row.type=="count",]-mean(Theta.last[row.type=="count",]))/sd(Theta.last[row.type=="count",])
	Theta.last[row.type=="count",][missing.mat[row.type=="count",]==1]<- 0
	#Theta.last[row.type=="count",]<- 0
	}


	if(sum(row.type=="ord")>2){
	Theta.last[row.type=="ord",][votes.mat[row.type=="ord",]==-2]<-dnorm(0)/pnorm(.5)
	Theta.last[row.type=="ord",][votes.mat[row.type=="ord",]==0]<- -dnorm(0)/pnorm(.5)
	Theta.last[row.type=="ord",][votes.mat[row.type=="ord",]==-1]<- 0

	Theta.last[row.type=="ord",][missing.mat[row.type=="ord",]==1]<- 0}


	
	svd0<-svd(svd.mat,nu=maxdim,nv=maxdim)
	UDV.all<-NULL
	UDV.all$V.next<-svd0$v
	lambda.shrink<-lambda.lasso<-rep(median((2/(svd0$d))^.5),maxdim)[1]


cat("Step 2 of 4: Beginning Burnin Period","\n")


	params.init<-c(0,0)#log(c(mean(votes.mat)))
	proposal.sd<-0.1
	step.size<-scale.sd<-1
#############################
#############################
#############################		
###########Start loop here
for(i.gibbs in 1:(gibbs+burnin)){
	#print(i.gibbs)
	#Switch off EM
	if(i.gibbs==max.optim) {
		test.gamma.pois<-test.gamma.pois_gibbs
		make.Z<-make.Z_gibbs
		update_UDV<-update_UDV_gibbs		
#		source('FunctionsInternal_Count.R')
	}

	#Step 1: Make Z-star
	Z.call<-make.Z(Theta.last.0=Theta.last, votes.mat.0=votes.mat, row.type.0=row.type,  n0=n ,k0=k, iter.curr=i.gibbs,empir=empir.cdf,cutoff.seq.0=cutoff.seq,missing.mat.0=missing.mat,lambda.lasso=lambda.lasso,params=params.init,proposal.sd=proposal.sd,scale.sd=scale.sd,max.optim=max.optim,step.size=step.size,maxdim.0=maxdim,tau.ord.0=tau.ord)

	Z.next<-Z.call$Z
	params.init<-Z.call$params.init
	accept.run[i.gibbs]<-Z.call$accept
	dev.run[i.gibbs]<-Z.call$prob
	params.init<-Z.call$params
	proposal.sd<-Z.call$proposal.sd
	step.size<-Z.call$step.size
	tau.ord<-Z.call$tau.ord
	#Step 2: Update UDV and components therein
	UDV.all<-	update_UDV(
	Z.next.0=Z.next,
	k0=k, n0=n, lambda.lasso.0=lambda.lasso,lambda.shrink.0=lambda.shrink,
	votes.mat.0=votes.mat, Dtau.0=Dtau, iter.curr=i.gibbs, row.type.0=row.type, missing.mat.0=missing.mat,maxdim.0=maxdim, V.last = UDV.all$V.next
	)
	
	last.Z<-Z.call
	last.UDV<-UDV.all
	#print(tau.ord)
	
	if(i.gibbs>(max.optim+50)&i.gibbs%%50==1){
		seq.check<-(i.gibbs-1):(i.gibbs-20)
		#print(accept.run)
		if(mean(accept.run[(i.gibbs-1):(i.gibbs-20)])>.9 ) step.size<-step.size*1.25
		if(mean(accept.run[(i.gibbs-1):(i.gibbs-20)])<.1 ) step.size<-step.size*.8
		#print("acceptance percentage")
		#print(mean(accept.run[(i.gibbs-1):(i.gibbs-20)]))
	}
	

	if(save.each){
          file.temp<-paste(save.curr,i.gibbs)
          if(i.gibbs%%thin.save ==1) {save(UDV.all,file=file.temp)}
	} else {
          ## file.temp<-paste("/scratch/insong/sfa_output/", save.curr,i.gibbs, sep="")
          file.temp <- paste(save.curr)
          save(UDV.all,file=file.temp)
	}
	#Update inputs to UDV.all
	Dtau=UDV.all$Dtau
	lambda.lasso<-UDV.all$lambda.lasso
	lambda.shrink<-UDV.all$lambda.shrink
	Theta.last<-UDV.all$Theta.last
	#Report output
	#if(i.gibbs%%100==0){
	#print(table(sign(Theta.last),votes.mat))
	#print("range theta")
	#print(range(Theta.last))
	#print("iter")
	#print(i.gibbs)
	#print("rms")
	#print(mean((Theta.last-Z.next)^2))
	#print(c("sigma",sigma))
	#print("trunc")
	#print(UDV.all$D.trunc)
	#print("post")
	#print(UDV.all$D.post[1:length(UDV.all$D.trunc)])
	#print("lambda.mean")
	#print(UDV.all$lambda.lasso)
	#print("lambda.mode")
	#print(UDV.all$lambda.shrink)
	#}##Closes printing of results
	
	
	if(i.gibbs==1) cat("    0% of Burnin Completed:", i.gibbs,"/",burnin,"\n")	
	if(i.gibbs==floor(burnin/4)) cat("   25% of Burnin Completed:", i.gibbs,"/",burnin,"\n")
	if(i.gibbs==floor(burnin/2)) cat("   50% of Burnin Completed:", i.gibbs,"/",burnin,"\n")
	if(i.gibbs==floor(3*burnin/4)) cat("   75% of Burnin Completed:", i.gibbs,"/",burnin,"\n")
	if(i.gibbs==floor(burnin)) cat("  100% of Burnin Completed:", i.gibbs,"/",burnin,"\n")
	if(i.gibbs == burnin+1) cat("Step 3 of 4: Burnin Completed; Saving Samples \n")
	if(i.gibbs==floor(burnin+1)) {
		cat("    0% of Samples Gathered:", i.gibbs-burnin,"/",gibbs,"\n")
		cat("Current dimensionality estimate:\n")
		print(UDV.all$D.trunc[1:maxdim])
		}
	if(i.gibbs==floor(burnin+(gibbs)/4)) {
		cat("   25% of Samples Gathered:", i.gibbs-burnin,"/",gibbs,"\n")
		cat("Current dimensionality estimate:\n")
		print(UDV.all$D.trunc[1:maxdim])
		}
	if(i.gibbs==floor(burnin+(gibbs)/2)) {
		cat("   50% of Samples Gathered:", i.gibbs-burnin,"/",gibbs,"\n")
		cat("Current dimensionality estimate:\n")
		print(UDV.all$D.trunc[1:maxdim])
		}
	if(i.gibbs==floor(burnin+3*(gibbs)/4)) {
		cat("   75% of Samples Gathered:", i.gibbs-burnin,"/",gibbs,"\n")
		cat("Current dimensionality estimate:\n")
		print(UDV.all$D.trunc[1:maxdim])
		}
	if(i.gibbs==floor(burnin+gibbs)) {
		cat("  100% of Samples Gathered:", i.gibbs-burnin,"/",burnin,"\n")
		cat("Current dimensionality estimate:\n")
		print(UDV.all$D.trunc[1:maxdim])
		}

##Gather results
	lambda.run[i.gibbs]<-UDV.all$lambda.lasso[1]
	U.run[,i.gibbs]<-UDV.all$U.next[,1]
	U.run.2[,i.gibbs]<-UDV.all$U.next[,2]
	V.run[,i.gibbs]<-UDV.all$V.next[,1]
	V.run.2[,i.gibbs]<-UDV.all$V.next[,2]
	D.run[,i.gibbs]<-UDV.all$D.trunc[1:maxdim]
	D.mean.run[,i.gibbs]<-UDV.all$D.post[1:maxdim]
	
	
	post.process.inner<-function(x,y){
			for(i.cor in 1:ncol(x)) y[,i.cor]<-y[,i.cor]*sign(cov(x[,i.cor],y[,i.cor]))
			return(y)
	}
	if(i.gibbs==burnin){
		U.order<-UDV.all$U.next
		V.order<-UDV.all$V.next
		U.average<-V.average<-0
	}

	if(i.gibbs>burnin & i.gibbs %% thin == 0) {
		Z.run<- Z.run + Theta.last
		Z.mode.run<- Z.mode.run + UDV.all$Theta.mode
		U.average<-U.average+post.process.inner(U.order,UDV.all$U.next)
		V.average<-V.average+post.process.inner(V.order,UDV.all$V.next)
		}
	#if(i.gibbs %%100 ==0) print(i.gibbs)
}#Closes out for loop
	cat("Step 4 of 4: Gathering Output\n")
	post.process<-function(x){apply(x,2,FUN=function(x2) x2*sign(cor(x2,x[,1])) ) }
	V.run<-post.process(V.run[,-c(1:burnin)])
	V.run.2<-post.process(V.run.2[,-c(1:burnin)])
	U.run <-post.process(U.run[,-c(1:burnin)])
	U.run.2 <-post.process(U.run.2[,-c(1:burnin)])


out<-list("dim.sparse"=D.run[,-c(1:burnin)], "dim.mean"=D.mean.run[,-c(1:burnin)], "rowdim1"=U.run, "rowdim2"=U.run.2, 
"coldim1"=V.run, "coldim2"=V.run.2, 
"lambda.lasso"=lambda.run[-c(1:burnin)], "Z" = Z.run/gibbs, "Z.mode"= Z.mode.run/gibbs,
"rowdims.all"=U.average/gibbs, "coldims.all"=V.average/gibbs)
class(out)<-"sfa"
invisible(out)
}
