probsvm=function(x, y, fold=5, kernel=c("linear","polynomial","radial"), kparam=NULL, Inum=20, type="ovr", lambdas=2^(-10:10))
{

kernel=match.arg(kernel)
this.call = match.call()

if (kernel=="polynomial" & is.null(kparam)) {kparam=1;cat("Polynomial kernel is used without kparam. Set kparam = 1.\n")}

if (kernel=="radial" & is.null(kparam)) {kparam=1;cat("Radial kernel is used without kparam. Set kparam = 1.\n")}

x.train = x
y.train = as.vector(y)


if (is.na(sum(x.train)) | is.nan(sum(x.train))) {stop("There should be no NA/NaN in x.")}

if (nrow(x.train)!=length(y.train)) {stop("The dimensions of x and y do not match.")}
if (fold<2) {stop("Fold number of CV should be greater than 1.")}
if (Inum<2) {stop("The number of partitions on [0,1] should be bigger than 1.")} 
if (type != "ovr" & type != "ovo") {stop("Only ovr and ovo are available.")}

lambdas = sort(lambdas)
if (any(lambdas<=0)) {stop("The lambdas should be positive.")}

flag=0
if (is.numeric(y.train)) {flag=1}

x.train=as.matrix(x.train)

y.train=as.factor(y.train)

#################################

	tempname = levels(factor(y.train))
	
	y.temp = rep(0,length(y.train))	
	
	kk=length(tempname)

	for (ii in 1:kk){ y.temp[which(y.train %in% tempname[ii])]=ii }
	
	

if (type=="ovr")
{

best.lambda=rep(0,kk)

	for (ii in 1:kk) ## kk differe OVA methods
	{	 
	index.pos = which(y.temp==ii)
	index.neg = which(y.temp!=ii)
	x.train.pos=x.train[index.pos,]
	x.train.neg=x.train[index.neg,]
	index.pos.cv=sample(1:length(index.pos))	
	index.neg.cv=sample(1:length(index.neg))

	loglik=matrix(0,length(lambdas),fold)
		for (jj in 1:fold) ## CV
		{	 
			index.pos.tune=index.pos.cv[ ((jj-1)*ceiling(length(index.pos)/fold)+1) : min(length(index.pos),jj*ceiling(length(index.pos)/fold))]
			index.neg.tune=index.neg.cv[ ((jj-1)*ceiling(length(index.neg)/fold)+1) : min(length(index.neg),jj*ceiling(length(index.neg)/fold))]
			index.pos.train=setdiff(index.pos.cv,index.pos.tune)
			index.neg.train=setdiff(index.neg.cv,index.neg.tune)
			x.train.temp=rbind( x.train.pos[index.pos.train,],x.train.neg[index.neg.train,] )
			x.tune.temp=rbind( x.train.pos[index.pos.tune,],x.train.neg[index.neg.tune,] )
			y.train.temp=c( rep(1,length(index.pos.train)), rep(-1,length(index.neg.train)) )
			y.tune.temp=c( rep(1,length(index.pos.tune)), rep(-1,length(index.neg.tune)) )
			
				for (pp in 1:length(lambdas))
				{	
				lambda=lambdas[pp]; 
	loglik[pp,jj]=loglik.svm(x.train.temp,y.train.temp,x.tune.temp,y.tune.temp,lambda=lambda,Inum=Inum, kernel = kernel, kparam = kparam)
	 			}			
		} ## CV
 	
	result=apply(loglik,1,mean,na.rm=T)
 
	best.lambda[ii]=lambdas[ max(which(result==min(result))) ]
	 
	} ## kk differe OVA methods

} # type=ova

if (type=="ovo")
{

best.lambda=rep(0,kk^2)


	for (ii in 1:kk) ## kk choose 2 differe OVO methods
	{ 
	for (mm in 1:kk)
	{ 
	if (ii<mm)
   {
	index.pos = which(y.temp==ii)
	index.neg = which(y.temp==mm)
	x.train.pos=x.train[index.pos,]
	x.train.neg=x.train[index.neg,]
	index.pos.cv=sample(1:length(index.pos))	
	index.neg.cv=sample(1:length(index.neg))

	loglik=matrix(0,length(lambdas),fold)
		for (jj in 1:fold) ## CV
		{  
			index.pos.tune=index.pos.cv[ ((jj-1)*ceiling(length(index.pos)/fold)+1) : min(length(index.pos),jj*ceiling(length(index.pos)/fold))]
			index.neg.tune=index.neg.cv[ ((jj-1)*ceiling(length(index.neg)/fold)+1) : min(length(index.neg),jj*ceiling(length(index.neg)/fold))]
			index.pos.train=setdiff(index.pos.cv,index.pos.tune)
			index.neg.train=setdiff(index.neg.cv,index.neg.tune)
			x.train.temp=rbind( x.train.pos[index.pos.train,],x.train.neg[index.neg.train,] )
			x.tune.temp=rbind( x.train.pos[index.pos.tune,],x.train.neg[index.neg.tune,] )
			y.train.temp=c( rep(1,length(index.pos.train)), rep(-1,length(index.neg.train)) )
			y.tune.temp=c( rep(1,length(index.pos.tune)), rep(-1,length(index.neg.tune)) )

				for (pp in 1:length(lambdas))
				{
				lambda=lambdas[pp]

	loglik[pp,jj]=loglik.svm(x.train.temp,y.train.temp,x.tune.temp,y.tune.temp,lambda=lambda,Inum=Inum, kernel = kernel, kparam = kparam)

					}			
		} ## CV
 
	result=apply(loglik,1,mean,na.rm=T)
 
	best.lambda[((ii-1)*kk+mm)]=lambdas[ max(which(result==min(result))) ]	

   } ## ii<mm
	} ## mm =1:kk
	} ## ii = kk choose 2 differe OVO methods
	

} # type=ovo

z=list(x=x,y=y,fold=fold, kernel=kernel, kparam=kparam, Inum=Inum, type=type, lambdas=lambdas, best.lambda=best.lambda, call=this.call)

class(z) = "probsvm"

return(z)

} #### probsvm



