predict.probsvm = function(object,new.x=NULL,...)
{
obj=object

if (!is.null(new.x)) {x.test = new.x} else {x.test=obj$x}


if (ncol(obj$x)!=ncol(x.test)) {stop("The dimensions of x and new.x do not match.")}

x.train = obj$x
y.train = obj$y

flag=0
if (is.numeric(y.train)) {flag=1}

x.train=as.matrix(x.train)
x.test=as.matrix(x.test)
best.lambda=obj$best.lambda
kernel=obj$kernel
kparam=obj$kparam
type=obj$type
Inum=obj$Inum

y.train=as.factor(y.train)

tempname = levels(factor(y.train))
	
y.temp = rep(0,length(y.train))	
	
kk=length(tempname)

for (ii in 1:kk){ y.temp[which(y.train %in% tempname[ii])]=ii }

est.prob=matrix(0,nrow(x.test),kk)

	if (type=="ovr")
	{
	
	for (ii in 1:kk) ## kk differe OVA methods
	{	 
	index.pos = which(y.temp==ii)
	index.neg = which(y.temp!=ii)

	ytemp=y.temp
	ytemp[index.pos]=1
	ytemp[index.neg]=-1

	est.prob[,ii]=prob.svm(x.train, ytemp, x.test, lambda=best.lambda[ii], Inum=Inum, kernel = kernel, kparam = kparam)
	} # ii in 1:kk
	} # if ovr

	if (type=="ovo")
	{

R=matrix(0,nrow(x.test),kk^2) ## r matrix

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
	
	ytemp=y.temp[which(y.temp==ii | y.temp==mm)]
	ytemp[which(ytemp==ii)]=1
	ytemp[which(ytemp==mm)]=-1
	x.train.ovo=x.train[which(y.temp==ii | y.temp==mm),]

R[,((ii-1)*kk+mm)]=prob.svm(x.train.ovo,ytemp,x.test,lambda=best.lambda[((ii-1)*kk+mm)],Inum=Inum, kernel = kernel, kparam = kparam); 

R[,((mm-1)*kk+ii)]=1-R[,((ii-1)*kk+mm)]

   } ## ii<mm
	} ## mm =1:kk
	} ## ii = kk choose 2 differe OVO methods


for (ii in 1:nrow(x.test)) 
	{
	### Q matrix in paper svmprob
	Q=matrix(0,kk,kk)
		for (pp in 1:kk)
		{
			for (qq in 1:kk)
			{
			if (pp==qq) {
						for (gg in setdiff(1:kk,pp))
						{Q[pp,qq]=Q[pp,qq]+R[ii,((gg-1)*kk+pp)]^2}
					}
			if (pp!=qq) {Q[pp,qq]= -R[ii,((pp-1)*kk+qq)]*R[ii,((qq-1)*kk+pp)]}
			}
		}
	### Q matrix
	p.old=rep(1/kk,kk)
	p.new=p.old
	for (iter in 1:20)
		{
			for (tt in 1:kk)
			{
			p.new[tt]=1/Q[tt,tt] * ( -sum((Q[,tt]*p.new)[-tt])+ as.numeric(p.new%*%Q%*%p.new))
			p.new=p.new/(sum(p.new)) 
			}
			
		if (any(is.nan(p.new)) | any(is.na(p.new)) | any(p.new==Inf)) {p.new=p.old; break}
		if (sum(abs(p.new-p.old))<0.01) {break}
		p.old=p.new
		}
	p.new=p.new/(sum(p.new))
	est.prob[ii,]=p.new
	} ## ii 1:nrow(x.test)

	} # if ovo

pred.y=apply(est.prob,1,pred.cz)

pred.y2=character(length(pred.y))

for (ii in 1:kk)
{
pred.y2[pred.y==ii] = tempname[ii]
}

colnames(est.prob)=tempname

est.prob[est.prob<=0] = 1e-6

pred.prob=est.prob/apply(est.prob,1,sum)

if (flag==1)
{
pred.y2 = as.numeric(pred.y2)
}

return(list(object=obj,new.x=x.test,pred.prob=pred.prob,pred.y=pred.y2))

} # predict function






