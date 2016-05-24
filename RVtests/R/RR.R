RR <-
function(x, y, z = NULL, scale = FALSE, weights = 1, lambda = 1, 
		npermutation = 1000, npermutation.max, min.nonsignificant.counts = 100)
{
#RR: ridge regression
#The adaptive permutation is implemented by blocks
#The number of blocks is ceiling(npermutation.max/npermutation)
#x: a matrix of genotypes
#y: a vector of phenotypes
#z: a matrix of covariates if including
#always centralizing x, y, and z
x<- x*rep(weights, each = nrow(x))
x<- scale(x, scale = scale)
y<- y-mean(y)
if (!is.null(z)) z<- scale(z, scale = scale)

if (missing(npermutation.max)) npermutation.max<- npermutation

#x=UDV': svd(x)
#svd has a bug if not use LINPACK!
UDV<- function(x){
	#a<- svd(x)
	useLINPACK = TRUE #but first try LAPACK.
	try({a<- svd(x); useLINPACK=FALSE}, silent = FALSE)
	if(useLINPACK) {a<- svd(x, LINPACK=TRUE); message("Note: the error due to LAPACK is disappeared using LINPACK.")}

	d<- a$d
	indx<- (abs(d) > 1e-8)  #remove zero eigenvalues.
	list(u=a$u[,indx, drop=FALSE], d=d[indx], v=a$v[,indx, drop=FALSE])
	}
	
UDVxz<- UDV(cbind(x, z))
if (!is.null(z)) UDVz<- UDV(z)

###############
#score function
#sd(y)=sqrt(sum((y-ym)^2)/(length(y)-1))
RRscore<- function(yp, UDVx){
	u<- UDVx$u
	d<- UDVx$d
	yp<- as.matrix(yp)  #yp: each column coresponding to a sample of y.
	n<- nrow(yp)
	ypsd<- apply(yp, 2, sd)
	uy<- t(u)%*%yp  #a matrix of p by m
	score<- matrix(NA, ncol(yp), length(lambda)) 
	colnames(score)<- paste("RR", lambda, sep="")
	for (i in 1:length(lambda)){
		lam<- lambda[i]
		ye<- u%*%((d^2/(d^2+lam))*uy)
		yesd<- apply(ye, 2, sd)
		score[, i]<- (colSums(yp*ye) - colSums(yp)*colSums(ye)/n)/(ypsd*yesd)/(n-1)
		}#i
	return(score)
	}


###############
#test score and pvalue (nominal pvalue)
{if (!is.null(z)) testscore<- abs(RRscore(y, UDVxz) - RRscore(y, UDVz))
else testscore<- abs(RRscore(y, UDVxz))}

###############
#permutation scores and pvalue (empirical pvalue)
#permutation
ys<- apply(matrix(rep(y, npermutation), nrow=length(y), ncol=npermutation), 2, sample)
{if (!is.null(z)) permscore<- abs(RRscore(ys, UDVxz) - RRscore(ys, UDVz))
else permscore<- abs(RRscore(ys, UDVxz))}
counts<- apply(permscore >= rep(testscore, each = npermutation), 2, sum)
jblock<- 1

#adaptive permutation by blocks
nblocks<- ceiling(npermutation.max/npermutation)

while ((jblock < nblocks) & (min(counts) < min.nonsignificant.counts)){
jblock<- jblock + 1
	ys<- apply(matrix(rep(y, npermutation), nrow=length(y), ncol=npermutation), 2, sample)
	{if (!is.null(z)) permscore<- abs(RRscore(ys, UDVxz) - RRscore(ys, UDVz))
	else permscore<- abs(RRscore(ys, UDVxz))}
	counts<- counts + apply(permscore >= rep(testscore, each = npermutation), 2, sum)
	}#while

total.permutation<- npermutation*jblock
permpvalue<- (1+ counts)/(1 + total.permutation)


list(nonsignificant.counts = counts, total.permutation = total.permutation, score = testscore, pvalue.empirical = permpvalue, pvalue.nominal = NA)
}#RR function

#tmp<- proc.time();RR(x,y,lambda=0:4, npermutation=1000); proc.time()-tmp
#
