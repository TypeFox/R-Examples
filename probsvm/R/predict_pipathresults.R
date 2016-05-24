predict_pipathresults = function(obj, new.x=NULL, pi=NULL)
{

obj.pi=obj$pi[1]
obj.alpha0=obj$alpha0[1]
obj.alpha=obj$alpha[1,]

for (ii in 2:length(obj$pi))
	{
	if (obj$pi[ii] != obj$pi[ii-1])
		{
		obj.pi = c(obj.pi, obj$pi[ii])
		obj.alpha0 = c(obj.alpha0, obj$alpha0[ii])
		obj.alpha = rbind(obj.alpha, obj$alpha[ii,])
		}
	}

kernel = obj$kernel 

kparam =  obj$kparam

if (is.null(new.x)) {	new.x=obj$x  }

if (is.null(pi)) {pi=obj.pi}

if (class(new.x)!="matrix" & class(new.x)!="data.frame") 
{stop("The new covariates must be either a matrix or a data.frame.")}

if (class(new.x)=="data.frame")
{ new.x=as.matrix(new.x) }

if (ncol(new.x)!=ncol(obj$x)) 
{stop("The new covariates matrix has a wrong dimension.")}

if (!is.numeric(pi)) {stop("The parameter pi must be numeric.")}
if (min(pi)<0 | max(pi)>1) {stop("The parameter pi must be in [0,1].")}

K <- Kmat(new.x, obj$x, kernel, kparam)

pred.y = numeric(0)
alpha0 = numeric(0)
alpha =  numeric(0)
f.hat =  numeric(0)

for (i in 1:length(pi))	
	{
	temp=pi[i]
	index=which(obj.pi==temp)
#############################################
	if (length(index)==1)
		{
		temp.alpha=obj.alpha[index,]
		temp.alpha0=obj.alpha0[index]
		new.y1=K %*% temp.alpha + temp.alpha0
		}
###################
	if (length(index)==0)
		{
		if (temp<(obj.pi[1])) 	{temp.alpha=obj.alpha[1,]
						temp.alpha0=obj.alpha0[1]
						new.y1=K %*% (temp.alpha*obj$y) + temp.alpha0 }
		if (temp>(obj.pi[length(obj.pi)])) 
						{temp.alpha=obj.alpha[length(obj.pi),]
						temp.alpha0=obj.alpha0[length(obj.pi)]
						new.y1=K %*% (temp.alpha*obj$y) + temp.alpha0 }
		if (temp<(obj.pi[length(obj.pi)]) & temp>(obj.pi[1]))
			{
			index2=max(which(obj.pi<temp))
			temp.alpha=obj.alpha[index2,]+(temp-obj.pi[index2])/(obj.pi[(index2+1)]-obj.pi[index2]) * (obj.alpha[(index2+1),]-obj.alpha[index2,])
			temp.alpha0=obj.alpha0[index2]+(temp-obj.pi[index2])/(obj.pi[(index2+1)]-obj.pi[index2]) * (obj.alpha0[(index2+1)]-obj.alpha0[index2])
			new.y1=K %*% (temp.alpha*obj$y) + temp.alpha0
			}
		}

f.hat=cbind(f.hat,new.y1)
new.y1=as.numeric(new.y1>0)*2-1
pred.y=cbind(pred.y,new.y1)
alpha=cbind(alpha,temp.alpha)
alpha0=c(alpha0,temp.alpha0)

}

colnames(f.hat)=NULL;colnames(alpha)=NULL;colnames(pred.y)=NULL;

z=list(pi=obj.pi,fitted.alpha0=alpha0,fitted.alpha=alpha,fitted.f=f.hat,predicted.y=pred.y)

return(z)

} # main function
