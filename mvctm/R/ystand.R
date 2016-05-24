ystand <-
function(fixed,clustername,data,method,y,N,p,nlevel,stand,affequiv){

# if(p>1){library("MNM")}
# if(method=="mixed"){library("nlme")}
# else if(method=="rank" & p==1){library("Rfit")}
# else if(method=="sign" & p==1){library("quantreg")}

sv=diag(p) 	

if(method=="ls"){method="identity"}
if(method=="identity" & p>1){affequiv=TRUE}

if(method=="mixed")
	{
	rando="~1|"
	for(i in 1:(nlevel-1))
	{
	if(i==(nlevel-1))
		{
		rando=paste0(rando,clustername[i])		
		}
	else if(i<(nlevel-1))
		rando=paste0(rando,clustername[i],"/")
	}
rando=formula(rando)
	}

if(p==1)	{if(stand=="location")	{
			if(method=="identity")	{
				yc=y-mean(y)
						}
			else if(method=="sign")	{
				yc=sign(y-median(y))
						}
			else if(method=="rank")	{
				ray=rank(y)	
				yc=ray-mean(ray)	
						}
					}
		else if(stand=="reg")	{
			if(method=="identity")	{
				yc=lm(fixed,data)$residuals
						}
			if(method=="mixed")	{
				yc=residuals(lme(fixed,data,rando),level=0)
						}
			else if(method=="sign")	{
				yc=sign(rq(fixed,data,tau=.5)$residuals)
						}
			else if(method=="rank")	{
				resr=rank(rfit(fixed,data)$residuals)
				yc=resr-mean(resr)
						}
					}
		} # if(p==1)

else if(p>1)	{
			if(stand=="location")	
			{
				estim=mv.1sample.est(y,score=method,stand="inner",maxiter=10000)
				u=estim$location
				yres=y-matrix(u,nrow=N,ncol=p,byrow=T)
			}
			else if(stand=="reg")
			{
				if(method=="identity")	{
					yres=lm(fixed,data=data)$residuals	
								}
				else if(method=="sign" | method=="rank")
								{
					yres=mv.l1lm(fixed,scores=method,stand="inner",maxiter=10000,data=data)$residuals							
								}			
			}
			if(affequiv==TRUE)
				{	
				if(method=="identity" | method=="sign")
					{	
					sha=mv.shape.est(yres,score=method,location=rep(0,p),estimate="inner")
					}
				else if(method=="rank")
					{
					sha=mv.shape.est(yres,score=method,estimate="inner")
					sha=sha*p/sum(diag(sha))
					}
				ev=eigen(sha)
				svi=ev$vectors %*% sqrt(solve(diag(ev$values))) %*% t(ev$vectors)
				sv=ev$vectors %*% sqrt(diag(ev$values)) %*% t(ev$vectors)
					}

			if(affequiv==FALSE)
				{
				sv=diag(p)
				svi=diag(p)	
				}
			
			if(method=="identity")
				{
				yc=yres %*% svi 
				}
			else if(method=="sign")
				{
				yc=spatial.signs(yres %*% svi,center=FALSE,shape=FALSE)
				}
			else if(method=="rank")
				{
				yc=spatial.rank(yres %*% svi,center=FALSE,shape=FALSE)
				}
		} # else if(p>1)

yc=as.matrix(yc)

list(yc,sv)

}

