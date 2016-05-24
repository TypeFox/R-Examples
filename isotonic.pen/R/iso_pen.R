iso_pen <-
function(y,xmat,wt=1,pen=TRUE,default=TRUE,lambda=0,nsim=0,alpha=.05){
	if(length(xmat)==length(y)){xmat=matrix(xmat,ncol=1)}
	k=dim(xmat)[2]
	sm=1e-4
	if(length(wt)==1){
		usewt=FALSE
	}else{
		usewt=TRUE
		if(length(wt)!=length(y)){print("ERROR: weight vector not correct length")}
		if(!all(wt>=0)){print("ERROR: must have non-negative weights")}
		if(!all(wt>0)){keep=wt>0;y=y[keep];wt=wt[keep];xmat=xmat[keep,]}
		wt=n*wt/sum(wt)
	}
	n=length(y)
	if(dim(xmat)[1]!=n){print("ERROR: number of rows of xmat must be length of y")}
	if(k>3){print("ERROR: xmat can have at most 3 columns")}
	if(rankMatrix(xmat)[1]<k-sm){print("ERROR: xmat is not full rank")}
## clean up xmat
	mm=matrix(nrow=2,ncol=k)
	for(i in 1:k){mm[1,i]=min(xmat[,i]);mm[2,i]=max(xmat[,i])}
	for(i in 1:k){xmat[,i]=(xmat[,i]-min(xmat[,i]))/(max(xmat[,i]-min(xmat[,i])))}
	xmat=round(xmat,10)
	x=matrix(as.numeric(xmat),ncol=k)
	if(k==1){
		if(pen&default){
			bmult=sum((x-mean(x))*(y-mean(y)))/sum((x-mean(x))^2)
			if(bmult>0){
				y=y/bmult
				lambda=n^(1/3)
			}else{
				lambda=n^(1/2)
				bmult=1
			}
		}else if(!pen){
			lambda=0
			bmult=1
		}
		delta=makedelta1(x)
		one=matrix(1:n*0+1,ncol=1)		
		if(usewt){
			yw=y*sqrt(wt);delw=delta%*%diag(sqrt(wt));onew=one*sqrt(wt)
		}else{
			yw=y;delw=delta;onew=one
		}
		if(pen){
			fit0=coneB(yw,delw,onew)			
			sighat=sqrt(sum((yw-fit0$yhat)^2)/(n-1.5*fit0$df))
			lambda=lambda*sighat
		}
		nsm=sum(x==min(x))
		ytil=yw;ytil[x==min(x)]=yw[x==min(x)]+lambda/2/nsm
		nlg=sum(x==max(x))
		ytil[x==max(x)]=yw[x==max(x)]-lambda/2/nlg
		wfit=coneB(ytil,delw,onew)
		fit1=wfit$yhat
		sigest=sqrt(sum((yw-fit1)^2)/(n-1.5*wfit$df))
		if(nsim>0){  ## get confidence interval in the transformed model
			fitmat=matrix(nrow=nsim,ncol=n)
			for(isim in 1:nsim){
				ysim=fit1+rnorm(n)*sigest
				ysim[x==min(x)]=ysim[x==min(x)]+lambda/2/nsm
				ysim[x==max(x)]=ysim[x==max(x)]-lambda/2/nlg
				fitmat[isim,]=coneB(ysim,delw,onew)$yhat
			}
			low=round(nsim*alpha/2);upp=round(nsim*(1-alpha/2))
			upper=1:n;lower=1:n
			for(i in 1:n){st=sort(fitmat[,i]);upper[i]=st[upp];lower[i]=st[low]}
		}	
		if(usewt){
			fit1=fit1/sqrt(wt)
			if(nsim>0){
				upper=upper/sqrt(wt)
				lower=lower/sqrt(wt)
			}
		}	
	}else{  ## multiple isotonic
		if(k==2){
			np=n+2
			if(pen&default){
				m1=lm(y~xmat)
				bmult=as.numeric(m1$coef[2]+m1$coef[3])
				if(bmult>0){y=y/bmult}else{bmult=1}
				lambda=sqrt(n)/2
			}else{bmult=1}
			ox=order(xmat[,1],xmat[,2])
			rk=order(ox)
			if(!usewt){wt=1:n*0+1}
			wp=1:np*0+sm;wp[2:(n+1)]=wt
			smat=matrix(nrow=np,ncol=2)
			smat[2:(n+1),]=xmat[ox,]
			smat[1,]=1:k*0-.01
			smat[np,]=1:k*0+1.01
			yp=1:np;yp[2:(n+1)]=y[ox];yp[1]=min(y);yp[np]=max(y)			
			au=getunique(cbind(smat,yp,wp),2)
			sumat=au$xmat
			yu=au$y
			wu=au$wt
			nu=length(yu)
			amat0=makeamat(sumat[2:(nu-1),],2)$amat
			m=dim(amat0)[1]
			obs=1:m
			n0=nu-2
			minx=1:n0<0
			for(i in 1:n0){if(sum(amat0[,i]!=0)>0){first=min(obs[amat0[,i]!=0]);if(amat0[first,i]==-1){minx[i]=TRUE}}}
			maxx=1:n0<0
			for(i in 1:n0){if(sum(amat0[,i]!=0)>0){last=max(obs[amat0[,i]!=0]);if(amat0[last,i]==1){maxx[i]=TRUE}}}
			if(sum(maxx)==1&sum(minx)==1){  ### unique max & min -- easy!
				m=dim(amat0)[1]
				nu=nu-2
				wu=wu[2:(nu+1)];yu=yu[2:(nu+1)];augroup=au$group[2:(n+1)]
				yw=yu*sqrt(wu);aw=amat0
				for(i in 1:nu){aw[,i]=amat0[,i]/sqrt(wu[i])}
				ans0=coneA(yw,aw)
				if(nu>1.5*ans0$df){  
					sigest=sqrt(sum((yw-ans0$thetahat)^2)/(nu-ans0$df*1.5))
	       	    }else{sigest=sqrt(sum((yw-ans0$thetahat)^2)/(nu-ans0$df+1))}
				if(pen){
					if(default){
	    	    	    lambda=lambda*sigest
	    	    	}
	    	        yw[1]=yw[1]+lambda/2/sqrt(wu[1]);yw[nu]=yw[nu]-lambda/2/sqrt(wu[nu])
	    	        ans0=coneA(yw,aw)
	    	    }
	    	    if(nsim>0){   ### do confidence intervals in the transformed model
	    	    	fitmat=matrix(nrow=nsim,ncol=nu)
	    	    	for(isim in 1:nsim){
	    	    		ysim=ans0$thetahat+sigest*rnorm(nu)
	    	    		if(pen){ysim[1]=ysim[1]+lambda/2/sqrt(wu[1]);ysim[nu]=ysim[nu]-lambda/2/sqrt(wu[nu])}
	    	    		fitmat[isim,]=coneA(ysim,aw)$thetahat
	    	    	}
	    	    	lq=round(nsim*alpha/2);uq=round(nsim*(1-alpha/2))
	    	    	upq=1:nu;lwq=1:nu
	    	    	for(i in 1:nu){vals=sort(fitmat[,i]);upq[i]=vals[uq];lwq[i]=vals[lq]}
	    	    }  ## now, untransform
		 		thetahat=ans0$thetahat/sqrt(wu)
		 		if(nsim>0){lwq=lwq/sqrt(wu);upq=upq/sqrt(wu)}
				agr=1:n   ## expand
				for(i in 1:nu){agr[augroup==unique(augroup)[i]]=i}
				thb=thetahat[agr]
				theta=thb[rk]
				if(nsim>0){
					lb=lwq[agr];lower=lb[rk]
					ub=upq[agr];upper=ub[rk]
				}
				sighat=sigest
			}else{  ## non-unique max & min -- make fake data with small weight
				amat=makeamat(sumat,2)$amat
				m=dim(amat)[1]
				yw=yu*sqrt(wu);aw=amat
				for(i in 1:nu){aw[,i]=amat[,i]/sqrt(wu[i])}
				ans0=coneA(yw,aw)
				if(pen){
					if(default){
						if(nu>1.5*ans0$df){  
							sigest=sqrt(sum((yw[2:(nu-1)]-ans0$thetahat[2:(nu-1)])^2)/(nu-ans0$df*1.5))
	    	    	    }else{sigest=sqrt(sum((yw[2:(nu-1)]-ans0$thetahat[2:(nu-1)])^2)/(nu-ans0$df+1))}
	    	    	    lambda=lambda*sigest
	    	    	}
	    	        yw[1]=yw[1]+lambda/2/sqrt(wu[1]);yw[nu]=yw[nu]-lambda/2/sqrt(wu[nu])
	 				ans=coneA(yw,aw)  ## weighted, penalized
	 				thetahat=ans$thetahat
					if(nu>1.5*ans$df){  
						sighat=sqrt(sum((yw[2:(nu-1)]-ans$thetahat[2:(nu-1)])^2)/(nu-ans$df*1.5))
	    	        }else{sighat=sqrt(sum((yw[2:(nu-1)]-ans$thetahat[2:(nu-1)])^2)/(nu-ans$df+1))}
	 			}else{
	 				thetahat=ans0$thetahat
					if(nu>1.5*ans0$df){  
						sighat=sqrt(sum((yw[2:(nu-1)]-ans0$thetahat[2:(nu-1)])^2)/(nu-ans0$df*1.5))
	    	        }else{sighat=sqrt(sum((yw[2:(nu-1)]-ans0$thetahat[2:(nu-1)])^2)/(nu-ans0$df+1))}
	    	   }
		 		if(nsim>0){  ## do confidence bounds in transformed model
					fitmat=matrix(nrow=nsim,ncol=nu)
					for(isim in 1:nsim){
						ysim=thetahat+rnorm(nu)*sighat
						if(pen){
							ysim[1]=ysim[1]+lambda/2/sqrt(wu[1])
							ysim[nu]=ysim[nu]-lambda/2/sqrt(wu[nu])
						}
						fitmat[isim,]=coneA(ysim,aw)$thetahat
					}
	    	    	lq=round(nsim*alpha/2);uq=round(nsim*(1-alpha/2))
	    	    	upq=1:nu;lwq=1:nu
	    	    	for(i in 1:nu){vals=sort(fitmat[,i]);upq[i]=vals[uq];lwq[i]=vals[lq]}
	 			}   ## now untransform
		 		thetahat=thetahat/sqrt(wu)
	 			if(nsim>0){lwq=lwq/sqrt(wu);upq=upq/sqrt(wu)}
				agr=1:np  ## expand
				for(i in 1:nu){agr[au$group==unique(au$group)[i]]=i}
				thb=thetahat[agr]
				thb=thb[2:(n+1)]
				theta=thb[rk]
				if(nsim>0){
					upp=upq[agr];low=lwq[agr]
					uppr=upp[2:(n+1)];lowr=low[2:(n+1)]
					upper=uppr[rk];lower=lowr[rk]
				}
			}
			xg1=0:20/20;xg2=0:20/20;xgmat=matrix(nrow=21,ncol=21)
			for(i1 in 1:21){
				for(i2 in 1:21){
					a1=FALSE;a2=FALSE
					if(sum(xmat[,1]>=xg1[i1]&xmat[,2]>=xg2[i2])>0){
						aup=min(theta[xmat[,1]>=xg1[i1]&xmat[,2]>=xg2[i2]])
						a1=TRUE
					}
					if(sum(xmat[,1]<=xg1[i1]&xmat[,2]<=xg2[i2])>0){
						alw=max(theta[xmat[,1]<=xg1[i1]&xmat[,2]<=xg2[i2]])
						a2=TRUE
					}
					if(a1&a2){
						xgmat[i1,i2]=(aup+alw)/2
					}else if(a1&!a2){
						xgmat[i1,i2]=aup
					}else if(a2&!a1){
						xgmat[i1,i2]=alw
					}else{print("ERROR in 2d plot")}
				}
			}			
		}else{print("ERROR -- too many columns in xmat")}
	}
	ans=new.env()
	if(k==1){
		ans$fit=fit1*bmult
		ans$lambda=lambda*bmult
		ans$sighat=sighat*bmult
		if(nsim>0){
			ans$upper=upper*bmult
			ans$lower=lower*bmult
		}
	}else if(k==2){
		ans$fit=theta*bmult
		if(nsim>0){
			ans$upper=upper*bmult
			ans$lower=lower*bmult
		}
		ans$xg1=xg1*(mm[2,1]-mm[1,1])+mm[1,1]
		ans$xg2=xg2*(mm[2,2]-mm[1,2])+mm[1,2]
		ans$xgmat=xgmat*bmult
	}else{
		ans$fit=theta*bmult
		if(nsim>0){
			ans$upper=upper*bmult
			ans$lower=lower*bmult
		}
	}
	ans
}
