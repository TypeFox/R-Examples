"basisfns" <-
function(x,f,pred,neigh,int,clo,keep,plot.f=FALSE,plot.bas=FALSE,separate=FALSE){

#produces plots of PRIMAL wavelet (psi) basis functions from a specified transform

#since multiplying Winv with a delta vector gives the columns of Winv, 
#the basis func. vectors are just the columns of winv.

out<-fwtnp(x,f,LocalPred=pred,neighbours=neigh,intercept=int,closest=clo,nkeep=keep,do.W=TRUE,varonly=FALSE)
w<-out$W

basmat<-matrix(0,length(f),length(f))
maxv<-NULL
maxb<-NULL
minb<-NULL
schemehist<-out$schemehist
interhist<-out$interhist
pointsin<-out$pointsin			#so that we know which basis functions are the scaling functions
removelist<-out$removelist

winv<-solve(w)

schhist<-matrix(0,1,length(f))
inthist<-matrix(0,1,length(f))

#reorders schemehist, interhist to be in order of x

minx<-min(x)
maxx<-max(x)
minf<-min(f)
maxf<-max(f)
rx<-maxx-minx
rf<-maxf-minf

if (plot.f==TRUE){
	cls<-class(getOption("device"))	#bug fix for R-2.8.0beta +  
	if(cls=="function"){
		getOption("device")()
	}
	else{
		get(getOption("device"))()
	}

	plot(x,seq(minf,maxf+(rf*.5),length=length(x)),type="n",xlab="",ylab="")
	lines(sort(x),f[order(x)],type="l")
	if(!is.null(schemehist)){
		leg<-c("LinearPred","QuadPred","CubicPred","Basis function")
		legend(.65*rx,maxf+(rf*.5),leg,fill=c(2,3,4,5))
	}
	orig<-dev.cur()
}

for (k in 1:length(f)){
	colour<-1
	basmat[k,]<-winv[,k]   #each ROW of basmat is basis fun  

	if(is.null(schemehist)){
		schhist<-NULL
	}
	else{
		q<-which(removelist==k)  #finds out the scheme of removelist in order of x
		if (length(q)!=0){
			schhist[k]<-schemehist[q]
			inthist[k]<-interhist[q]
		}

		if(schhist[k]=="Linear"){
			colour<-2
		}
		if(schhist[k]=="Quad"){
			colour<-3
		}
		if(schhist[k]=="Cubic"){
			colour<-4
		}
		if(schhist[k]=="0"){
			colour<-5
		}

	}  

	minb[k]<-min(basmat[k,])
	maxb[k]<-max(basmat[k,])
	maxv[k]<-order(abs(basmat[k,]))[length(f)]

	if (plot.f==TRUE){
		dev.set(orig)
		if(!is.null(schhist)){
			lines(rep(x[k],times=10),seq(minf,maxf,length=10),type="l",col=colour)   #plots line at datapoint (x axis)
		}
	}
}  

# Note: basmat values will be in the "order" of x[1:length(f)] , so still need to sort

if (plot.f==TRUE){
	dev.set(orig)
	lines(sort(x),f[order(x)])
}

c<-setdiff(plot.bas,0)

if (length(c)!=0){
	m<-min(minb[c])
	M<-max(maxb[c])
	if (any(is.na(match(plot.bas,1:length(f))))){
		stop("can't plot requested basis functions")
	}
	else{
		if (separate==FALSE){
			cls<-class(getOption("device"))	#bug fix for R-2.8.0beta +  
			if(cls=="function"){
				getOption("device")()
			}
			else{
				get(getOption("device"))()
			}

			plot(x,seq(m,M+((M-m)*.5),length=length(x)),type="n",xlab="x",ylab="basis function")
			newdev<-dev.cur()
		}

		for (i in plot.bas){

			if(is.null(schhist[i])){
				colour<-i
			}
			else{
				if(schhist[i]=="Linear"){
					colour<-2
				}
				if(schhist[i]=="Quad"){
					colour<-3
				}
				if(schhist[i]=="Cubic"){
					colour<-4
				}
				if(schhist[i]=="0"){
					colour<-5
				}
			}

			if (separate==TRUE){
				cls<-class(getOption("device"))	#bug fix for R-2.8.0beta +  
				if(cls=="function"){
					getOption("device")()
				}
				else{
					get(getOption("device"))()
				}

				plot(x,seq(minb[i],maxb[i]+((maxb[i]-minb[i])*.5),length=length(x)),type="n",xlab="x",ylab="basis function")
				lines(sort(x),basmat[i,][order(x)],type="l",col=colour)
				if(!is.null(schhist[i])){
					leg<-c("LinearPred","QuadPred","CubicPred","Basis function")
					legend(.65*rx,maxb[i]+((maxb[i]-minb[i])*.5),leg,fill=c(2,3,4,5))
				}
			} 
			else{
				dev.set(newdev)
				lines(sort(x),basmat[i,][order(x)],type="l",col=colour)
				if(!is.null(schhist[i])){
					leg<-c("LinearPred","QuadPred","CubicPred","Basis function")
					legend(.65*rx,M+((M-m)*.5),leg,fill=c(2,3,4,5))
				}
			}
		}	#end for
	}	#end else
}	#end if

return(list(out=out,maxv=maxv,schhist=schhist,inthist=inthist,basmat=basmat))

}

