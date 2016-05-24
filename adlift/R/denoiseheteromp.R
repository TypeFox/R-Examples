"denoiseheteromp" <-
function(x,f,pred,neigh,int,clo,keep,rule="median",mpdet="ave",returnall=FALSE){

nonorcoeff<-NULL
nordetlist<-NULL
nortclist<-NULL
al<-list()
temp<-list()
sdvec<-NULL
sdvec1<-NULL
sdvec2<-NULL
newcoeff<-NULL
newcoeff1<-NULL
newcoeff2<-NULL
fhat<-list()
fhat1<-list()
fhat2<-list()
newcoefflist<-list()
newcoefflist1<-list()
newcoefflist2<-list()

out<-fwtnpmp(input=x,f=f,LocalPred=pred,neighbours=neigh,intercept=int,closest=clo,nkeep=keep,mpdet=mpdet)

xnew<-adjustx(x,f,"mean")$sepx


nonorcoeff<-out$coeff
lr<-out$lengthsremove    #vector deciding how to divide up coefficients into artificial levels
rem<-out$removelist      #used to convert output to original lr,rem)

al<-artlev(lr,rem)      #the list of indices of removelist separated into levels
levno<-length(al)

y<-matrix(0,1,(length(xnew)-keep))

detail<-matrix(0,1,(length(xnew)-keep))

y<-xnew[setdiff((1:length(xnew)),out$pointsin)]
detail<-nonorcoeff[setdiff((1:length(xnew)),out$pointsin)]
h<-heterovar(y,detail,al)


sdvec<-h$varvec
sdvec1<-h$varvec1
sdvec2<-h$varvec2

sd<-NULL
for (i in 1:length(xnew)){
	if (i<min(out$pointsin)){
		sd[i]<-sdvec[i]
	}
	if (i==min(out$pointsin)){
		sd[i]<-NA
	}
	if ((i>min(out$pointsin)) & (i<max(out$pointsin))){
		sd[i]<-sdvec[i-1]
	}
	if(i==max(out$pointsin)){
		sd[i]<-NA
	}
	if(i>max(out$pointsin)){
		sd[i]<-sdvec[i-2]
	}
}

sd1<-NULL
for (i in 1:length(xnew)){
	if (i<min(out$pointsin)){
		sd1[i]<-sdvec1[i]
	}
	if (i==min(out$pointsin)){
		sd1[i]<-NA
	}
	if ((i>min(out$pointsin)) & (i<max(out$pointsin))){
		sd1[i]<-sdvec1[i-1]
	}
	if(i==max(out$pointsin)){
		sd1[i]<-NA
	}
	if(i>max(out$pointsin)){
		sd1[i]<-sdvec1[i-2]
	}
}

sd2<-NULL
for (i in 1:length(xnew)){
	if (i<min(out$pointsin)){
		sd2[i]<-sdvec2[i]
	}
	if (i==min(out$pointsin)){
		sd2[i]<-NA
	}
	if ((i>min(out$pointsin)) & (i<max(out$pointsin))){
		sd2[i]<-sdvec2[i-1]
	}
	if(i==max(out$pointsin)){
		sd2[i]<-NA
	}
	if(i>max(out$pointsin)){
		sd2[i]<-sdvec2[i-2]
	}
}


for (i in 1:levno){
	nordetlist<-nonorcoeff[al[[i]]]/(sd[al[[i]]])
	nortclist<-ebayesthresh(nordetlist,prior="cauchy",a=NA,sdev=1,threshrule=rule)
	newcoeff[al[[i]]]<-nortclist*(sd[al[[i]]])

	nordetlist<-nonorcoeff[al[[i]]]/(sd1[al[[i]]])
	nortclist<-ebayesthresh(nordetlist,prior="cauchy",a=NA,sdev=1,threshrule=rule)
	newcoeff1[al[[i]]]<-nortclist*(sd[al[[i]]])

	nordetlist<-nonorcoeff[al[[i]]]/(sd2[al[[i]]])
	nortclist<-ebayesthresh(nordetlist,prior="cauchy",a=NA,sdev=1,threshrule=rule)
	newcoeff2[al[[i]]]<-nortclist*(sd[al[[i]]])
}

newcoeff[out$pointsin]<-out$coeff[out$pointsin]
newcoeff1[out$pointsin]<-out$coeff[out$pointsin]
newcoeff2[out$pointsin]<-out$coeff[out$pointsin]

for (i in 1:length(newcoeff)){
	newcoefflist[[i]]<-newcoeff[i]
	newcoefflist1[[i]]<-newcoeff1[i]
	newcoefflist2[[i]]<-newcoeff2[i]
}

fhat<-invtnpmp(x,newcoefflist,newcoeff,out$lengths,lr,out$pointsin,rem,out$neighbrs,out$newneighbrs,out$schemehist,out$interhist,length(x)-keep,int,neigh,clo,pred,mpdet)
fhat1<-invtnpmp(x,newcoefflist,newcoeff1,out$lengths,lr,out$pointsin,rem,out$neighbrs,out$newneighbrs,out$schemehist,out$interhist,length(x)-keep,int,neigh,clo,pred,mpdet)
fhat2<-invtnpmp(x,newcoefflist,newcoeff2,out$lengths,lr,out$pointsin,rem,out$neighbrs,out$newneighbrs,out$schemehist,out$interhist,length(x)-keep,int,neigh,clo,pred,mpdet)

if(returnall){
	return(list(fhat=fhat,fhat1=fhat1,fhat2=fhat2,al=al,sd=sd))
}
else{
	return(fhat$coeff)
}


}







