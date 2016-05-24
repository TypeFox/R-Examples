"denoisehetero" <-function(x,f,pred=LinearPred,neigh=1,int=TRUE,clo=FALSE,keep=2,rule="median",returnall=FALSE){

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

n<-length(x)

out<-fwtnp(x,f,LocalPred=pred,neighbours=neigh,intercept=int,closest=clo,nkeep=keep,do.W=FALSE,varonly=FALSE)

po<-out$pointsin
nonorcoeff<-out$coeff
lr<-out$lengthsremove    #vector deciding how to divide up coefficients into artificial levels
rem<-out$removelist      #used to convert output to original lr,rem)

al<-artlev(lr,rem)      #the list of indices of removelist separated into levels
levno<-length(al)


detail<-y<-matrix(0,1,n-keep)

y<-x[setdiff(1:n,po)]
detail<-nonorcoeff[setdiff(1:n,po)]
h<-heterovar(y,detail,al)

sdvec<-h$varvec
sdvec1<-h$varvec1
sdvec2<-h$varvec2

sd<-sd1<-sd2<-NULL

sd[setdiff(1:n,po)]<-sdvec
sd1[setdiff(1:n,po)]<-sdvec1
sd2[setdiff(1:n,po)]<-sdvec2

sd[po]<-NA
sd1[po]<-NA
sd2[po]<-NA

for (i in 1:levno){
	nordetlist<-nonorcoeff[al[[i]]]/(sd[al[[i]]])
	nortclist<-ebayesthresh(nordetlist,prior="cauchy",a=NA,sdev=1,threshrule=rule)
	newcoeff[al[[i]]]<-nortclist*(sd[al[[i]]])

	nordetlist<-nonorcoeff[al[[i]]]/(sd1[al[[i]]])
	nortclist<-ebayesthresh(nordetlist,prior="cauchy",a=NA,sdev=1,threshrule=rule)
	newcoeff1[al[[i]]]<-nortclist*(sd1[al[[i]]])

	nordetlist<-nonorcoeff[al[[i]]]/(sd2[al[[i]]])
	nortclist<-ebayesthresh(nordetlist,prior="cauchy",a=NA,sdev=1,threshrule=rule)
	newcoeff2[al[[i]]]<-nortclist*(sd2[al[[i]]])
}

newcoeff[po]<-out$coeff[po]
newcoeff1[po]<-out$coeff[po]
newcoeff2[po]<-out$coeff[po]

fhat<-invtnp(x,newcoeff,out$lengths,lr,po,rem,out$neighbrs,out$schemehist,out$interhist,n-keep,int,neigh,clo,pred)
fhat1<-invtnp(x,newcoeff1,out$lengths,lr,po,rem,out$neighbrs,out$schemehist,out$interhist,n-keep,int,neigh,clo,pred)
fhat2<-invtnp(x,newcoeff2,out$lengths,lr,po,rem,out$neighbrs,out$schemehist,out$interhist,n-keep,int,neigh,clo,pred)

if(returnall){
	return(list(fhat=fhat,fhat1=fhat1,fhat2=fhat2,w=out$W,al=al,sd=sd))
}
else{
	return(fhat$coeff)
}

}
