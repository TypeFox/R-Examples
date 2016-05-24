"denoisehetero2"<-function(x,f,pred=1,neigh=1,int=TRUE,clo=FALSE,keep=2,rule="median",verbose=TRUE,index=1:length(x),returnall=FALSE){

nonorcoeff<-NULL
nordetlist<-NULL
nortclist<-NULL

al<-list()
temp<-NULL
sdvec<-NULL
sdvec1<-NULL
sdvec2<-NULL
newcoeff<-NULL
newcoeff1<-NULL
newcoeff2<-NULL
fhat<-NULL
fhat1<-NULL
fhat2<-NULL

n<-length(x)
out<-fwtnp2(x,f,pred,neigh,int,clo,keep)
x<-out$x
lca<-out$lca

po<-out$pointsin
nonorcoeff<-out$coeff
lr<-out$lengthsremove    #vector deciding how to divide up coefficients into artificial levels
rem<-lca[,1] 		     #used to convert output to original lr,rem)

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
	newcoeff1[al[[i]]]<-nortclist*(sd[al[[i]]])

	nordetlist<-nonorcoeff[al[[i]]]/(sd2[al[[i]]])
	nortclist<-ebayesthresh(nordetlist,prior="cauchy",a=NA,sdev=1,threshrule=rule)
	newcoeff2[al[[i]]]<-nortclist*(sd[al[[i]]])
}

newcoeff[po]<-out$coeff[po]
newcoeff1[po]<-out$coeff[po]
newcoeff2[po]<-out$coeff[po]

neighbrs<-list()
for(kk in 1:length(rem)){
        neighbrs[[kk]]<-.C("nbrsfromlca",as.double(t(lca)),as.integer(ncol(lca)),as.integer(kk),
        n=as.integer(rep(0,times=lca[kk,2])),PACKAGE="adlift")$n
}
                                                                                
adds <- findadds(rem, neighbrs, po, index)
add <- max(adds)
if (verbose) {
	cat("doing ", add, " steps...\n", sep = "")
}

fhat <- invtnp2(x, newcoeff, out$lengths, lr, po,lca, add)
fhat1 <- invtnp2(x, newcoeff1, out$lengths, lr, po,lca, add)
fhat2 <- invtnp2(x, newcoeff2, out$lengths, lr, po,lca, add)

if(returnall){                                                                                                                    
	return(list(fhat=fhat,fhat1=fhat1,fhat2=fhat2,al=al,sd=sd))
}
else{
	return(fhat$coeff)
}

}
