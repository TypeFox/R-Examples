superpc.listfeatures<- function(data, train.obj, fit.red,  fitred.cv=NULL,
num.features=NULL, component.number=1){

ii=component.number
total.num=sum(abs(fit.red$import[,ii])>0)

if(is.null(num.features)){ num.features=total.num}

if(num.features< 1 | num.features > total.num){

    stop("Error: num.features   argument out of range")

}

featurenames.short<- substring(data$featurenames,1,40)



oo=rank(abs(fit.red$import[,ii]))> nrow(data$x)-num.features

res<-cbind(round(fit.red$import[oo,ii],3), round(train.obj$feature.scores[oo],3),
#round(fit.red$wt[oo,ii],3),
 featurenames.short[oo])

collabs=c("Importance-score", "Raw-score" , "Name")

if(!is.null(data$featureid)){
  res=cbind(res, data$featureid[oo])
 collabs=c(collabs, "ID")
}

if(!is.null(fitred.cv)){
nfold=ncol(fitred.cv$import.cv)

ind=matrix(F,nrow=nrow(data$x),ncol=nfold)
ranks=NULL
for( j in 1:nfold){
          r <- fitred.cv$import.cv[,j,component.number]
        ranks=cbind(ranks,rank(-abs(r)))

  junk=  fitred.cv$import.cv[,j, component.number]!=0
        ind[junk,j]=T
}

av.rank=apply(ranks,1,median)
av.rank=round(av.rank[oo],2)
prop=apply(ind[oo,,drop=F],1,sum)/nfold
res=cbind(res,av.rank,prop)
collabs=c(collabs,"median-rank-in-CV","prop-selected-in-CV")
}


o<-order(-abs(fit.red$import[oo,ii]))
res<-res[o,]
dimnames(res)<-list(NULL,collabs)

return(res)
}
