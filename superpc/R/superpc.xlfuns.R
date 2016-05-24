superpc.xl.get.threshold.range  <- function(train.obj) {


  cur.tt <- train.obj$feature.scores
p<-length(cur.tt)

min.features=5
max.features=p


  lower <- quantile(abs(cur.tt), 1 - (max.features/p))
  upper <- quantile(abs(cur.tt), 1 - (min.features/p))

  return(c(lower,upper))
}


superpc.xl.listgenes.compute<- function(data, train.obj, fit.red, fitred.cv=NULL,
num.genes=NULL, component.number=1){

ii=component.number
total.num=sum(abs(fit.red$import[,ii])>0)
if(is.null(num.genes)){num.genes=total.num}

if(num.genes< 1 | num.genes > total.num){
    stop("Error: num.genes   argument out of range")

}

genenames.short<- substring(data$geneid,1,40)


oo=rank(abs(fit.red$import[,ii]))> nrow(data$x)-num.genes

res<-cbind(data$genenames[oo], genenames.short[oo], round(fit.red$import[oo,ii],3), round(train.obj$feature.scores[oo],3))

collabs=c("GeneID", "Genename","Importance-score", "Raw-score" )


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

av.rank=apply(ranks,1,mean)
av.rank=round(av.rank[oo],2)
prop=apply(ind[oo,,drop=F],1,sum)/nfold
res=cbind(res,av.rank,prop)
collabs=c(collabs,"av-rank-in-CV","prop-selected-in-CV")
}


o<-order(-abs(fit.red$import[oo,ii]))
res<-res[o,]
dimnames(res)<-list(NULL,collabs)

if(is.null(fitred.cv)){foo=list(gene.ids=res[,1],gene.names=res[,2],gene.scores=res[,3:4],
    gene.headings=c("ID", "Name"," Importance score", "Raw score"))
}
if(!is.null(fitred.cv)){foo=list(gene.ids=res[,1],gene.names=res[,2],gene.scores=res[,3:6],
    gene.headings=c("ID", "Name"," Importance score", "Raw score","av-rank-in-CV", "prop-selected-in-CV"))
}

return(foo)
}



superpc.xl.fit.to.clin<- function(fit, data.test,score, pamr.xl.test.sample.labels, pamr.xl.clindata, pamr.xl.clinsamplabels, pamr.xl.clinprednames, pamr.xl.clinpredtype ){

# strip off extra first element, row, col that rServer seems to include
# put everything in same order as test expression data

errorflag=FALSE
temp= pamr.xl.clindata[-1,-1,drop=F]
clinsamplabels=pamr.xl.clinsamplabels[-1]
clinpredtype=pamr.xl.clinpredtype[-1]
clinprednames=pamr.xl.clinprednames[-1]

o=match(pamr.xl.test.sample.labels, clinsamplabels)
if(sum(is.na(o))>0){ errorflag=TRUE
 return(list(errorflag=errorflag))
}
clindata=temp[,o,drop=F]
npreds=nrow(clindata)
clinlist=vector("list",npreds)
for(i in 1:npreds){
  if(clinpredtype[[i]]=="continuous") {clinlist[[i]]=clindata[i,]}
  if(clinpredtype[[i]]=="discrete") {clinlist[[i]]=as.factor(clindata[i,])}
}

names(clinlist)= clinprednames

foo=superpc.fit.to.outcome(fit, data.test,score, competing.predictors=clinlist, print=FALSE)

return(list(results=foo$coeftable, teststat.table=foo$teststat.table,  errorflag=errorflag))
}


superpc.xl.decorrelate<- function(data, pamr.xl.train.sample.labels, pamr.xl.clindata, pamr.xl.clinsamplabels, pamr.xl.clinprednames, pamr.xl.clinpredtype ){

# strip off extra first element, row, col that rServer seems to include
# put everything in same order as training expression data
errorflag=FALSE

temp= pamr.xl.clindata[-1,-1]
clinsamplabels=pamr.xl.clinsamplabels[-1]
clinpredtype=pamr.xl.clinpredtype[-1]
clinprednames=pamr.xl.clinprednames[-1]

o=match(pamr.xl.train.sample.labels, clinsamplabels)
if(sum(is.na(o))>0){ errorflag=TRUE
 return(list(errorflag=errorflag))
}

clindata=temp[,o]
npreds=nrow(clindata)
clinlist=vector("list",npreds)
for(i in 1:npreds){
  if(clinpredtype[[i]]=="continuous") {clinlist[[i]]=clindata[i,]}
  if(clinpredtype[[i]]=="discrete") {clinlist[[i]]=as.factor(clindata[i,])}
}

names(clinlist)= clinprednames

foo=superpc.decorrelate(data$x, competing.predictors=clinlist)

return(list(results=foo, errorflag=errorflag))

}

superpc.xl.decorrelate.test<- function(object.decorr, xtest, pamr.xl.train.sample.labels, pamr.xl.clindata, pamr.xl.clinsamplabels, pamr.xl.clinprednames, pamr.xl.clinpredtype ){

# strip off extra first element, row, col that rServer seems to include
# put everything in same order as training expression data

errorflag=FALSE
temp= pamr.xl.clindata[-1,-1]
clinsamplabels=pamr.xl.clinsamplabels[-1]
clinpredtype=pamr.xl.clinpredtype[-1]
clinprednames=pamr.xl.clinprednames[-1]

o=match(pamr.xl.train.sample.labels, clinsamplabels)
if(sum(is.na(o))>0){ errorflag=TRUE
 return(list(errorflag=errorflag))
}

clindata=temp[,o]
npreds=nrow(clindata)
clinlist=vector("list",npreds)
for(i in 1:npreds){
  if(clinpredtype[[i]]=="continuous") {clinlist[[i]]=clindata[i,]}
  if(clinpredtype[[i]]=="discrete") {clinlist[[i]]=as.factor(clindata[i,])}
}

names(clinlist)= clinprednames

xtest<- xtest-t(predict(object.decorr, clinlist))

return(list(results=xtest, errorflag=errorflag))

}



superpc.xl.rainbowplot=function(data, pred,  pamr.xl.test.sample.labels, 
    pamr.xl.clindata, pamr.xl.clinsamplabels, pamr.xl.clinprednames, 
    pamr.xl.clinpredtype, call.win.metafile=FALSE){


 temp = pamr.xl.clindata[-1, -1, drop = F]
    clinsamplabels = pamr.xl.clinsamplabels[-1]
    clinpredtype = pamr.xl.clinpredtype[-1]
    clinprednames = pamr.xl.clinprednames[-1]
    o = match(pamr.xl.test.sample.labels, clinsamplabels)
    if (sum(is.na(o)) > 0) {
        errorflag = TRUE
        return(list(errorflag = errorflag))
    }
    clindata = temp[, o, drop = F]
    npreds = nrow(clindata)
    clinlist = vector("list", npreds)
    for (i in 1:npreds) {
        if (clinpredtype[[i]] == "continuous") {
            clinlist[[i]] = clindata[i, ]
        }
        if (clinpredtype[[i]] == "discrete") {
            clinlist[[i]] = as.factor(clindata[i, ])
        }
    }
    names(clinlist) = clinprednames



junk<- superpc.rainbowplot(data, pred, pamr.xl.test.sample.labels,  clinlist, call.win.metafile=call.win.metafile)
return()
}

  
