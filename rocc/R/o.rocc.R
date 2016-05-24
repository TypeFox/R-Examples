o.rocc <-
function (g,out,xgenes=200) 
{

  if(!is.factor(out))
    stop("out should be a factor (with labels 0 and 1 and exactly in this order)")
  if(levels(out)[1]!="0")
    stop("levels of out have to be 0 and 1 and exactly in this order")  
  if(levels(out)[2]!="1")
    stop("levels of out have to be 0 and 1 and exactly in this order")  
  if(!is.matrix(g))
    warning("g should be a matrix (with genes as rows and samples as columns)")  
  if(is.null(colnames(g)))
    warning("Colnames for g with sample names are missing")
  if(is.null(rownames(g)))
    stop("Rownames for g with gene names are missing")    
  if(0%in%xgenes)
    stop("xgenes must not contain 0")
  if(!is.vector(xgenes))
    stop("xgenes must be a vector")   

require(ROCR)





#############   CROSSVALIDATION LOOP

g<-as.matrix(g)
out<-factor(out,levels=c(0,1))

confusion<-matrix(nrow=17,ncol=length(xgenes))
confusion<-as.data.frame(confusion)
rownames(confusion)[1:17]<-c("Accuracy","AccuracyLower95CI","AccuracyUpper95CI","AccurayNull", "AccuracyNull_pValue","AccuracyRandomAssignment","AccuracyRandomAssignment_pValue","Sensitivity","Specificity","PPV",
"NPV","Prevalence","predicted1-true1","predicted0-true1","predicted1-true0","predicted0-true0","Balanced Accuracy")

concordance<-matrix(nrow=length(out),ncol=length(xgenes))
concordance<-as.data.frame(concordance)
rownames(concordance)<-colnames(g)


for (size in 1:length(xgenes)){


colnames(confusion)[size]<-xgenes[size]
colnames(concordance)[size]<-xgenes[size]

pr<-as.numeric(rep(NA,length(colnames(g))))
pr<-factor(pr,level=c(0,1))
names(pr)<-colnames(g)

for(v in 1:length(colnames(g))){

## Prepare Dataset and Outcome for LOOCV
f<-g[,-v]
outf<-out[-v]


####### Feature Selection in each loop

aucv<-vector(mode = "numeric", length = length(rownames(f)))
names(aucv)<-rownames(f)

for (v1 in 1:length(rownames(f))){
x1 = f[v1,outf==1]
n1 = length(x1)
x2 = f[v1,outf==0]
n2 = length(x2)
r = rank(c(x1,x2))
aucv[v1] = (sum(r[1:n1]) - n1*(n1+1)/2) / (n1*n2)
}


auciv<-as.data.frame(aucv)

posneg<-rep("NA",length(rownames(auciv)))
posneg[which(auciv[,1]<0.5)]<-"neg"
posneg[posneg=="NA"]<-"pos"
auciv$posneg<-posneg

allpos<-vector(mode = "numeric", length = length(rownames(auciv)))
allpos[which(auciv$posneg=="pos")]<-auciv$aucv[which(auciv$posneg=="pos")]
allpos[which(auciv$posneg=="neg")]<-1-auciv$aucv[which(auciv$posneg=="neg")]
auciv$allpos<-allpos

aucivsort<-auciv[order(auciv$allpos,decreasing = T),]
aucivsort<-aucivsort[,c(3,2,1)]

genelistf<-rownames(aucivsort)[1:xgenes[size]]


#### combined expression

e<-f[genelistf,]
ifelse(length(genelistf)==1,e<-t(as.matrix(e)),e<-e)  #### To secure for xgenes = 1
minus<-e[which(aucivsort$posneg[1:length(genelistf)]=="neg"),]
plus<- e[which(aucivsort$posneg[1:length(genelistf)]=="pos"),]
minustoplus<-minus*-1
allplus<-rbind(plus,minustoplus)
meanvector<-apply(allplus,2,mean)


###### find best accuracy cuttoff with according value
pred<- prediction(meanvector, outf)
perf.acc<-performance(pred,"acc")
cutoffpoint<-(which.max(perf.acc@y.values[[1]]))-1   ##### shift back -1 because first value is Inf
ord<-as.matrix(unlist(perf.acc@x.values))
ifelse(ord[1,1]==Inf,ord<-ord[-1,1],ord[,1])


ifelse(cutoffpoint==0,cutoffvalue<-ord[1],
       ifelse(cutoffpoint==dim(e)[2],cutoffvalue<-ord[dim(e)[2]],
              cutoffvalue<-(ord[cutoffpoint]+ord[cutoffpoint+1])/2    ))


positiv <- rownames(aucivsort)[which(aucivsort$posneg[1:length(genelistf)]=="pos")]
negativ <- rownames(aucivsort)[which(aucivsort$posneg[1:length(genelistf)]=="neg")]


#####  Classification
plusn<-g[positiv,v]
minusn<-g[negativ,v]
minustoplusn<-minusn*-1
newsample<-mean(c(plusn,minustoplusn))

ifelse(newsample>cutoffvalue,pr[v]<-1,pr[v]<-0)

print(v)
}

concordance[,size]<-pr
tab<-table(pr,out)

confusion[1,size]<-(tab[1,1]+tab[2,2])/sum(tab)
confusion[2,size]<-binom.test(tab[1,1]+tab[2,2],sum(tab))$conf.int[1]
confusion[3,size]<-binom.test(tab[1,1]+tab[2,2],sum(tab))$conf.int[2]
anu<-max(table(out)/sum(table(out)))
confusion[4,size]<-anu
confusion[5,size]<-binom.test(tab[1,1]+tab[2,2],sum(tab),p=anu,alternative="greater")$p.value
ara<-anu^2+((1-anu)^2)
confusion[6,size]<-ara
confusion[7,size]<-binom.test((tab[1,1]+tab[2,2]),sum(tab),p=ara,alternative="greater")$p.value
confusion[8,size]<-tab[2,2]/(tab[2,2]+tab[1,2])
confusion[9,size]<-tab[1,1]/(tab[1,1]+tab[2,1])
confusion[10,size]<-tab[2,2]/(tab[2,2]+tab[2,1])
confusion[11,size]<-tab[1,1]/(tab[1,1]+tab[1,2])
confusion[12,size]<-(tab[2,2]+tab[1,2])/sum(tab)
confusion[13,size] <- tab[2,2]
confusion[14,size] <- tab[1,2]
confusion[15,size] <- tab[2,1]
confusion[16,size] <- tab[1,1]
confusion[17,size]<-(confusion[8,size]+confusion[9,size])/2 

print(paste("Size=",xgenes[size],"done"))
print(date())
}

concordance_out<-cbind(concordance,out)

JVAL<- list(confusion = confusion, concordance = concordance_out, method = "ROC.based.predictor")
class(JVAL) <- "orocc"
return(JVAL)

}

