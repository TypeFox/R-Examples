tr.rocc <-
function (g,out,xgenes=200) 
{

  if(length(xgenes)!=1)
    stop("For tr.rocc function the length of xgenes has to be 1")
  if(!is.vector(xgenes))
    stop("xgenes must be a vector")  
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

require(ROCR)


####### Feature Selection according to highest AUC
g<-as.matrix(g)
out<-factor(out,levels=c(0,1))
f<-g
outf<-out
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

genelistf<-rownames(aucivsort)[1:xgenes]
                            

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
names(cutoffvalue)<-NULL

positiv <- rownames(aucivsort)[which(aucivsort$posneg[1:length(genelistf)]=="pos")]
negativ <- rownames(aucivsort)[which(aucivsort$posneg[1:length(genelistf)]=="neg")]


DUVAL<- list(AUCs = aucivsort[genelistf,], genes = genelistf, positiv = positiv, negativ = negativ, metagene.expression = meanvector, metagene.expression.ranked = ord, cutoffvalue = cutoffvalue, method = "ROC.based.predictor")
class(DUVAL) <- "trocc"
return(DUVAL)

}

