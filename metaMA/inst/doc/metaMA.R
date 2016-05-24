### R code from vignette source 'metaMA.Rnw'

###################################################
### code chunk number 1: metaMA.Rnw:41-42
###################################################
options(width=60)


###################################################
### code chunk number 2: loadparameters (eval = FALSE)
###################################################
## library(GEOquery)


###################################################
### code chunk number 3: loaddata (eval = FALSE)
###################################################
## data1 = getGEO('GSE9844')
## data2 = getGEO('GSE3524')
## data3 = getGEO('GSE13601')


###################################################
### code chunk number 4: boxplot (eval = FALSE)
###################################################
## par(mfrow=c(2,2))
## boxplot(data.frame(exprs(data1[[1]])),main="data1",outline=FALSE)
## boxplot(data.frame(exprs(data2[[1]])),main="data2",outline=FALSE)
## boxplot(data.frame(exprs(data3[[1]])),main="data3",outline=FALSE)


###################################################
### code chunk number 5: normalizedata (eval = FALSE)
###################################################
## exprs(data3[[1]])<-log2(exprs(data3[[1]]))
## boxplot(data.frame(exprs(data3[[1]])),main="data3",outline=FALSE)


###################################################
### code chunk number 6: classes (eval = FALSE)
###################################################
## c1=as.numeric(pData(data1[[1]])["source_name_ch1"]==
##                 "Oral Tongue Squamous Cell Carcinoma")
## c2=as.numeric(apply(pData(data2[[1]])["description"],
##                     1,toupper)=="SERIES OF 16 TUMORS")
## c3=as.numeric(pData(data3[[1]])["source_name_ch1"]=="Tumor")
## classes=list(c1,c2,c3)


###################################################
### code chunk number 7: classesman (eval = FALSE)
###################################################
## c1=c(1,0,1,1,0,1,0,0,0,1)
## c2=c(1,1,1,0,0,0)
## classes2=list(c1,c2)


###################################################
### code chunk number 8: install bioconductor package (eval = FALSE)
###################################################
## #source("http://bioconductor.org/biocLite.R")
## #biocLite("org.Hs.eg.db")


###################################################
### code chunk number 9: Human Gene database (eval = FALSE)
###################################################
## require("org.Hs.eg.db")
## x <- org.Hs.egUNIGENE
## mapped_genes <- mappedkeys(x)
## link <- as.list(x[mapped_genes])


###################################################
### code chunk number 10: map probes id with unigenes
###################################################
probe2unigene<-function(expset)
{
  #construction of the map probe->unigene
  probes=rownames(exprs(expset))
  gene_id=fData(expset)[probes,"ENTREZ_GENE_ID"]
  unigene=link[gene_id]
  names(unigene)<-probes
  probe_unigene=unigene
}


###################################################
### code chunk number 11: permutation unigene->probe
###################################################
unigene2probe<-function(map)
{
  suppressWarnings(x <- cbind(unlist(map), names(map)))
  unigene_probe=split(x[,2], x[,1]) 
}


###################################################
### code chunk number 12: totalconversion (eval = FALSE)
###################################################
## convert2metaMA<-function(listStudies,mergemeth=mean)
## {
##   if (!(class(listStudies) %in% c("list"))) {
##     stop("listStudies must be a list")
##   } 
##   conv_unigene=lapply(listStudies, 
##                       FUN=function(x) unigene2probe(probe2unigene(x)))
##   id=lapply(conv_unigene,names)
##   
##   inter=Reduce(intersect,id)
##   if(length(inter)<=0){stop("no common genes")}
##   print(paste(length(inter),"genes in common"))
##   esets=lapply(1:length(listStudies),FUN=function(i){
##     l=lapply(conv_unigene[[i]][inter],
##              FUN=function(x) exprs(listStudies[[i]])[x,,drop=TRUE])
##     esetsgr=t(sapply(l,FUN=function(ll) if(is.null(dim(ll))){ll}
##                      else{apply(ll,2,mergemeth)}))
##     esetsgr
##   })
##   return(list(esets=esets,conv.unigene=conv_unigene))
## }
## conv=convert2metaMA(list(data1[[1]],data2[[1]],data3[[1]]))
## esets=conv$esets
## conv_unigene=conv$conv.unigene


###################################################
### code chunk number 13: nb (eval = FALSE)
###################################################
## nb=length(rownames(esets[[1]]))


###################################################
### code chunk number 14: metama (eval = FALSE)
###################################################
## library(metaMA)
## res=pvalcombination(esets=esets,classes=classes)
## length(res$Meta)
## Hs.Meta=rownames(esets[[1]])[res$Meta]


###################################################
### code chunk number 15: IDDIRR (eval = FALSE)
###################################################
## x=IDDIRR(res$Meta,res$AllIndStudies)


###################################################
### code chunk number 16: venndiagram (eval = FALSE)
###################################################
## library(VennDiagram)
## venn.plot<-venn.diagram(x = list(study1=res$study1,
##                                  study2=res$study2,
##                                  study3=res$study3,
##                                  meta=res$Meta),
##              filename = NULL, col = "black",
##              fill = c("blue", "red", "purple","green"),
##              margin=0.05, alpha = 0.6)
## jpeg("venn_jpeg.jpg")
## grid.draw(venn.plot)
## dev.off()


###################################################
### code chunk number 17: annotationdb (eval = FALSE)
###################################################
## getanndb<-function(expset)
## {
##   gpl_name=annotation(expset)
##   gpl=getGEO(gpl_name)
##   title=Meta(gpl)$title
##   db=paste(gsub("[-|_]","",
##                 tolower(strsplit(title,"[[]|[]]")
##                         [[1]][[2]])),"db",sep=".")
##   return(db)
## }
## db1=getanndb(data1[[1]])
## db1
## db2=getanndb(data2[[1]])
## db2
## db3=getanndb(data3[[1]])
## db3


###################################################
### code chunk number 18: getannotationpackages (eval = FALSE)
###################################################
## #source("http://bioconductor.org/biocLite.R")
## #biocLite(db1)
## #biocLite(db2)
## #biocLite(db3)
## library(db1,character.only=TRUE)
## library(db2,character.only=TRUE)
## library(db3,character.only=TRUE)


###################################################
### code chunk number 19: getannotation (eval = FALSE)
###################################################
## origId.Meta=lapply(conv_unigene, 
##                    FUN=function(vec) as.vector(unlist(vec[Hs.Meta])))
## library(annaffy)
## annlist=lapply(1:length(origId.Meta), 
##                FUN=function(i) aafTableAnn(origId.Meta[[i]],chip=get(paste0("db",i)),colnames=aaf.handler()))  
## annot=do.call(rbind,annlist)
## saveHTML(annot,file="annotation.html",title="Responder genes")


###################################################
### code chunk number 20: sessionInfo
###################################################
sessionInfo()


