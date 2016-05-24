


cleanclust<- function(marker, nafrac=0.2, mafb=0.1, corbnd=0.5, method="complete"){

## This is a simplified version of the Hclust method described in the paper
## "Characterization of Multilocus Linkage Disequilibrium" by Rinald, Bacanu, Devlin, Sonpar, Wasserman and Roeder
## The R code for the original Hclust package can be find at http://www.wpic.pitt.edu/WPICCompGen/hclust/hclust.htm,
## which provides more functionalities.

## Functions:
## remove markers with proportions of missing values higher than nafrac
## impute missing entries with population average
## remove markers with minor allele frequency below mafb
## if necessary, re-encode the marker so that minor allele is encoded as 1 and major allele as 0.


if (!is.matrix(marker) & !is.data.frame(marker)) stop("marker must be a matrix or data.frame")

if (is.null(colnames(marker))) stop("Please give unique names to each marker")

## remove markers with high proportion of missing values
## impute missing values with column mean
nna<- sum(is.na(marker))
if (nna>0){

  cna<- apply(marker,2, .sumna)
  marker1<- marker[cna/(dim(marker)[2])< nafrac]
  marker1<- apply(marker1, 2, .cimput)


} else {marker1<- marker}


if (min(marker1)<0 | max(marker1)>1) stop("marker should be encoded as 0 or 1 according to the presence of the minor alleles")


af<- apply(marker1,2, mean)
if (sum(af<0 | af>1)>0) stop("each marker must taken value between 0 and 1")

marker1<- marker1[, af>= mafb & af<= 1-mafb]

af1 <- af[ af>= mafb & af<= 1-mafb]

if (sum(af1>0.5)>0){
flip<- colnames(marker1)[af1>0.5]
marker1[af1>0.5]<- 1-marker1[af1>0.5]
cat("Warning: Some markers have minor frequency allele encoded as 0, these markers have been re-encoded such that the minor allele is encoded as 1 and major allele as 0.  The names of re-encoded markers are stored in the vector flip.","\n") 
}else{
flip<-NULL}

##calculate correlation matrix
cormat1<- cor(marker1, method="pearson")

# Perform Hierarchical clustering with distance = 1-cor^2
modhc = hclust(as.dist(1-cormat1^2),method=method)

###  Cut tree based on hcbound on cor^2
clusters = cutree(modhc, h=1-corbnd)
sequ<- 1:length(clusters)
block.ID<- rep(NA,length(clusters))
CorMean <- rep(NA,length(clusters))


allc = table(clusters)
singclust = as.numeric(names(allc)[allc==1])

numsing <- length(singclust)

if(numsing == 0){
  tagSNP.sing = NULL
}else{
  tagSNP.sing = sequ[is.element(clusters,singclust)]
  block.ID[tagSNP.sing] = seq(numsing)
  CorMean[tagSNP.sing] = 1
 
}



### Find tag SNPs for each block cluster
blocks = setdiff(as.numeric(names(allc)),singclust)
tagSNP.bl = NULL      # the list of tagSNPs from the block clusters

if(length(blocks)>0){
  for(i in 1:length(blocks)){
    blSNP <- sequ[is.element(clusters,blocks[i])]
    cor1 <- cormat1[blSNP,blSNP]^2
    a = apply(cor1,1,mean)
    block.ID[blSNP] = i + numsing
    CorMean[blSNP] = a

     # find best correlated SNP
      indexSNP =  blSNP[a==max(a)]

      # choose the middle SNP if there are several optimal SNPs
      indexSNP = indexSNP[floor(length(indexSNP)/2 +.5)]

      tagSNP.bl = c(tagSNP.bl, sequ[indexSNP])
      
  }
}



tagSNP <- c(tagSNP.sing,tagSNP.bl)

newmarker<- marker1[,tagSNP]

res<- list(newmarker=newmarker, flip=flip, tagged=tagSNP)
return(res)

}





.cimput<- function(x){
## calculate column mean for imputation
if(sum(is.na(x))>1){
  cmean<- mean(x,na.rm=T)
  x[is.na(x)]<- cmean
}
return(x)
}

.sumna<- function(x){
## count the number of missing values
cna<- sum(is.na(x))
cna
}






