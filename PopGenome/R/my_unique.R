my_unique <- function(hapldata){

if(is.data.frame(hapldata)){hapldata <- as.matrix(hapldata)}

dub       <- duplicated(hapldata,fromLast=FALSE)
coresets  <- hapldata[!dub,,drop=FALSE]
temp      <- which(!dub) # ids of haplotypes

idx <- vector(,dim(hapldata)[1])

for(xx in 1:dim(coresets)[1]){
   comp1   <- coresets[xx,]
   ss      <- apply(hapldata,1,function(x){return(all(comp1==x))})
   idx[ss] <- xx
}
numHap <- dim(coresets)[1]
sizHap <- vector(,numHap)
for(xx in 1:numHap){
sizHap[xx] <- sum(idx==xx)
}

return(list(Matrix=hapldata,uniquematrix=coresets,ids=temp,idx=idx,numHap=numHap,sizHap=sizHap))

}
