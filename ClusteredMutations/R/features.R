features <-
function(data=NULL,chr=NULL,position,refbase,mutantbase,min=6,max=5000){

arguments <- as.list(match.call())

if(min > max ) stop("ERROR: min > max")
if(min < 2 ) stop("ERROR: min < 2 ")
if(max < 2 ) stop("ERROR: max < 2 ")

position = eval(arguments$position, data)
chr = eval(arguments$chr, data)
refbase = eval(arguments$refbase, data)
mutantbase = eval(arguments$mutantbase, data)

if (is.null(chr)==TRUE){
chr<-c(rep("N",length(position)))
}
y<-data.frame(cluster=integer(),chr=character(),position=integer(),
  refbase=character(),mutantbase=character(),stringsAsFactors=FALSE)
chr<-factor(chr,levels=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","Y","M","N"),ordered=TRUE)
z<-data.frame(chr,position,refbase,mutantbase)
z<-z[order(chr,position),]

    clustered <- 1
for (k in 1:length(levels(z$chr))){
x<-subset(z, chr==levels(z$chr)[k], select = c(chr,position,refbase,mutantbase))
i <- 1
l <- 1
fdij<-max+1
while (i<=(nrow(x)-min+1)){
dij <- x[i+min-1,2]-x[i,2]
if (dij <= max){
y<-rbind(y,cbind(clustered,x[i,]))
l <- 0
}
else if (dij > max){
l <- l + 1
if ((fdij <= max) & (l <= min-1)){
y<-rbind(y,cbind(clustered,x[i,]))
}
if ((fdij <= max) & (l == min -1)){
clustered <- clustered +1
}
}
if (l == 0){
fdij<-dij
}
i <- i+1
if ((i>(nrow(x)-min+1)) & (dij <= max)){
for (p in (nrow(x)-min+2): nrow(x)){
y<-rbind(y,cbind(clustered,x[p,]))
}
}
}
   }
   return(y)
}
