showers <-
function(data=NULL,chr=NULL,position,min=6,max=5000){

arguments <- as.list(match.call())

if(min > max ) stop("ERROR: min > max")
if(min < 2 ) stop("ERROR: min < 2 ")
if(max < 2 ) stop("ERROR: max < 2 ")

position = eval(arguments$position, data)
chr = eval(arguments$chr, data)

if (is.null(chr)==TRUE){
chr<-c(rep("N",length(position)))
}

chr<-factor(chr,levels=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","Y","M","N"),ordered=TRUE)
z<-data.frame(chr,position)
z<-z[order(chr,position),]
y<-data.frame(chr=character(),pend=integer(),pstart=integer(),nend=integer(), 
  nstart=integer(),stringsAsFactors=FALSE) 
for (k in 1:length(levels(z$chr))){
x<-subset(z, chr==levels(z$chr)[k], select = c(position))
ii <- 1
while (ii<=(nrow(x)-min+1)){
m <- 0
i <- ii
ii <- ii+1
w<-data.frame(chr=character(),pend=integer(),pstart=integer(),nend=integer(), 
  nstart=integer(),stringsAsFactors=FALSE)
for (j in (i+min-1):nrow(x)){
dij <- x[j,1]-x[i,1]
if (dij <= max){
if (m==0){
w[1,1]<-levels(z$chr)[k]
w[1,2]<-x[j,1]
w[1,3]<-x[i,1]
w[1,4]<-j
w[1,5]<-i
m<-1
}
else if (m==1){
w[1,2]<-x[j,1]
w[1,4]<-j
}
}
else if (dij>max){
if (j > i+min-1){
if (nrow(y)==0){
y<-rbind(y,w)
}
else{
if (y[nrow(y),1]==w[1,1]){
if (y[nrow(y),4]>=w[1,5]){
y[nrow(y),2]<-w[1,2]
y[nrow(y),4]<-w[1,4]
}
else {
y<-rbind(y,w)
}
}
else if (y[nrow(y),1]!=w[1,1]){
y<-rbind(y,w)
}
}
ii<-j-min+1
}
break
}
}
}
   }
   y$distance=y[,2]-y[,3]
   y$number=y[,4]-y[,5]+1
   return(y)
}
