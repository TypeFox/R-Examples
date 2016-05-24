dissmutmatrix <-
function(data=NULL,chr=NULL,position,subset=NULL,upper=FALSE){

arguments <- as.list(match.call())

position = eval(arguments$position, data)
chr = eval(arguments$chr, data)

if (is.null(chr)==TRUE){
data <- sort(position)
}
else if (is.null(chr)==FALSE){
data <- data.frame(chr,position)
if (is.null(subset)==FALSE){
data<-subset(data, chr==subset, select = c(position))
}
data<- sort(data$position)
}

mydist<-dist(data,method="euclidean",upper=upper)
mydist<-log10(mydist)
return(mydist)
}
