imd <-
function(data=NULL,chr=NULL,position,extra=NULL){

arguments <- as.list(match.call())

chr = eval(arguments$chr, data)
position = eval(arguments$position, data)
extra = eval(arguments$extra, data)

if (is.null(chr)==TRUE){
chr<-c(rep("N",length(position)))
}
if (is.null(extra) == TRUE){
data <- data.frame(chr,position)
}
else if (is.null(extra) == FALSE){
data <- data.frame(chr,position,extra)
}

data$chr <- factor(data$chr,levels=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","Y","M","N"))
data <- data[order(data$chr,data$position),]
data$number <- 1:nrow(data)
data$distance <- ave(data$position,factor(data$chr), FUN=function(x) c(diff(x),NA))
data$log10distance<-round(log10(data$distance),digits=3)
data<-data[complete.cases(data),]
return(data)
}
