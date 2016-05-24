commonInterval=function(rwl){
years<-as.numeric(rownames(rwl))
TimeSpan<-range(years)
cat(paste("\nTime span: ",TimeSpan[1], "-", TimeSpan[2], ".\n", sep=""))
rwl.common<-na.omit(rwl)
last<-max(years)
if (nrow(rwl.common)>0){
years<-as.numeric(rownames(rwl.common))
first<-as.numeric(min(years))
cat("\nCommon interval: ",first, "-", last, ".\n", sep="")
}
if(!nrow(rwl.common)>0) {cat("NO COMMON INTERVAL!!!!!\n")} 
}
# common.interval(rwl)
