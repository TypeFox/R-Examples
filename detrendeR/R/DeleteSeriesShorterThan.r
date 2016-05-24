DeleteSeriesShorterThan = function(rwl, filename, rwlStats=NULL, YEAR=100, SaveShort=FALSE, Rwl=TRUE)
{

nSeries1<-nSeries<-ncol(rwl)
    if (is.null(rwlStats)) {
yr.range = function(x) {
        yr.vec = as.numeric(names(x))
        mask = !is.na(x)
        range(yr.vec[mask])
    }
    series.stats = data.frame(Series = colnames(rwl))
    series.stats$First = apply(rwl, 2, yr.range)[1, ]
    series.stats$Last = apply(rwl, 2, yr.range)[2, ]
    series.stats$Year = series.stats$Last - series.stats$First + 1
	rwlStats <- series.stats
}

       flag<-rwlStats[,4]>=YEAR

        rwl1<-subset(rwl,select=c(flag))
if (ncol(rwl1)==0){ cat("All series are shorter than ", YEAR, ".\n", sep="")
SaveShort=FALSE
Rwl=FALSE
}
if (Rwl){
        samp.depth = apply(rwl1, 1, function(y) sum(!is.na(y)))
        rwl1<-subset(rwl1, samp.depth > 0)
nSeries1<-ncol(rwl1)}
{if (nSeries!=nSeries1){
cat("The following series were shorter than " , YEAR, ":\n", sep="")
Flag<-rwlStats[,4]<YEAR
Flag=colnames(rwl)[Flag]
for (i in 1:length(Flag)) cat(Flag[i],"\n")

      b<-paste(filename, "+", YEAR,".rwl", sep="")
        write.rwl(rwl1, b)

if (SaveShort)
  {
        rwl2<-subset(rwl,select=c(!flag))                               #To save series with less than YEAR
        samp.depth = apply(rwl2, 1, function(y) sum(!is.na(y)))
        rwl2<-subset(rwl2, samp.depth > 0)

        b<-paste(filename, "-", YEAR,".rwl", sep="")            #
        write.rwl(rwl2, b)
  }
}                                              #
else{ if (Rwl)cat(paste("All series are longer than ", YEAR, ".\n", sep=""))
      
}}
cat(rep("=",98) , "\n",sep="")
if(ncol(rwl1)>=1)return(rwl1)
}
#DeleteSeriesShorterThan(rwl, filename="teste", YEAR=41, Rwl=TRUE)->a

