aggreg<-function(x,ENV)  # x .. data
{
   # aggregate data - calculate counts for patterns
   notnaidx<-!is.na(x[1,])             # which columns are not missing
   cdat<-as.matrix(x[,notnaidx])       # complete data - remove NA columns from data
   cdatStr<-apply(cdat,1, paste,sep="",collapse="")
   YStr<-apply(as.matrix(ENV$Y[,notnaidx]),1,paste,sep="",collapse="")
   counts<-as.vector(table(factor(cdatStr,levels=unique.default(YStr))))

   # calculate CL vector s
   mpStr<-apply(as.matrix(ENV$Y[,notnaidx]),1,paste,sep="",collapse="")
   s<-as.numeric(unclass(factor(mpStr,levels=unique.default(mpStr))))

   list(counts=counts,notnaidx=notnaidx,s=s)

}
