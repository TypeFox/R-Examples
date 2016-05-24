aggregR<-function(x,ENV)  # x .. data
{
   notnaidx<-!is.na(x[1,])             # which columns are not missing
   cdat<-as.matrix(x[,notnaidx])       # complete data - remove NA columns from data
   cdatStr<-apply(cdat,1,function(x) paste(x,sep="",collapse=""))

   # generate all PC patterns from all ranks patterns
   nobj<-ENV$Rnobj
   Rpatt<-perm(nobj,nobj)
   pc.patt<-NULL
   for (j in 2:nobj){
     for (i in 1:(j-1)){
       pc.patt<-cbind(pc.patt,as.numeric(Rpatt[,i]>Rpatt[,j]))
     }
   }
   pc.patt<-as.matrix(pc.patt[,notnaidx])      # corresponding complete patterns - remove NA columns from patterns

   # counts
   RpattStr<-apply(pc.patt,1,function(x) paste(x,sep="",collapse=""))
   counts<-tabulate(match(cdatStr,unique(RpattStr)),length(unique(RpattStr)))

   # CL vector
   s<-as.numeric(unclass(factor(RpattStr,levels=unique.default(RpattStr))))

   list(counts=counts,notnaidx=notnaidx,s=s)

}
