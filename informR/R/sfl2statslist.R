sfl2statslist<-function(sformlist,type=1){
   outpm<-vector("list",length=length(sformlist))
   if(type==1){
      for(i in 1:length(sformlist)){
         outpm[[i]]$global<-sformlist[[i]]
      }
   }
   if(type==2){
      for(i in 1:length(sformlist)){
         outpm[[i]]$local<-sformlist[[i]]
      }
   }
   names(outpm)<-names(sformlist)
   return(outpm)
}
