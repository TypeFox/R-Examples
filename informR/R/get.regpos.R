`get.regpos` <-
function(gregx){
   outregpos<-list()
   start.tmp<-lapply(gregx,function(x) c(x))
   len.tmp<-lapply(gregx,function(x) attr(x,"",""))
   outregpos$start<-start.tmp
   outregpos$stop<-mapply(function(x,y) x+y-1,start.tmp,len.tmp,SIMPLIFY=FALSE)
   outregpos
}

