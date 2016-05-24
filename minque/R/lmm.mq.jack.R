lmm.mq.jack <-
function(formula,data=list(),JacNum=NULL,JacRep=NULL){
   if(is.null(JacNum))JacNum=10
   if(is.null(JacRep))JacRep=1
   if(is.null(data))gdata=mixed.data(formula)
   else gdata=mixed.data(formula,data)
   res=genmod.jack(gdata,JacNum,JacRep)
   return(res)
}
