lmm.mq.perm <-
function(formula,data=list(),PermNum=NULL){
   if(is.null(PermNum))PermNum=500
   ##if(is.null(JacRep))JacRep=1
   if(is.null(data))gdata=mixed.data(formula)
   else gdata=mixed.data(formula,data)
   res=genmod.perm(gdata=gdata,au=NULL,PermNum=PermNum)
   return(res)
}
