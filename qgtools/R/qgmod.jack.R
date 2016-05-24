qgmod.jack <-
function(Y,Ped,Model,Cross,JacNum=NULL,JacRep=NULL){
   if(is.null(JacNum))JacNum=10
   if(is.null(JacRep))JacRep=1
   gdata=genmod.data(Y,Ped,Model,Cross)
   res=genmod.jack(gdata,JacNum=JacNum,JacRep=JacRep)
   return(res)
}
