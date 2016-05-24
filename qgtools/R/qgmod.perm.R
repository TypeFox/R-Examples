qgmod.perm <-
function(Y,Ped,Model,Cross,PermNum=NULL){
   if(is.null(PermNum))PermNum=500
   gdata=genmod.data(Y,Ped,Model,Cross)
   res=genmod.perm(gdata=gdata,au=NULL,PermNum=PermNum)
   return(res)
}
