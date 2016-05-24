qgmod <-
function(Y,Ped,Model,Cross){
   gdata=genmod.data(Y,Ped,Model,Cross)
   res=mq0(gdata)
   #res$FixedEffect=res$FixEffect
   return(res)
}
