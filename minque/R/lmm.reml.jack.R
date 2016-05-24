lmm.reml.jack <-
function(formula,data=list(),criterion=NULL){
   if(is.null(criterion))criterion=1e-3
   if(is.null(data))gdata=mixed.data(formula)
   else gdata=mixed.data(formula,data)
   res=genmod.reml.jack(gdata,criterion=criterion)
   #res$FixedEffect=res$FixEffect
   return(res)
}
