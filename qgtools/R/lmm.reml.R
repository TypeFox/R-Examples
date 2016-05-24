lmm.reml <-
function(formula,data=list(),criterion=NULL){
   if(is.null(criterion))criterion=1e-3
   if(is.null(data))gdata=mixed.data(formula)
   else gdata=mixed.data(formula,data)
   res=reml0(gdata,criterion)
   #res$FixedEffect=res$FixEffect
   return(res)
}
