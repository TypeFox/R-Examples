lmm.mq <-
function(formula,data=list()){
   if(is.null(data))gdata=mixed.data(formula)
   else gdata=mixed.data(formula,data)
   res=mq0(gdata)
   #res$FixedEffect=res$FixEffect
   return(res)
}
