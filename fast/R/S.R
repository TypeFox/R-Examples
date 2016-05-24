`S` <-
function(m, factor=1,cukier=TRUE, par.names=NULL, reorder=1:m){
#calculate the array of S-values
#usually tabled values from
#cukier 75 are used
#if cukier=F table values from
#McRae82 are used
if(cukier){
  omega=freq_cukier(m)
}else{
  omega=freq_mcrae82(m)
}
tab<-s(m,factor,cukier) %o% omega
toreturn <- asin(sin(tab))/pi
toreturn <- data.frame(toreturn[,reorder])
if(!is.null(par.names[1]))
   names(toreturn)<- par.names
return(toreturn)
}

