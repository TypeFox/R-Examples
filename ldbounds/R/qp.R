"qp" <-
function(xq,last,nints,yam1,ybm1,stdv){
   hlast <- (ybm1-yam1)/nints
   grid <- seq(yam1,ybm1,length=nints+1)
   fun <- last*pnorm(grid,mean=xq,sd=stdv)
   qp <- 0.5*hlast*(2*sum(fun)-fun[1]-fun[length(fun)]) # This is "trap"
   return(qp)
 }

