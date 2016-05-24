`fast_parameters`<-function(minimum, maximum, names=paste(sep="", "P", 1:n),factor=1,logscale=rep(FALSE, n), cukier= TRUE, reorder=1:n){
  n<-NROW(minimum)
  if(n!=NROW(maximum)) stop("Expecting minimum and maximum of same size")
  if(n!=NROW(names)) stop("Expecting minimum and names of same size")
  if(n!=NROW(logscale)) stop("Expecting minimum and logscale of same size")

  toreturn <- S(m=n, par.names = names, cukier=cukier, factor=factor, reorder=reorder)
  for(i in 1:n){
      toreturn[,i]<-rerange(toreturn[,i], minimum[i],maximum[i])
      if(logscale[i]){
         toreturn[,i]<-10^toreturn[,i]   
      }
  }
  return(toreturn)
      
}
