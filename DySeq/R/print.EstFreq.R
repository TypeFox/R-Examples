#'print.state.trans
#'
#'Generates output for EstFrag object, see: \code{\link{EstFreq}}
#'
#'
#'@param x a EstFrag object, that should be printed. See help(EstFreq)
#'@param ... further arguments passed to or from other methods.
#'@export


print.EstFreq<-function(x, ...){

  if(!(class(x)=="EstFreq")) warning("x should be a  EstFreq object!")

  zero.abs<-sum(x[[1]]==0)
  zero.freq<-1-zero.abs/x[[4]]
  zero.m<-mean(x[[1]])
  zero.max<-max(x[[1]])
  zero.min<-min(x[[1]])

  below.abs<-sum(x[[2]]==0)
  below.freq<-1-below.abs/x[[4]]
  below.m<-mean(x[[2]])
  below.max<-max(x[[2]])
  below.min<-min(x[[2]])



  r.names<-paste("Below ", x[[3]], " Entrys", sep="", collapse="")

  output<-matrix(c(round(zero.freq*100,2),round(below.freq*100,2),round(zero.m,2), round(below.m,2), zero.max,below.max, zero.min, below.min),2,4)
  dimnames(output)[[1]]<-c("Zero Entrys" , r.names)
  dimnames(output)[[2]]<-c("percentage", "||num.of cells: mean", "max", "min")



  cat("\n\n Number of cases with expected frequency problems\n\n")
  print(output)

  return(output)

}




