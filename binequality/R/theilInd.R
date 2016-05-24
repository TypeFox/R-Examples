theilInd <-
function(samps){
  theils<-c()
  for(i in samps){
    t.i<-i/mean(samps)*log(i/mean(samps))
    theils<-c(theils,t.i)
  }#end for i
  t.out<-sum(theils)/length(theils)
  #t.out<-t.out/log(length(theils))
  return(t.out)
}
