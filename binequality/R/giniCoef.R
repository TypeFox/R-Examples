giniCoef <-
function(seq.i,samps.i){
  ginTops<-c()
  Ns<-seq.i*100
  for(i in 1:length(Ns)){
    ginTop.i<-(2*Ns[i]-max(Ns)-1)*samps.i[i]
    ginTops<-c(ginTops,ginTop.i)
  }#end for i
  gin.out<-sum(ginTops)/((max(Ns)^2)*mean(samps.i))
  return(gin.out)
}
