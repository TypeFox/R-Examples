makeFitComb <-
function(distFitsList){
  for(i in 1:length(distFitsList)){
    if(i==1){
      out<-distFitsList[[i]]$datOut
    }else{
    out<-rbind(out, distFitsList[[i]]$datOut)
    }#end if/else i==1
  }#end for i
  return(out)
}
