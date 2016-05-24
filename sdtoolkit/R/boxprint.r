
boxprint <- function(boxseq,boxes=1:(length(boxseq)-2),reldim=FALSE,outfile=NA,append=FALSE){

  
  if(!is.na(outfile)){
    sink(outfile, split=TRUE,append=append)
  }
  
  
  for (i in boxes){
    
    if(reldim){
      bbounds <- boxseq[[i]]$relbox
    } else { 
      bbounds <- boxseq[[i]]$box
    }
    
    print(boxformat(bbounds,morestats=boxseq[[i]]$morestats,dimlist=boxseq[[i]]$dimlist))
    cat("\n","\n")
  
  }
  
  
  if(!is.na(outfile)){
    sink()
  }

}
