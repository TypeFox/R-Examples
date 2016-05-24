##' @export
`WriteMultiSV.default` <- function(MultiData, File)
{
  if(length(MultiData[,1]) >0){
    cat(paste(   "#chr","Source","SVType","Start","End","Pval",".",".","Info","\n",sep="\t"),file=File)   
    for(i in 1:length(MultiData[,1])){
      
      if(MultiData[i,6] <= 0.05){
            if(MultiData[i,5] > 0){
                  cat(MultiData[i,1],paste("MultiSV"), paste("DUP"),MultiData[i,2],MultiData[i,3],format.pval(c(MultiData[i,6]),digits=3),paste(".",".",sep="\t"),paste("SVSize=",MultiData[i,4],";SVlog2=",sprintf("%0.3g",MultiData[i,5] ), sep=""),sep="\t",file=File, fill=TRUE, append=TRUE)     
            } else if (MultiData[i,5] < 0){
                  cat(MultiData[i,1],paste("MultiSV"), paste("DEL"),MultiData[i,2],MultiData[i,3],format.pval(c(MultiData[i,6]),digits=3),paste(".",".",sep="\t"),paste("SVSize=",MultiData[i,4],";SVlog2=",sprintf("%0.3g",MultiData[i,5] ), sep=""),sep="\t",file=File, fill=TRUE, append=TRUE)     
            } 
        
      }
    }
  }  
}
