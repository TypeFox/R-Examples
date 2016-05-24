table_import <-
function(file,sep=F,format_corrector=F,...){
  x<-data.frame(1:2)
   if(sep==F){sep=c(";",",","\t",":"," ") }
  for(i in 1:length(sep)){
    
    x<-read.table(file,sep=sep[i],...)
    i<-i+1
    if(ncol(x)>1){break}
  }
  if(ncol(x)==1){cat("specify separator")}
  if(format_corrector){x<-format_corrector(x)}
  return(x)
}
