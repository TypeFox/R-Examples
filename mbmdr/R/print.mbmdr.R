print.mbmdr <-
function(x,...){
 if(!is.null(x$call$output)){
 	 cat("Output readed from file: ",x$call$output,"\n\n")
 	 x$result <- read.table(x$call$output,header=TRUE,sep=";",stringsAsFactors=FALSE)
 }
 
 if(nrow(x$result)==0){
 	cat("No interaction models found:\n")
 	return(NULL)
 }

 ind <- order(pmax(x$result$WH,x$result$WL,na.rm=TRUE),decreasing=TRUE)
 aux <- x$result[ind,]
 print(aux[order(aux$MIN.P),],digits=4,row.names=FALSE)
}

