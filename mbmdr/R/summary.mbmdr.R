summary.mbmdr <-
function(object,sig.level=0.05,...){
 if(!is.null(object$call$output)){
 	 cat("Output readed from file: ",object$call$output,"\n")
 	 object$result <- read.table(object$call$output,header=TRUE,sep=";",stringsAsFactors=FALSE)
 }

 if(nrow(object$result)==0) return(NA)

 ind <- order(pmax(object$result$WH,object$result$WL,na.rm=TRUE),decreasing=TRUE)
 aux <- object$result[ind,]
 ind <- order(aux$MIN.P)
 sig <- aux$MIN.P<=sig.level
 if(sum(sig)==0) out <- aux[ind[1:5],] 
 else out <- aux[ind[1:sum(sig)],]
 class(out) <- c("summary.mbmdr","data.frame")
 return(out)
}

