print.mbmdr.PermTest <-
function(x,...){
 if(is.null(x$PermTest)){
 	cat("No interaction models found.\n")
 	return(NULL)
 }
 
 precision <- as.integer(log10(x$n))+2
 out <- x$PermTest[order(x$PermTest$Wmax,decreasing=TRUE),]
 out <- out[order(out$Perm.P),]
 out$Perm.P <- as.character(round(out$Perm.P,digits=precision))
 out[out$Perm.P==0,"Perm.P"] <- paste("< ", round(1/(x$n+1),digits=precision) ,sep="")
 print(out,digits=4,row.names=FALSE)
}

