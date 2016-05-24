`print.WGassociation` <-
function (x, digits = 5, ...) 
{
    if (!inherits(x, "WGassociation")) 
        stop("x must be an object of class 'WGassociation'")
    ans <- attr(x, "pvalues")
    if(ncol(ans)>1){
    if(!is.numeric(ans[,1])){
       ans[, -1] <- round(ans[, -1,drop=FALSE], digits)
       out <- as.matrix(ans)
       out[,1]<-gsub("\\\\","",out[,1,drop=FALSE])
    } else {
          
       out <- as.matrix(round(ans, digits))
    }
    } else {
       out<-gsub("\\\\","",as.matrix(ans))
    }
   
    print(out, quote = FALSE, na.print = "-", ...)
    invisible(ans)
        
}

