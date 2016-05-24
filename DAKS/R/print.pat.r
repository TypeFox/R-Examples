print.pat<-function(x, ...){		
cat(x$n, "largest response patterns in the data:", " ")
print(x$response.patterns)
if(is.null(x$states) == FALSE){
cat("\nNumber of times a state occurs in the data:\n")
print(x$states)	
}
}
