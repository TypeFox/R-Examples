`print.re` <- function(x,...){
 X <- list()
 if(x$inputN>1){
   for(i in 1:x$inputN){
     X[[i]] <- as.numeric(x[[i]]$sigTests)
   }
 } else {
   X$sigTests <- as.numeric(x$sigTests)
 }
 print(X,...)
}