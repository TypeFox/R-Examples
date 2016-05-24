is.empty <-
function(x, mode=NULL){
     if (is.null(mode)) mode <- class(x)
     identical(vector(mode,1),c(x,vector(class(x),1)))
}
