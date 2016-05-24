listm<-function(x){
#
# Store the data in a matrix or data frame in a new
# R variable having list mode.
# Col 1 will be stored in y[[1]], col 2 in y[[2]], and so on.
#
if(is.null(dim(x)))stop("The argument x must be a matrix or data frame")
y<-list()
for(j in 1:ncol(x))y[[j]]<-x[,j]
y
}
