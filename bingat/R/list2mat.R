list2mat <-
function(x, type="adjMatrix"){
if(tolower(type) == "adjmatrixlt"){
y <- matrix(0, (((dim(x[[1]])[1]*dim(x[[1]])[1]) - dim(x[[1]])[1])/2), length(x))
for(i in 1:length(x))
y[,i] <- x[[i]][lower.tri(x[[i]])]
}else{
y <- matrix(0, dim(x[[1]])[1]*dim(x[[1]])[1], length(x))
for(i in 1:length(x))
y[,i] <- as.vector(x[[i]])
}
return(y)
}
