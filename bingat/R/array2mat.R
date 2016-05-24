array2mat <-
function(x, type="adjMatrix"){
if(tolower(type) == "adjmatrixlt"){
y <- matrix(0, (((dim(x)[1]*dim(x)[1])-dim(x)[1])/2), dim(x)[3])
for(i in 1:dim(x)[3])
y[,i] <- x[,,i][lower.tri(x[,,i])]
}else{
y <- matrix(0, dim(x)[1]*dim(x)[1], dim(x)[3])
for(i in 1:dim(x)[3])
y[,i] <- x[,,i]
}
return(y)
}
