vec2mat <-
function(x, type="adjMatrix"){
y <- NULL
nn <- getNumNodes(x, type)

if(tolower(type) == "adjmatrixlt"){
y <- matrix(0, nn, nn)
y[lower.tri(y)] <- x
y <- y + t(y)
}else if(tolower(type) == "diag"){
y <- matrix(0, nn, nn)
diag(y) <- x
}else{
y <- matrix(x, nn, nn)
}

return(y)
}
