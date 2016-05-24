getNumNodes <-
function(data, type){
if(missing(data) || missing(type))
stop("data and/or type is missing.")

data <- as.matrix(data)

if(tolower(type) == "adjmatrix"){
nodes <- sqrt(nrow(data))
}else if(tolower(type) == "adjmatrixlt"){
c <- nrow(data)
nodes <- (1+sqrt(1+8*c))/2
}else if(tolower(type) == "diag"){
nodes <- nrow(data)
}else{
stop(sprintf("%s is unknown.", as.character(type)))
}

return(nodes)
}
