# compute the cosine between two vectors
vcos <- function(x,y){
    if(length(x)!=length(y)) stop(" x and y have different lengths")
     vnorm <- function(x) as.matrix(x/as.numeric(sqrt(t(x) %*% x)))
     t(vnorm(x)) %*% vnorm(y)
}
