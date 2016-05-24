`sq.array` <- 
function(x) 
  array(x <- as.matrix(x), c(n<-sqrt(dim(x)[1]), n, dim(x)[2]))
