zeropad <-
function(x){

x <- cbind(x, rep(0, length(x)))

as.vector(t(x))

}
