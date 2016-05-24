estGStar <-
function(data){
if(missing(data))
stop("data is missing.")

mean <- rowSums(data)/ncol(data)
gstar <- 1*(mean > .5)

return(gstar)
}
