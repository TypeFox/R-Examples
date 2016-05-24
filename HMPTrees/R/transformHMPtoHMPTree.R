transformHMPtoHMPTree <-
function(data){
if(missing(data))
stop("A valid data set is required.")

data <- as.data.frame(t(data))

return(data)
}
