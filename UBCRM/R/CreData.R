CreData <-
function(ndose = 3, dosenames = paste("dose", 1:ndose, sep=" ")){
data <- data.frame(dose = 1:ndose, npt = rep(0,ndose), ndlt = rep(0,ndose), row.names = dosenames)
return(data)
}
