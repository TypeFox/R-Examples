Ll <-
function(a, data, sgl){
npt<- data$npt
ndlt<- data$ndlt
sapply(a,FUN=function(a){prod((psil(sgl,a)^ndlt)*((1-psil(sgl,a))^(npt-ndlt)))})
}
