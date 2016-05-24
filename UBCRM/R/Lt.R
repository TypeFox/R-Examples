Lt <-
function(a, data, sgl){
npt<- data$npt
ndlt<- data$ndlt
sapply(a,FUN=function(a){prod((psit(sgl,a)^ndlt)*((1-psit(sgl,a))^(npt-ndlt)))})
}
