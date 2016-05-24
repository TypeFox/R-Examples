Lp <-
function(a, data, sgl){
npt<- data$npt
ndlt<- data$ndlt
sapply(a,FUN=function(a){prod((psip(sgl,a)^ndlt)*((1-psip(sgl,a))^(npt-ndlt)))})
}
