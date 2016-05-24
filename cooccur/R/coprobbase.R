coprobbase <-
function(max_inc,j,min_inc,nsite){
    as.numeric(round(choose(max_inc,j) * choose(nsite - max_inc, min_inc - j),0) / round(choose(nsite,min_inc),0))
}
