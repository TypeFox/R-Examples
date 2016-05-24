coprob <-
function(max_inc,j,min_inc,nsite){
    as.numeric(round(chooseZ(max_inc,j) * chooseZ(nsite - max_inc, min_inc - j),0) / round(chooseZ(nsite,min_inc),0))
}
