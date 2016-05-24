get.w.remove.j <-
function(w,j){
trans.w=t(w)
trans.w[,j]=0
return(as.vector(trans.w))

}
