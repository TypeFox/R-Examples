tunemat.fn <-
function(l.forget,l.noforget,d){
tune.mat<-matrix(unlist(unique(combn(rep(c(l.noforget,l.forget),d),d,simplify=FALSE))),(2^d),d,byrow=TRUE)
return(tune.mat)}

