tr <-
function(a){

p=nrow(a)
q=ncol(a)
if (p==q){
	t=sum(diag(a))
}
return(t)
}
