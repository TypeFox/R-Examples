phi <-
function(a,b){

if ((is.matrix(a)) & (is.matrix(b))){
	p=solve(diag(diag(t(a)%*%a),nrow=ncol(a))^.5)%*%t(a)%*%b%*%solve(diag(diag(t(b)%*%b),nrow=ncol(b))^.5)
} else{
	p=solve((t(a)%*%a)^.5)%*%t(a)%*%b%*%solve((t(b)%*%b)^.5)
}
return(p)
}
