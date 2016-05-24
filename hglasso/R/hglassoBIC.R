hglassoBIC <- function(x,S,c=0.2){
	
n <- x$n	
v <- length(x$hubind)	
Zcard <- (sum(abs(x$Z)!=0)-x$p)/2
Vcard <- x$V+t(x$V)
diag(Vcard) <- 0
Vcard <- sum(Vcard!=0)/2

return(list(BIC=n*(-log(det(x$Theta))+sum(diag(S%*%x$Theta)))+log(n)*Zcard + log(n)*(v+c*(Vcard-v))))

}