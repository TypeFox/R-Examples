negent <-
function(x,K=ceiling(log2(length(x))+1)) {
	p <- table(cut(x,breaks=K))/length(x)
	sigma2 <- sum((1:K)^2*p)-sum((1:K)*p)^2
	p <- p[(p>0)]
	(1+log(2*pi*(sigma2+1/12)))/2+sum(p*log(p))
}
