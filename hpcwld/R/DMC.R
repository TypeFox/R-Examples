DMC <-
function(X, Y) {
if(is.null(X) || is.null(Y))
	stop("'X', 'Y' must be defined!")
if(!is.numeric(X) || !is.numeric(Y))
	stop("Only numeric vectors allowed!")
n=min(length(X),length(Y))
Yl=Y[X<median(X)]
Yu=Y[X>=median(X)]
sum=0
Y=sort(Y)
for(i in 1:n)
	if(abs(ecdf(Yu)(Y[i])-ecdf(Yl)(Y[i]))>2*sqrt(i)/n)
		sum=sum+1
res=sum/n
res
}
