depth.radius = function(data, alpha, beta, weight)
{
	data = t(data$y)
	p = dim(data)[2]
	n = dim(data)[1]
	if(class(colnames(data))=="NULL")
	{	
		h = as.matrix(c(0.5, rep(1,(p-2)), 0.5))
	}
	else
	{
		sequence = as.numeric(colnames(data))
		h = as.matrix(c((sequence[2]-sequence[1])/2, (sequence[3:p]-sequence[1:(p-2)])/2, (sequence[p]-sequence[p-1])/2))
	}
	na = ceiling(alpha * n)
	nb = n - floor(beta * n)
	
	r = vector(,n)
    for(i in 1:n)
	{	
		d = sqrt((data - t(matrix(rep(data[i,], n),,n)))^2 %*% h)
		d = sort(d)
		r[i] = d[na]
	}
	
	w = vector(,n)
	d = sort(r)
	rnk_n = apply(matrix(rep(1,n),n,1)%*%matrix(r,1,n) <= matrix(r,n,1)%*%matrix(rep(1,n),1,n),1,mean)
	if(weight == "hard")
	{
		w = ifelse(r<d[nb], 1, 0)
	}
	else	
	{
		a = 0.5
		b = 1 - beta
		w = ifelse(rnk_n <= a, 1, 0)
		index = which(rnk_n >a & rnk_n<=b)
		w[index] = (rnk_n[index]-b)*((1/(a-b))+(rnk_n[index]-a)*(2*rnk_n[index]-(a+b))/(b-a)^3)
	}
	mu = colSums(matrix(rep(w,p),,p)*data)/sum(w)
	robustmedian = apply(data[which(w==1),],2,median)	
	return(list(median = robustmedian, mtrim = mu, weight = w))
}
	
	
	
	