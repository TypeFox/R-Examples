residu = function(xp, data_x, data_y)
{
    dm = ncol(data_x)
    data_num = nrow(data_x)
    hn = data_num^(-1/(4 + dm))
    tau2 = h = vector(, dm)
	for(i in 1:dm)
	{
		tau2[i] = exp(xp[i])
		h[i] = sqrt(tau2[i]) * hn
	}
	hprod = prod(h)
	cont = exp(-0.5 * dm * log(2.0 * pi))
	data_e = vector(,data_num)
	for(i in 1:data_num)
	{
		temp = (sweep(data_x[-i,], 2, data_x[i,])/h)^2
		weight = cont * exp(-0.5 * apply(temp, 1, sum))/hprod
		suma = sum(weight * data_y[-i])
		sumb = sum(weight)
		data_e[i] = data_y[i] - suma/sumb
	}
    return(data_e)
}