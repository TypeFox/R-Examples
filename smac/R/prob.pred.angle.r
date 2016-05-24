prob.pred.angle = function (new.x, beta, beta0, k, loss)
{

XI=XI.gen(k)

f.matrix = t(t(new.x %*% beta)+beta0)

inner.matrix=matrix(0,nrow(new.x),k)

for (ii in 1:k)
{
inner.matrix[,ii] = apply(t(t(f.matrix)*XI[,ii]),1,sum)
}

if (loss=="logi") {fir.der = firlogi}
if (loss=="boost") {fir.der = firboost}
if (loss=="psvm") {fir.der = firpsvm}

fir.der.matrix = 1 / matrix(as.vector( fir.der(inner.matrix)),nrow(new.x),k)

temp=apply( fir.der.matrix,1,sum )

z =  fir.der.matrix / temp

if (loss=="psvm")
	{
		for (i in 1:nrow(z))
			{
			if (any(z[i,]<0)) {z[i,] = (z[i,] - min(z[i,])) / sum((z[i,] - min(z[i,])))}
			}
	}

return(z)

}
