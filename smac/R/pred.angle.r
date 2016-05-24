pred.angle = function (new.x, beta, beta0, k)
{

XI=XI.gen(k)


f.matrix = t(t(new.x%*% beta)+beta0)

inner.matrix=matrix(0,nrow(new.x),k)

for (ii in 1:k)
{
inner.matrix[,ii] = apply(t(t(f.matrix)*XI[,ii]),1,sum)
}

z=apply(inner.matrix,1,pred)

return(z)

}
