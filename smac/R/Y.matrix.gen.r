Y.matrix.gen=function(y.train, k)
{
nobs=length(y.train)

Y.matrix = matrix(0,nobs,k-1)

XI=XI.gen(k)

for (ii in 1:nobs)
{
Y.matrix[ii,] = XI[,y.train[ii]]
}
return(Y.matrix)

}
