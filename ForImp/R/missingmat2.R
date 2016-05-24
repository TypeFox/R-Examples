missingmat2<-function(mat, missing)
{
# element i of vector missing contains the desired number of rows
# with i missing values
n<-dim(mat)[[1]]
m<-dim(mat)[[2]]
# random permutation
s<-sample(n,n,replace=FALSE)
# cumnum[z] gives the number of rows with a number of missing values
# between 1 and z
cumnum<-cumsum(missing)
cumnum<-c(0,cumnum)
for(k in 1:length(missing))
{
num<-missing[k]
if(num>0)
{
for(h in (cumnum[k]+1):cumnum[k+1])
{
# choose of the columns where to put the missing values
ss<-sample(m,k,replace=FALSE)
# s[h] contains the indexes of the rows where to put the missing values
mat[s[h],ss]<-NA
}
}
}
mat
}