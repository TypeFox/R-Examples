ESpredict <-
function(object,c){
EPS <- 100*(.Machine$double.eps)
outputbeta <- matrix(0,length(c),ncol(object$beta))
for (i in 1:length(c))
{
if (c[i] < object$c1[1] && c[i] > EPS)
{
if(min(object$c1) >  c[i])
{
index<- c(max(which(object$c1 >= c[i])),length(object$c1))
}
else
{
index<- c(max(which(object$c1 >= c[i])),min(which(object$c1 <= c[i])))
}
first <- object$c1[index[1]]
second <-object$c1[index[2]]

diff <- (first-second)
if (diff < EPS)
{
prop <- 0
}
else
{
prop <- (object$c1[index[1]]-c[i])/diff
}
outputbeta[i,] <- object$beta[index[1],] + prop*(object$beta[index[2],] - object$beta[index[1],])

}
if (c[i] >= object$c1[1])
{
outputbeta[i,] <- object$beta[1,]
}
if (c[i] < EPS)
{
outputbeta[i,] <- object$beta[length(object$c1),]
}
}
return(outputbeta)
}
