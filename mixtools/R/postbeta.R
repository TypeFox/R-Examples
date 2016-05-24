post.beta <- function(y, x, p.beta, p.z)
{
N = length(y)
k = ncol(p.z)
g = apply(p.z, 1, function(i) (1:length(i))[i == max(i)])


min.int=min(sapply(1:N, function(i) min(p.beta[[i]][1,])))
max.int=max(sapply(1:N, function(i) max(p.beta[[i]][1,])))
min.slp=min(sapply(1:N, function(i) min(p.beta[[i]][2,])))
max.slp=max(sapply(1:N, function(i) max(p.beta[[i]][2,])))


par(mfrow=c(2,2))



#All posterior regression lines.
plot(c(min(sapply(x,min)),max(sapply(x,max))),c(min(sapply(y,min)),max(sapply(y,max))),type="n",xlab="x-values",ylab="y-values",main="Data and Posterior Regression Lines")
sapply(1:length(x),function(i) points(x[[i]], y[[i]]))
for(j in 1:k){
sapply(1:length(p.beta), function(i) abline(coef=p.beta[[i]][,j], col=(j+1)))
}


#Posterior regression lines chosen according to the membership probabilities.
plot(c(min(sapply(x,min)),max(sapply(x,max))),c(min(sapply(y,min)),max(sapply(y,max))),type="n",xlab="x-values",ylab="y-values",main="Data and Most Probable Posterior Regression Lines")
sapply(1:length(x),function(i) points(x[[i]], y[[i]]))
sapply(1:length(p.beta), function(i) abline(coef=p.beta[[i]][,g[i]], col=(g[i]+1)))



#All posterior beta values.
plot(c(min.int,max.int),c(min.slp,max.slp),type="n",xlab="Posterior Intercepts",ylab="Posterior Slopes",main="All Posterior Regression Coefficients")
for(j in 1:k){
sapply(1:length(p.beta), function(i) points(t(p.beta[[i]][,j]), col=(j+1)))
}


#Posterior beta values chosen according to the membership probabilities.
min.max=c(min(sapply(p.beta,min)),max(sapply(p.beta,max)))
plot(c(min.int,max.int),c(min.slp,max.slp),type="n",xlab="Posterior Intercepts",ylab="Posterior Slopes",main="Most Probable Posterior Regression Coefficients")
sapply(1:length(p.beta), function(i) points(t(p.beta[[i]][,g[i]]), col=(g[i]+1)))


}


