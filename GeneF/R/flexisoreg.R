`flexisoreg` <-
function(y, x, lambda=0, alpha.location=1, alpha.adjacency=0.5){
o <- order(x, decreasing=FALSE)
x <- x[o]
y <- y[o]

#sum of squares under the null
mu <- mean(y)
ss0 <- sum( (y-mu)^2 )

#estimate under the alternative (reduced isotonic)
g <- get.numbers(x)
score <- double(length(g)-1)
link <- integer(length(g)-1)
mu <- double(length(g)-1)
#print(g)

mu[1] <- mean(y[(g[1]+1):g[2]])
score[1] <- sum( (y[(g[1]+1):g[2]]-mu[1])^2 )
link[1] <- 0
for(i in 2:length(score)){
mu[i] <- mean(y[(g[1]+1):g[i+1]])
score[i] <- sum( (y[(g[1]+1):g[i+1]]-mu[i])^2 )
link[i] <- 0
for(j in 1:(i-1)){
temp1 <- mean(y[(g[j+1]+1):g[i+1]])
temp2 <- score[j] + sum( (y[(g[j+1]+1):g[i+1]]-temp1)^2 )
if(temp2 < score[i]){
if( t1p1(y[(g[j+1]+1):g[i+1]]-lambda, g[i+1]-g[j+1]) < alpha.location){
if( t1p2(y[(g[link[j]+1]+1):g[i+1]], g[j+1]-g[link[j]+1], g[i+1]-g[j+1]) < alpha.adjacency ){
mu[i] <- temp1
score[i] <- temp2
link[i] <- j
}
}
}
}
}
ss1 <- score[length(score)]
statistic <- (ss0/ss1 - 1)*(length(y)-length(g)+1)/(length(g)-1-1)
#print(link)

i <- length(link)
partition <- g[i+1]
while(link[i]>0){
i <- link[i]
partition <- c(g[i+1], partition)
}
partition <- c(g[0+1], partition)
#print(partition)

estimates <- double(length(y))
for(i in 2:length(partition)){
estimates[o][(partition[i-1]+1):partition[i]] <- mean( y[(partition[i-1]+1):partition[i]] )
}

groups <- integer(length(y))
for(i in 2:length(partition)){
groups[o][(partition[i-1]+1):partition[i]] <- i-1
}

return(list(groups=groups, estimates=estimates, statistic=statistic))
}

