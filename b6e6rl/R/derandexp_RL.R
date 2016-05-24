derandexp_RL <-
function(P,F,CR,expt){
N <- length(P[,1])
d <- length(P[1,])-1
y <- P[expt[1],1:d]
vyb <- nahvyb_expt(N,3,expt)
r123 <- P[vyb,]
fmin <- min(r123[,d+1])
indmin <- which.min(r123[,d+1]) 
r1 <- r123[indmin,1:d]
r123 <- r123[-indmin,]
r2 <- r123[1,1:d]
r3 <- r123[2,1:d]
v <- r1+F*(r2-r3)
L <- trunc(d*runif(1))
change <- L
position <- L
while ((runif(1) < CR) && (length(change) < d)){
	position <- position + 1
	if (position <= d){
		change[length(change)+1] <- position
	}else{
		change[length(change)+1] <- position %% d
	}
}	

y[change]=v[change]
return(y)		
}
