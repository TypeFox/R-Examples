roulete <-
function(cutpoints){
h <- length(cutpoints)
ss <- sum(cutpoints)
p_min <- min(cutpoints)/ss
cp <- cutpoints[1]
for(i in 2:h){
	cp[i] <- cp[i-1]+cutpoints[i]
}
cp <- cp/ss
res <- 1 + trunc(sum(cp<runif(1)))		
return(list(res=res,p_min=p_min))
}
