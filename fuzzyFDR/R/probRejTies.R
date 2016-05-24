"probRejTies" <-
function(p1,p2,r.minus,len,ntests,alpha=0.05){

#############################
#
###  p2, p1 are scalars
#
###  calculates T_(len-s),len = T_k,len
###  returns pi and eta for the interval def. by p2,p1
#
#############################

tsum <- 0
expnum <- 0

t <- alpha/(ntests*(p2-p1))
q <- ((r.minus:(r.minus+len-1))*(alpha/ntests)-p1)/(p2-p1)
#print(paste("q = ",q))
#print(paste("len = ",len))
#print(paste("r.min = ",r.minus))

## k0 is the first k with alpha_k > p1
##k0 <- max(1,ceiling( ntests*p1/alpha - (r.minus+1) ))
#### CHANGED CALC. (need to find 1st index where q>=0)
k0 <- sum(q<=0)+1
#print(paste("k0 = ",k0))

############################# loop over k

## only do this loop if k0<=len 
## this fn may be called with no alpha_k > p1, then prob.rej=0
if(k0<=len){

for(k in k0:len){

s <- len-k

if(s==0){
	term <- q[len]**len
}
else if(s==1){
	term <- len*(1-q[len])*q[len-1]**(len-1)
}
else{
	## alloc matrices must already exist
	text <- paste("nn <- alloc",s,sep="")
	eval(parse(text=text))
	nfact <- factorial(nn)
	nalloc <- dim(nn)[1]

	#print(nn)
	#print(nfact)
	#print(paste("s=",s))
	#print(paste("nalloc=",nalloc))

	### i labels allocations
	term <- 0
	for(i in 1:nalloc){
	  logterm <- s*log(t) + nn[i,s]*log((1-q[len])/t) - sum(log(nfact[i,]))
	  term <- term + exp(logterm)
	}
	term <- term * q[len-s]**(len-s) * 
              factorial(len)/factorial(len-s)
}

#print(paste("l=",len))
#print(paste("k=",k))
#print(paste("term=",round(term,d=3)))

tsum <- tsum + term
expnum <- expnum + k*term
}

tzero <- 1-tsum
piprob <- expnum/len
}
else{
tzero <- 1
piprob <- 0
}

return(list(tzero=tzero,piprob=piprob))

}

