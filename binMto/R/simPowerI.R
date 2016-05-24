"simPowerI" <-
function(H0diff, pH1, n, n.sim=1000 , conf.level=0.95, alternative="two.sided", method="Add4", 
adj="Dunnett")

{


###### interval calculation

CIfun<-function(n, x, quant, method)
{

xC <- x[1]
xT <- x[-1]

nC <- n[1]
nT <- n[-1]

upper <- numeric(length=k)
lower <- numeric(length=k)


if (method=="Add4")
{
 for (i in 1:k) 
  {
   temp <- Add4(nx=nT[i], X=xT[i], ny=nC, Y=xC, quantile=quant, alternative=alternative)
   lower[i] <- temp$conf.int[1]
   upper[i] <- temp$conf.int[2]
   
  }
}

if (method=="Add2")
{
 for (i in 1:k) 
  {
   temp <- Add2(nx=nT[i], X=xT[i], ny=nC, Y=xC, quantile=quant, alternative=alternative)
   lower[i] <- temp$conf.int[1]
   upper[i] <- temp$conf.int[2]
  
  }
}

if (method=="NHS")
{
 for (i in 1:k) 
  {
   temp <- NHS(nx=nT[i], X=xT[i], ny=nC, Y=xC, quantile=quant, alternative=alternative)
   lower[i] <- temp$conf.int[1]
   upper[i] <- temp$conf.int[2]
   
  }
}

if (method=="Wald")
{
 for (i in 1:k) 
  {
   temp <- Wald(nx=nT[i], X=xT[i], ny=nC, Y=xC, quantile=quant, alternative=alternative)
   lower[i] <- temp$conf.int[1]
   upper[i] <- temp$conf.int[2]
   
  }
}

return(list(lower=lower, upper=upper))
}

# end CIfun


# calculated intervals for n.sim random settings
# an save it in the matrix CImat

H1diff <- pH1[-1]-pH1[1]
k <- length(H0diff)

quant <- Mtoquant(nc=n[1], nx=n[-1], pc=pH1[1], px=pH1[-1], conf.level=conf.level, adj=adj, alternative=alternative)

CImatlo<-matrix(ncol=k, nrow=n.sim)
CImatup<-matrix(ncol=k, nrow=n.sim)

for(i in 1:n.sim)
 {
 data<-rbinom(n=k+1, size=n, prob=pH1)

 temp<-CIfun(n=n, x=data, quant=quant, method=method)
 
 CImatlo[i, 1:k] <- temp$lower
 CImatup[i, 1:k] <- temp$upper

 }



if(alternative=="greater")
{

covvec<-apply(CImatlo, MARGIN=1, FUN=function(x){all(x <= H1diff)} )
coverage<-sum(covvec)/n.sim

anyppvec <-apply(CImatlo, MARGIN=1, FUN=function(x){any(x > H0diff)} )
anypp<-sum(anyppvec)/n.sim 

}



if(alternative=="less")
{

covvec <- apply(CImatup, MARGIN=1, FUN=function(x){all(x >= H1diff)} )
coverage <- sum(covvec)/n.sim

anyppvec <- apply(CImatup, MARGIN=1, FUN=function(x){any(x < H0diff)} )
anypp <- sum(anyppvec)/n.sim 

}



if(alternative=="two.sided")
{

covvecup <- apply(CImatup, MARGIN=1, FUN=function(x){all(x >= H1diff)} )
covveclo <- apply(CImatlo, MARGIN=1, FUN=function(x){all(x <= H1diff)} )

coverage <- sum(covvecup*covveclo)/n.sim

CImatts <- cbind(CImatlo, CImatup)

# a simple two-sided power without considering wrong direction 

anyppvec <- apply(CImatts, MARGIN=1,
 FUN=function(x){any(x[1:k] > H0diff) | any(x[(k+1):(2*k)]< H0diff)}
 )
anypp<-sum(anyppvec)/n.sim 

}



out<-list(coverage=coverage,
 any_pair_power=anypp,
 n.sim=n.sim,
 H0diff=H0diff,
 pH1=pH1 )

return(out)
}

