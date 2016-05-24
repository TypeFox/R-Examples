perm.paired.loc <-
function(x, y, parameter, variable = NULL,
                            alternative = c("two.sided","less","greater"),
                            R=9999)
{
statistic <- parameter
pop.1 <- all.names(substitute(x))
if (length(pop.1)>1) pop.1 <- pop.1[[3]]
pop.2 <- all.names(substitute(y))
if (length(pop.2)>1) pop.2 <- pop.2[[3]]
obs <- statistic(x)-statistic(y)
n <- length(x)
z <- vector(length=R)
for (i in 1:R)
{
 b <- rbinom(n,1,0.5)
 u <- b*x+(1-b)*y
 v <- (1-b)*x+b*y
 z[i] <- statistic(u)-statistic(v)
}
if (identical(alternative,c("two.sided","less","greater")))
  alternative <- "two.sided"
ltp <- (sum(z<=obs)+1)/(R+1)
rtp <- (sum(z>=obs)+1)/(R+1)
tc <- c("two.sided","less","greater")
pc <- c(2*min(ltp,rtp),ltp,rtp)
p <- signif(pc[tc==alternative],digits=3)
pv <- c((p>=0.001)&(p<=0.999),(p<0.001),(p>0.999))
pt <- c(p,"P < 0.001","P > 0.999")
p.value <- pt[pv]
ac <- c("shifted","shifted.left","shifted.right")
alt <- ac[tc==alternative]
# FOR THE RESULTS
stat.name <- paste("diff.",all.names(substitute(parameter)),sep="")
results <- 
 list(Perm.values=z,
      Header=c("RESULTS OF PERMUTATION PAIRED LOCATION TEST\n",
      paste("BASED ON",R,"REPLICATIONS")),
      Variable=variable,Pop.1=pop.1,Pop.2=pop.2,n=n,Statistic=stat.name,
      Observed=obs,Null="identical",Alternative=alt,P.value=p.value,
      p.value=p)
class(results) <- "perm.paired.loc"  # permutation, paired, location.
results
}
