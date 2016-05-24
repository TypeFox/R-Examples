perm.ind.spread <-
function(x, y, parameter, stacked = TRUE, variable = NULL,
                            alternative = c("two.sided","less","greater"),
                            R=9999)
{
statistic <- parameter
if (stacked)
  {
   variable <- all.names(substitute(x))
   if (length(variable)>1) variable <- variable[[3]]
   y <- as.factor(y)
   l <- levels(y)
   pop.1 <- l[1]
   pop.2 <- l[2]
   u <- x[y==l[1]]
   v <- x[y==l[2]]
   x <- u
   y <- v
  } 
 else
  {
   pop.1 <- all.names(substitute(x))
   if (length(pop.1)>1) pop.1 <- pop.1[[3]]
   pop.2 <- all.names(substitute(y))
   if (length(pop.2)>1) pop.2 <- pop.2[[3]]
  }
x <- x
y <- y
m <- length(x)
n <- length(y)
obs <- statistic(x)/statistic(y)
u <- stack(list(x=x,y=y))
s <- u$values
t <- u$ind
z <- vector(length=R)
for (i in 1:R)
{
 v <- sample(t)
 z[i] <- statistic(s[v=="x"])/statistic(s[v=="y"])
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
ac <- c("different.spread","smaller.spread","larger.spread")
alt <- ac[tc==alternative]
# FOR THE RESULTS
stat.name <- paste("ratio.",all.names(substitute(parameter)),sep="");
results <- 
 list(Stacked=stacked,Perm.values=z,
      Header=c("RESULTS OF PERMUTATION INDEPENDENT TWO-SAMPLE SPREAD TEST\n",
      paste("BASED ON",R,"REPLICATIONS")),
      Variable=variable,Pop.1=pop.1,Pop.2=pop.2,n.1=m,n.2=n,Statistic=stat.name,
      Observed=obs,Null="identical",Alternative=alt,P.value=p.value,p.value=p)
class(results) <- "perm.ts.ind"  # permutation, two-sample, independent.
results
}
