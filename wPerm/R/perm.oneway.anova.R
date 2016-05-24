perm.oneway.anova <-
function(x, y, trim = 0, ford = NULL, R = 9999)
{
response <- all.names(substitute(x))
if (length(response)>1) response <- response[[3]]
factor <- all.names(substitute(y))
if (length(factor)>1) factor <- factor[[3]]
y <- as.factor(y)
if (!is.null(ford)) y <- factor(y,levels(y)[ford])
Mean <- tapply(x,y,mean)
n <- tapply(x,y,length)
SD <- tapply(x,y,sd)
trim.vector <- function(x)
 {
  n <- length(x)
  lo <- floor(n * trim) + 1
  hi <- n + 1 - lo
  x <- sort.int(x, partial = unique(c(lo, hi)))[lo:hi]
  x
 }
F.trim <- function(x,y) 
 {
  l <- tapply(x,y,trim.vector)
  d <- stack(l)
  oneway.test(d$values~d$ind,var.equal=TRUE)$statistic
 }
obs <- F.trim(x,y)
z <- vector(length=R)
for (i in 1:R)
{
 u <- sample(y)
 z[i] <- F.trim(x,u)
}
p <- signif((sum(z >= obs)+1)/(R+1),digits=3)
pv <- c((p>=0.001)&(p<=0.999),(p<0.001),(p>0.999))
pt <- c(p,"P < 0.001","P > 0.999")
p.value <- pt[pv]
# FOR THE RESULTS
stat.name <- "F.trim"
results <- 
 list(Perm.values=z,
      Header=c(paste("RESULTS OF PERMUTATION ",100*trim,"% TRIMMED ONE-WAY ANOVA\n",sep=""),
      paste("BASED ON",R,"REPLICATIONS ")),
      Response=response,Factor=factor,Levels=levels(y),n=n,Mean=Mean,SD=SD,
      Statistic=stat.name,Observed=obs,P.value=p.value,p.value=p,Trim=trim)
class(results) <- "perm.oneway.anova"  # permutation, one-way ANOVA.
results
}
