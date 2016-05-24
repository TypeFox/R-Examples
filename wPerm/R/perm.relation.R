perm.relation <-
function(x, y, method=c("pearson", "kendall", "spearman"),
                          alternative = c("two.sided","less","greater"), R = 9999)
{
Cor <- function(x,y) cor(x,y,method=method)
var1.name <- all.names(substitute(x))
if (length(var1.name)>1) var1.name <- var1.name[[3]]
var2.name <- all.names(substitute(y))
if (length(var2.name)>1) var2.name <- var2.name[[3]]
obs <- Cor(x,y)
n <- length(x)
z <- vector(length=R)
for (i in 1:R)
{
 w <- sample(y)
 z[i] <- Cor(x,w)
}
ltp <- (sum(z<=obs)+1)/(R+1)
rtp <- (sum(z>=obs)+1)/(R+1)
if (identical(method,c("pearson", "kendall", "spearman"))) method <- "pearson"
if (identical(alternative,c("two.sided","less","greater"))) alternative <- "two.sided"
tc <- c("two.sided","less","greater")
pc <- c(2*min(ltp,rtp),ltp,rtp)
p <- signif(pc[tc==alternative],digits=3)
pv <- c((p>=0.001)&(p<=0.999),(p<0.001),(p>0.999))
pt <- c(p,"P < 0.001","P > 0.999")
p.value <- pt[pv]
ac <- c("relation","neg.relation","pos.relation")
alt <- ac[tc==alternative]
# FOR THE RESULTS
stat.name <- paste(method,".","cor",sep="")
results <- 
 list(Perm.values=z,Header=c("RESULTS OF PERMUTATION RELATIONSHIP TEST\n",
      paste("BASED ON",R,"REPLICATIONS")),Variable.1=var1.name,
      Variable.2=var2.name,n=n,Statistic=stat.name,Observed=obs,
      Null="no.relation",Alternative=alt,P.value=p.value,p.value=p)
class(results) <- "perm.two.var";  # permutation, two variables.
results
}
