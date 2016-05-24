# To avoid NOTE in 'R CMD check' for an unknown variable, specifically,
# the note "no visible binding for global variable 'Freq'".
globalVariables("Freq")

perm.ind.test <-
function(x, type = c("cont", "flat", "raw"), var.names = NULL, R = 9999)
{
# Note: unfct and unct are due to Marc Schwartz.
# unfct converts a "flat" contingency table to a raw data frame. 
unfct <- function(x) 
{
z <- sapply(1:nrow(x), function(i) x[rep(i, each = x$Freq[i]), ],simplify = FALSE)
z <- subset(do.call("rbind", z), select = -Freq)
for (i in 1:ncol(z)) {z[[i]] <- type.convert(as.character(z[[i]]))}
data.frame(z,row.names=NULL)
}
# unct converts a contingency table (of class table) to a raw data frame.
unct <- function(x) 
{
y <- as.data.frame(x)
unfct(y)
}
if (identical(type,c("cont", "flat", "raw"))) type <- "cont"
if (type=="cont") x <- unct(as.table(as.matrix(x[-1])))
else
 (if (type=="flat") {names(x)[3]="Freq"; x <- unfct(x)})
if (is.null(var.names)) var.names <- c("Var.1", "Var.2")
obs <- suppressWarnings(chisq.test(x[[1]],x[[2]],correct=FALSE)$statistic)
n <- length(x[[1]])
z <- vector(length=R)
for (i in 1:R)
{
 u <- sample(x[[2]])
 z[i] <- suppressWarnings(chisq.test(x[[1]],u,correct=FALSE)$statistic)
} 
p <- signif((sum(z >= obs)+1)/(R+1),digits=3)
pv <- c((p>=0.001)&(p<=0.999),(p<0.001),(p>0.999))
pt <- c(p,"P < 0.001","P > 0.999")
p.value <- pt[pv]
# FOR THE RESULTS
stat.name <- "chi.square"
results <- 
 list(Perm.values=z,Header=c("RESULTS OF PERMUTATION INDEPENDENCE TEST\n",
      paste("BASED ON",R,"REPLICATIONS")),Variable.1=var.names[1],
      Variable.2=var.names[2],Statistic=stat.name,Observed=obs,n=n,
      Null="nonassociated",Alternative="associated",P.value=p.value,
      p.value=p)
class(results) <- "perm.two.var"  # permutation, two variables.
results
}
