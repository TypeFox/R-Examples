pois.rel.clust <-
function(x,lambda.min,lambda.max,nclust = 3,len=200,plot=TRUE,seed=132)
{
if(!is.list(x)) {x <- lapply(seq_len(nrow(x)), function(i) x[i,])}

# Constructing the Distance matrix
mat.km <- matrix(NA,nrow = length(x),ncol = len)
rownames(mat.km) <- names(x)
colnames(mat.km) <- names(x)

for (i in 1:length(x)) 
{
y <- x[[i]]
mat.km[i,] <- wt.rel.pois(y,lambda.min,lambda.max,plot=FALSE,len=len)$Val
colnames(mat.km) <- wt.rel.pois(y,lambda.min,lambda.max,plot=FALSE,len=len)$Lambda
}

#---------------------------------------
# K means Clustering 
#---------------------------------------
set.seed(seed)
fit.kmeans <- kmeans(mat.km,nclust)
groups <- fit.kmeans$cluster
sort(groups);table(groups)

if(plot==TRUE)matplot(as.numeric(colnames(mat.km)),t(mat.km),type="l",col = groups,lty = groups,ylab = "Wt. Rel. Likld",xlab = expression(lambda))

out <- list("Wt. Rel. Likld" = mat.km,"Cluster.Assignment"= groups,"table" = table(cluster = groups))
return(out)
}
