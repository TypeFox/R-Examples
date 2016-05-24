negbin.rel.pca <-
function(x,mu.min,mu.max,len=200,plot=TRUE,seed=132)
{
if(!is.list(x)) {x <- lapply(seq_len(nrow(x)), function(i) x[i,])}
if(len > length(x)) stop("len exceeds row-length of data entered: For PCA number of columns cannot exceed number of rows.")

# Constructing the Distance matrix
mat.km <- matrix(NA,nrow = length(x),ncol = len)
rownames(mat.km) <- names(x)
colnames(mat.km) <- names(x)

for (i in 1:length(x)) 
{
y <- x[[i]]
mat.km[i,] <- wt.rel.negbin(y,mu.min,mu.max,plot = FALSE, len = len)$Val
colnames(mat.km) <- wt.rel.negbin(y,mu.min,mu.max,plot = FALSE, len = len)$mu
}

#---------------------------------------
# PCA
#---------------------------------------
lilkd.princomp <- princomp(mat.km)
if(plot==TRUE)biplot(lilkd.princomp,var.axes=FALSE,cex=0.7,col=c(1,0))
out <- list("PCA.Output" = lilkd.princomp)
return(out)
}
