sim.pois <-
function(lambda,n.obs,n.val)
{
#---------------------
# Error Check
#---------------------
n.integer <- function(x, tol = .Machine$double.eps) 
{ abs(x - round(x)) < tol }
if( mean(n.integer(n.obs)) !=1 || n.obs<= 0) stop("Num of Obs must be positive integer")
if( mean(n.integer(n.val)) !=1 || n.val<= 0) stop("Length of Each Obs must be positive integer")
if(any(lambda < 0)) stop("Poisson mean is negative")

dat.sim <- vector("list",length(lambda))
for (i in 1:length(lambda))
{
set.seed(132);set.seed(sample(.Random.seed,1))
X <- matrix(NA,nrow = n.obs, ncol = n.val)
for (j in 1:n.obs)
{
X[j,] <- rpois(n.val,lambda[i])
}
dat.sim[[i]] <- X
}

dat.sim <- do.call(rbind,dat.sim)
rownames(dat.sim) <- 1:dim(dat.sim)[1]
colnames(dat.sim) <- 1:dim(dat.sim)[2]
dat.sim <- dat.sim[sample(rownames(dat.sim),replace=FALSE),]
return(dat.sim)
}
