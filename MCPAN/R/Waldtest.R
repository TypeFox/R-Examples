`Waldtest` <-
function(estp, varp, cmat, alternative="greater", dist="MVN")
{

 k <- ncol(cmat)
 m <- nrow(cmat) 

if(any(c(length(estp), length(varp))!=k))
 {stop("estp, varp and ncol(cmat) must be of equal length")}
 
dist<-match.arg(dist, choices=c("MVN", "N"))

# point estimate:
 estC <- cmat %*% estp

# Variance estimates for the contrasts using unpooled (groupwise) variance estimates:
 varC <- (cmat^2) %*% (varp)  

# Teststatistic:

 teststatv <- estC/sqrt(varC)
 
switch(dist,
"MVN"={

CorrMat<-corrMatgen(CM=cmat, varp=varp)

p.val.adj <- numeric(length=m) 

# adjusted p.values for the single hypotheses

for( i in 1:m)
 {

if(alternative=="two.sided")
 {teststati <- abs(teststatv[i])
 p.val.adj[i] <- 1 - pmvnorm(lower=rep(-teststati,m), upper=rep(teststati,m), mean=rep(0,m), corr=CorrMat)}

if(alternative=="greater")
 {teststati <- teststatv[i]
  p.val.adj[i] <- 1 - pmvnorm(lower=rep(-Inf,m), upper=rep(teststati,m), mean=rep(0,m), corr=CorrMat)}

if(alternative=="less")
 {teststati <- teststatv[i]
  p.val.adj[i] <- 1 - pmvnorm(lower=rep(teststati,m), upper=rep(Inf,m), mean=rep(0,m), corr=CorrMat)}
 
 }

pmaxtest<-min(p.val.adj)

},
"N"={
p.val.adj <- numeric(length=m) 

# adjusted p.values for the single hypotheses

for( i in 1:m)
 {

if(alternative=="two.sided")
 {teststati <- abs(teststatv[i])
 p.val.adj[i] <- min(1, 2*(1 - pnorm(q=teststati)))}

if(alternative=="greater")
 {teststati <- teststatv[i]
  p.val.adj[i] <- 1 - pnorm(q=teststati)}

if(alternative=="less")
 {teststati <- teststatv[i]
  p.val.adj[i] <- pnorm(q=teststati)}
 
 }

pmaxtest<-NA

})

return(
list(teststat=teststatv, pval=pmaxtest, p.val.adj=p.val.adj, alternative=alternative, dist=dist)
)
}

