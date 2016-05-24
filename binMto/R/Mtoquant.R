"Mtoquant" <-
function(nc, nx, pc, px, conf.level=0.95, adj="Dunnett", alternative="two.sided")
{

# # calculate the correlation matrix R:
# # see Piegorsch, Biometrics (1991) 47(1), 45-52:
# # (3.3)
#
# n0 = sample size control
# nx = sample sizes of the treatment
# p0 = observed proportion in control
# px = observed proportions of treatments

corrmat <- function(n0, nx, p0, px)
{

k=length(nx)
# here, k=the number of treatment groups
biv <- numeric(length=k)
 for(i in 1:k)
  {
  biv[i] <- 1/sqrt(1 + n0/nx[i] * ( ( px[i]*(1-px[i]) ) / ( p0*(1-p0) ) )  ) 
  }

corr <- diag(x=1, nrow=k)

# the columns:

for(e in 1:k)
{

 for(a in 1:k)
  {   
  if( e !=a ) {corr[e,a] <- biv[e]*biv[a]}
  }
}

# check for NA (if plug in estimates are used p0=0, or px=1) 
# NA are replaces by 0

for(i in 1:length(corr))
 {
  if(is.na(corr[i])){corr[i]<-0}
    }
return(corr)

}
# end of corrmat



corrmatappr <- function(n0, nx, p0, px)
{

k=length(nx)
biv <- numeric(length=k)
 for(i in 1:k)
  {
  biv[i] <- 1/sqrt( 1 + n0/nx[i] ) 
  }

corr <- diag(x=1, nrow=k)

# the columns:

for(e in 1:k)
{

 for(a in 1:k)
  {   
  if( e !=a ) {corr[e,a] <- biv[e]*biv[a]}
  }
}

return(corr)
}
# end of corrmatappr

k <- length(nx)
# k= number of treatment groups

if(adj=="Dunnett" && alternative=="two.sided")
 {
  corrm <- corrmat(n0=nc, nx=nx, p0=pc, px=px)
  quant <- qmvnorm(p=conf.level,corr=corrm, tail="both.tails")$quantile
 }

if(adj=="Dunnett" && alternative=="less")
 {
  corrm <- corrmat(n0=nc, nx=nx, p0=pc, px=px)
  quant <- qmvnorm(p=conf.level,corr=corrm, tail="lower.tail")$quantile
 }


if(adj=="Dunnett" && alternative=="greater")
 {
  corrm <- corrmat(n0=nc, nx=nx, p0=pc, px=px)
  quant <- qmvnorm(p=conf.level,corr=corrm, tail="upper.tail")$quantile
 }


# Bonf-adjustment

if(adj=="Bonf")
 { 
 if(alternative=="two.sided")
  {
   quant <- qnorm(p=1-(1-conf.level)/(2*k))
  }

 if(alternative=="less")
  {
   quant <- qnorm(p=1-(1-conf.level)/(k), lower.tail=TRUE) 
  }

 if(alternative=="greater")
  {
   quant <- qnorm(p=1-(1-conf.level)/(k), lower.tail=FALSE) 
  }
 }

if(adj=="Unadj")
 { 
 if(alternative=="two.sided")
  {
   quant <- qnorm(p=1-(1-conf.level)/(2))
  }

 if(alternative=="less")
  {
   quant <- qnorm(p=1-(1-conf.level), lower.tail=TRUE) 
  }

 if(alternative=="greater")
  {
   quant <- qnorm(p=1-(1-conf.level), lower.tail=FALSE) 
  }
 }



return(quant)
}

