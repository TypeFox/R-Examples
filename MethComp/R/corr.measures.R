# A function that returns the values of some of the crap
# association measures proposed in the literature
corr.measures <-
function( x, y )
{
compl <- complete.cases( data.frame( x, y ) )
  x <- x[compl]
  y <- y[compl]
  
 Vx <- var(x)
 Mx <- mean(x)
 Vy <- var(y)
 My <- mean(y)
 MD <- mean(x-y)
 VD <- var(x-y)
Cxy <- cov(x,y)
MSD <- MD^2 + VD
corr <- cor( x, y )
CCC <- 2*Cxy / (MSD+2*Cxy)
AcC <- 2 / ( sqrt(Vx/Vy) + sqrt(Vy/Vx) + (Mx-My)^2/sqrt(Vx*Vy) )
res <- c( corr, MSD, CCC, AcC )
names( res ) <- c( "Corr", "MSD", "CCC", "Acc.C" )
res
}

# Functions to point at the middle or the extremes
middle <-
function( w, rm=1/3 )
{
qnt <- quantile( w, probs=0:1 + c(1,-1)*rm/2 )
w > qnt[1] & w < qnt[2]
}

ends <-
function( w, rm=1/3 )
{
!middle( w, 1-rm )
}
