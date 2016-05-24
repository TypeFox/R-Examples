
rearrange <- function( RS, RB, CS, CB, r ) { 

# r = rank of bootstrap matrix, so if < sample rank 
#     will ignore last sample axis
# RS = sample axes for row points (as columns)
# RB = bootstrap axes for row points (as columns)
# CS = sample axes for column points (as columns)
# CB = bootstrap axes for column points (as columns)
# T = matrix to rearrange xB so is equivalent to xS, i.e.
# xS ~= xB * T       
# find rearrangement of columns of RB and CB to maximise
# match = tr( abs(RS'*RB) + abs(CS'*CB) )
# Literal, uses all axes of bootstrap matrix up to its rank,
#  up to a maximum of maxrearrange=6

if (r>=1) {
maxrearrange <- 6
numrearranged <- min(r,maxrearrange)
switch(numrearranged,
  per <- matrix(1,1,1), 
  per <- rbind( c(1,2), c(2,1) ),
  per <- cbind( rep(1:3,each=2), c(2,3,1,3,1,2), c(3,2,3,1,2,1) ),
  { p <- cbind( rep(1:3,each=2), c(2,3,1,3,1,2), c(3,2,3,1,2,1), 4 )
    per <- rbind(p,p+3*(p==1)-3*(p==4),p+2*(p==2)-2*(p==4),p+(p==3)-(p==4) ) },
  { p <- cbind( rep(1:3,each=2), c(2,3,1,3,1,2), c(3,2,3,1,2,1), 4 )
    p <- rbind(p,p+3*(p==1)-3*(p==4),p+2*(p==2)-2*(p==4),p+(p==3)-(p==4) )
    p <- cbind(p,5)
    per <- rbind(p,p+4*(p==1)-4*(p==5),p+3*(p==2)-3*(p==5),p+2*(p==3)-2*(p==5),p+(p==4)-(p==5) ) },
  { p <- cbind( rep(1:3,each=2), c(2,3,1,3,1,2), c(3,2,3,1,2,1), 4 )
    p <- rbind(p,p+3*(p==1)-3*(p==4),p+2*(p==2)-2*(p==4),p+(p==3)-(p==4) )
    p <- cbind(p,5)
    p <- rbind(p,p+4*(p==1)-4*(p==5),p+3*(p==2)-3*(p==5),p+2*(p==3)-2*(p==5),p+(p==4)-(p==5) ) 
    p <- cbind(p,6)
    per <- rbind(p,p+5*(p==1)-5*(p==6),p+4*(p==2)-4*(p==6),p+3*(p==3)-3*(p==6),p+2*(p==4)-2*(p==6),p+(p==5)-(p==6) ) }
)

nper <- dim(per)[1]
match <- matrix(0,nper,1)
for (i in 1:nper) { 
match[i] = sum( diag( abs( t(RS[ ,1:numrearranged]) %*% RB[ ,per[i, ] ] + t(CS[ ,1:numrearranged]) %*% CB[ ,per[i, ] ] ) ) )
}

posn <- which.max(match)
same <- posn==1
I <- diag( rep(1,numrearranged) )
T <- I[ ,per[posn, ] ]
t <- diag( t(RS[ ,1:numrearranged]) %*% RB[ ,per[posn, ] ] + t(CS[ ,1:numrearranged]) %*% CB[ ,per[posn, ] ] )
T <- T %*% diag( (t>=0)-(t<0), nrow=numrearranged, ncol=numrearranged )
} else {  # r=0, i.e. raw bootstrap matrix is rank 1
T <- matrix(1,1,1)
numrearranged <- 1
match <- 0
same <- 0
}

list(T=T,numrearranged=numrearranged,match=match,same=same)

}


