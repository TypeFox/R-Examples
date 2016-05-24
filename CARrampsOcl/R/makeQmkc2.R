function(nr, nc) {
# makes CAR(1) structure matrix for nr by nc lattice data

N <- nr * nc

if( nr * nc  == 1)
   Q <- matrix( 1, nrow=1, ncol=1)
else
{
Q <- matrix( 0, nrow=N, ncol=N )
if(nr > 1)
for( i in 1:nc )
   for( j in 1:(nr-1) )
   {
      ind <- (i-1) * nr + j
      Q[ ind, ind+1  ] <- Q[ ind+1, ind ] <- -1
   }

if(nc > 1)
for( i in 1:(nc-1) )
   for( j in 1:nr )
   {
      ind <- (i-1) * nr + j
      Q[ ind, ind+nr  ] <- Q[ ind+nr, ind ] <- -1
   }

diag(Q) <- - apply( Q, 1, sum)
}

Q

}
