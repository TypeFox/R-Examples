abconv <-
function( a1, b1=1:4, a2=NULL, b2=NULL,
          col.names=c("alpha.2.1","beta.2.1","id.2.1") )
{
if( ( inherits( a1, "data.frame" ) |
      inherits( a1, "matrix" ) )
    & length( b1 )==4 )
  {
  cols <- a1
  wh <- b1
  a1 <- cols[,wh[1]]
  a2 <- cols[,wh[2]]
  b1 <- cols[,wh[3]]
  b2 <- cols[,wh[4]]
  }
a2.1 <- a2 - a1 * b2 / b1
b2.1 <- b2 / b1
id2.1 <- a2.1 / ( 1-b2.1 )
dfr <- data.frame( a2.1, b2.1, id2.1 )
names( dfr ) <- col.names
dfr
}

