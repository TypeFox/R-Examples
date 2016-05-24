sample.Meth <-
function( x,
        how = "random",
          N = if( how=="items" ) nlevels( x$item ) else nrow(x),
        ... )
{
if( !inherits( x, "Meth" ) ) x <- Meth( x )

if( tolower(substr(how,1,1))=="r" )
  {
  Nr <- nrow( x )
  new.x <- x[sample(1:Nr,Nr,replace=T),]
  new.x <- Meth( make.repl( new.x ), print=FALSE )
  }
  
if( tolower(substr(how,1,1))=="l" )
  {
  Nm <- nlevels( x$meth )
  IxR <- interaction( x$item, x$repl )
  ir.id <- sample( levels( IxR ), N/Nm, replace=T )
  new.x <- x[NULL,]
  for( i in 1:length(ir.id) )
     {
     new.x <- rbind( new.x, x[IxR==ir.id[i],] )
     }
  new.x <- Meth( make.repl( new.x ), print=FALSE )
  }
  
if( tolower(substr(how,1,1))=="i" )
  {
  i.id <- sample( levels( x$item ), N, replace=T )
  new.x <- cbind(x,new.item=0)[NULL,]
  for( i in 1:length(i.id) )
     {
     new.x <- rbind( new.x,
                     cbind( x[x$item==i.id[i],],
                            new.item=i ) )
     }
  new.x$old.item <- new.x$item
  new.x$item     <- new.x$new.item
  new.x <- new.x[,-grep("new.item",names(new.x))]
  new.x <- Meth( make.repl( new.x ), print=FALSE )
  }

new.x
}
