perm.repl <-
function( data )
{
# Check that data has item, method and repl
rq.nam <- c("meth","item","y")
if( sum( !is.na( wh <- match( rq.nam, names( data ) ) ) ) < 3 ) stop(
"\nPR: The supplied dataframe misses column(s) named ", rq.nam[is.na(wh)], ".\n" )

# Reorder data arbitrarily within (method,item)
random <- runif( nrow( data ) )
new.order <- order( data$meth, data$item, random )
data <- data[new.order,]

# Make the replicates in the new ordering
Meth( make.repl( data ), print=FALSE )
}

make.repl <-
function( data )
{
# Check that data has item, method and repl
rq.nam <- c("meth","item")
if( sum( !is.na( wh <- match( rq.nam, names( data ) ) ) ) < 2 ) stop(
"\nThe supplied dataframe misses column(s) named ", rq.nam[is.na(wh)], ".\n" )

# Was it a Meth object?
was.Meth <- inherits( data, "Meth" )

# For this to work we must sort the dataframe first
ord <- order(data$meth,data$item)
data <- data[order(data$meth,data$item),]

# 0: (xx <- ) uniquely number all combinations of (method,item).
# 1: (tapply) find the smallest sequence number within each (method,item).
# 2: subtract this from the sequence number, to get 0,1,... within each (m,i).
# 3: add 1 to get new replicate numbers.
# 4: (as.vector) turn it into a vector.
# Note the use of factor(interaction...) to make sure that empty levels
#      of method or item don't screw up.
xx <- as.integer( factor( interaction( data$meth, data$item ) ) )
data$repl <- as.vector( 1:nrow(data) -
                        tapply( 1:nrow(data), xx, min )[xx] + 1 )
# Then reestablish the sequence:
data <- data[order(ord),]
if( was.Meth ) Meth( data, print=FALSE ) else data
}

has.repl <-
function( data )
{
# This functions tests whether there actually are replicates
# in the dataframe x
if( !inherits( data, "data.frame" ) )
  stop( "argument must be a data frame (preferably a Meth object)")
if( !inherits( data, "Meth" ) ) data <- Meth( data, print=FALSE )
any( table( data$meth, data$item ) > 1 )
}
