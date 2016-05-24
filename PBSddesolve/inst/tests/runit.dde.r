
test.simple <- function()
{
	#simple horizontal line
	grad <- function( t, y )
	{
		return( 0 );
	}
	
	#test an int
	yinit <- 1
	out <- dde( y=yinit, func=grad, times=seq( 0, 1,length=100 ), hbsize=0)
	checkTrue( all( out$y1 == yinit ) )
	
	#test a float
	yinit <- 10.5
	out <- dde( y=yinit, func=grad, times=seq( 0, 1,length=100 ), hbsize=0)
	checkTrue( all( out$y1 == yinit ) )

	#test negative float
	yinit <- -10.5
	out <- dde( y=yinit, func=grad, times=seq( 0, 1,length=100 ), hbsize=0)
	checkTrue( all( out$y1 == yinit ) )

}

#check yinit names are properly passed to grad function
test.y_init_names <- function()
{
	grad <- function( t, y )
	{
		checkTrue( all( n == names( y ) ) )
		return( c(0, 1) );
	}
	
	yinit <- c( y1 = 1, y2 = 2 )
	n <- names( yinit )
	out <- dde( y=yinit, func=grad, times=seq( 0, 1,length=100 ), hbsize=0)
}
