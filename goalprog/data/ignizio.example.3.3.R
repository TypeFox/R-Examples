###
### ignizio.example.3.3
###
### The data in this file define the coefficients, target values and achievements
### for Example 3-3 in Ignizio (1976) Goal Programming and Extensions.
###
coefficients <- matrix( c( 1, 1, 1, 0, 5, 3, 1, 1 ), nrow=4, byrow=TRUE )
targets <- c( 10, 4, 56, 12 )
achievements <- data.frame( matrix( 
  c( 1, 1, 2, 0, 
     2, 1, 3, 0, 
     3, 2, 0, 1, 
     4, 3, 1, 0), nrow=4, byrow=TRUE ) )
names( achievements ) <- c( "objective", "priority", "p", "n" )
