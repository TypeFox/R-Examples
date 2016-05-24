###
### ignizio.example.3.2
###
### The data in this file define the coefficients, target values and achievements
### for Example 3-2 in Ignizio (1976) Goal Programming and Extensions.
###
coefficients <- matrix( 
  c( 1, 1,
     1, 1, 
     1, 0, 
     0, 1 ), 
  nrow=4, byrow=TRUE )
targets <- c( 40, 100, 30, 15 )
achievements <- data.frame( matrix( 
  c( 3, 1, 1.0, 0.0, 
     4, 1, 1.0, 0.0, 
     1, 2, 1.0, 0.0, 
     2, 3, 0.0, 1.0, 
     3, 4, 0.0, 1.0, 
     4, 4, 0.0, 1.5),
  nrow=6, byrow=TRUE ) )
names( achievements ) <- c( "objective", "priority", "p", "n" )
