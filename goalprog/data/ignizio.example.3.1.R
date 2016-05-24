###
### ignizio.example.3.1
###
### The data in this file define the coefficients, target values and achievements
### for Example 3-1 in Ignizio (1976) Goal Programming and Extensions.
###
coefficients <- matrix( 
  c(  10,  15,
     100, 100, 
       0,   1 ), 
  nrow=3, byrow=TRUE )
targets <- c( 40, 1000, 7 )
achievements <- data.frame( matrix( 
  c( 1, 1, 1, 0, 
     2, 2, 0, 1, 
     3, 3, 0, 1), 
  nrow=3, byrow=TRUE ) )
names( achievements ) <- c( "objective", "priority", "p", "n" )
