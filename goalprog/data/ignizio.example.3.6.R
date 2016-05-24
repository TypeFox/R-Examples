###
### ignizio.example.3.6
###
### The data in this file define the coefficients, target values and achievements
### for Example 3-6 in Ignizio (1976) Goal Programming and Extensions.
###
coefficients <- matrix( 
  c( 1.00, 1.00, 1.00, 1.00,
     1.00, 0.00, 0.00, 0.00,
     0.00, 1.00, 0.00, 0.00,
     0.00, 1.00, 0.00, 0.00,
     1.00, 0.00, 1.00, 0.00,
     0.00, 0.00, 0.00, 1.00,
     0.06, 0.05, 0.08, 0.07 ), 
  nrow=7, byrow=TRUE )
targets <- c( 50000, 20000, 5000, 15000, 10000, 30000, 4000 )
achievements <- data.frame( matrix( 
  c( 1, 1, 1.0, 0.0, 
     2, 2, 0.0, 1.0, 
     3, 2, 0.0, 2.0, 
     4, 2, 2.0, 0.0, 
     6, 3, 0.0, 1.0, 
     5, 4, 1.0, 0.0, 
     7, 4, 0.0, 1.0),
  nrow=7, byrow=TRUE ) )
names( achievements ) <- c( "objective", "priority", "p", "n" )
