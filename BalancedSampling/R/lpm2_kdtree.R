lpm2_kdtree <- function(
  prob,
  x,                      # data to mode seek on
  m=40
) {
  
  
  n <- NROW(x)
  K <- length(x) / n
  m <- min( m, n)
  
  
  # send our data to the C program
  r.result <- .C("R_lpm3",
                 as.double( t(x) ),                 # 1 data we query
                 as.double( prob ),                 # 2 probability
                 as.integer( n ),              
                 as.integer( K ),                   
                 as.integer( m )                    # max leaves per node
  )
  
  return( (1:n)[ r.result[[2]] > .5 ] )
  
}