library(hash)
library(testthat)

h0 <- hash()
h <- hash( letters, 1:26 )

# OBJECT CREATION / TYPE
  for( h in list( h0, h) ) expect_is( h, "hash" )


# EMPTY HASH
  expect_that( length(h0), equals(0) )
  expect_that( length(h0), equals(0) )

# POPULATED HASH  
  expect_that( length(h), equals(26) )
  expect_that( keys(h), is_identical_to(letters) )
               
               
# ALL HASHES
  for( h in list( h0, h) ) {
    expect_that( 
      h[['missing']], is_identical_to(NULL)
      , label="Attempt to retrieve missing key" 
    )  
  }
             
             
# TEST [[             
  for( h in list( h0, h ) ) {
    expect_error( h[[NULL]] )
    expect_error( h[[NA] ])
    expect_error( h[[]] )
    expect_error( h[[letters]] )
  }

  for( n in 1:26 ) expect_that( h[[ letters[n] ]], equals(n) )  


               
# TEST $                
  expect_that( h$a, equals(1)  )
  expect_that( h$z, equals(26) )             
        
               
# TEST [
  for( h in list( h0, h ) ) {
    expect_error( h[[NULL]] )
    expect_error( h[[NA] ])
    expect_error( h[[]] )
    expect_error( h[[letters]] )
  }
  expect_that( h[ letters ], equals(h) )             
  
               
# TEST [[ <- 
  h[['a']] <- -1             
  expect_that( h[['a']], equals(-1) )             
               
# TEST $<-             
  h$b <- -2
  expect_that( h$b, equals(-2) )

# TEST [ <- 
