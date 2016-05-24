context('make_names')


# -----
str <- c("foo and bar", "foo-and-bar")

expect_equal( 
    make_names( str, unique = FALSE )
  , c( "foo_and_bar", "foo_and_bar" ) 
)

expect_equal( 
    make_names( str, unique = TRUE )
    , c( "foo_and_bar", "foo_and_bar_1" ) 
)


# -----
str <- c("foo and bar", "foo.and_bar")
expect_equal(  
    make_names( str, unique = FALSE ) 
  , c("foo_and_bar", "foo_and_bar" )
)  


expect_equal(  
    make_names( str, unique = TRUE ) 
  , c("foo_and_bar", "foo_and_bar_1" )
)  


# -----
str <- c(".foo", "_bar")

expect_equal( 
    make_names( str )  
  , c( "foo", "bar" )
)

expect_equal( 
    make_names( str, leading="." )  
  , c( ".foo", ".bar" )
)

