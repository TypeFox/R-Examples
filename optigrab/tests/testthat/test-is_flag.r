# library(testthat)
# library(optigrab)
library(magrittr)

# No options
context( 'is.flag' )
message('ISFLAG')


# No option - empty string
opts <- optigrab:::str_to_opts('')

optigrab:::is.flag(opts)  %>%  expect_false

# Value 
opts <- optigrab:::str_to_opts( "value")
expect_false( optigrab:::is.flag( opts ) )

# Flag (BOOLEAN )
opts <- optigrab:::str_to_opts( "--one" ) 
expect_true( optigrab:::is.flag( opts) )
 
# Flag, short (BOOLEAN)
opts <- optigrab:::str_to_opts( "-o" )
expect_true( optigrab:::is.flag( opts) )

# Flag value 
opts <- optigrab:::str_to_opts( "--one 1")
expect_identical( optigrab:::is.flag(opts), c( T, F))

# Flag flag value
opts <- optigrab:::str_to_opts( "--one -o 1 ")
expect_identical( optigrab:::is.flag(opts), c( T, T, F))

# Flag value value 
opts <- optigrab:::str_to_opts( "--one 1 2 ")
expect_identical( optigrab:::is.flag(opts), c( T, F, F))
