library(magrittr)

context('opt_expand')

# Simple tests 

# No Options 
opts <- optigrab:::str_to_opts("")

t <- optigrab:::opt_expand(opts)
t  %>% expect_is("character")
t  %>% length  %>% expect_equal(1)


# Value
t <- optigrab:::opt_expand( optigrab:::str_to_opts( "one") )
expect_equal( t, "one" )


# Flag
t <- optigrab:::opt_expand( optigrab:::str_to_opts( "--flag" ))
expect_identical( t, "--flag" )


# Long-Flag
t <- optigrab:::opt_expand( optigrab:::str_to_opts( "--long-flag" ))
expect_identical( t, "--long-flag" )


# Short-flag 
t <- optigrab:::opt_expand( optigrab:::str_to_opts( "-f" ))
expect_identical( t, "-f" )


# Value Value 
t <- optigrab:::opt_expand( optigrab:::str_to_opts( "value1 value2" ) )
t %>% length %>% expect_equal(2) # expect_equal( length(t), 2 )
expect_identical( t, c( 'value1', 'value2' ) )


# Flag Value
t <- optigrab:::opt_expand( optigrab:::str_to_opts( "--flag value" ) ) 
expect_is( t, "character" )
expect_equal( length(t), 2 )
expect_equal( t[[1]], "--flag" )
expect_equal( t[[2]], "value" )

# Short-flag Value
t <- optigrab:::opt_expand( optigrab:::str_to_opts( "-f value" ) ) 
expect_is( t, "character" )
expect_equal( length(t), 2 )
expect_equal( t[[1]], "-f" )
expect_equal( t[[2]], "value" )




# Flag=Value
t <- optigrab:::opt_expand( optigrab:::str_to_opts( "--flag=value" ) ) 
expect_is( t, "character" )
expect_equal( length(t), 2 )
expect_equal( t[[1]], "--flag" )
expect_equal( t[[2]], "value" )





# Flag Flag 
t <- optigrab:::opt_expand( optigrab:::str_to_opts( "--flag1 --flag2" ) ) 
expect_is( t, "character" )
expect_equal( length(t), 2 )
expect_equal( t[[1]], "--flag1" ) 
expect_equal( t[[2]], "--flag2" )

# Flag Flag Value 

# Flag Value Valeu




# Complex 
opts <- optigrab:::str_to_opts( '/opt/r/R-2.13.0-default/lib/R/bin/exec/R --slave --no-restore  
         --file=./test.r --args --args --name fred --date 2011-05-17 -b=1 
         --end-date=2011-05-20 -a'
)


compare <- c(  "--args", "--name", "fred", "--date", "2011-05-17", "-b", "1" 
          , "--end-date", "2011-05-20", "-a" )


expanded <- optigrab:::opt_expand( opts=opts )



expect_that( is.character(expanded), is_true() )
expect_equal( length(expanded), 10 )
expect_identical( expanded, compare ) 
