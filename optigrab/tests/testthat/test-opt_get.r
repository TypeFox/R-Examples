library(magrittr)


context( "opt_get")
message("OPT_GET" )

flags <- c( "f", "flag", "long-flag" )

# No options # 
#  Opt Stings that contain neither values or flags

for( f in flags ) {
  t <- opt_get( f, opts=optigrab:::str_to_opts('') )
  expect_identical( t, NA )
}

for( f in flags ) {
  t <- opt_get( f, opts=optigrab:::str_to_opts("") )
  expect_identical( t, NA )
}


# Value(s) #
#  Opt Strings that do not contain any flags
#  These should not produce any values since there are no flags.
opt_strings <- c( 'v', 'value', 'val1 val2', 'val1 val2 val3' )
opts = optigrab:::str_to_opts(opt_strings)

for( str in opt_strings ) {
  for ( f in flags ) {
    t <- opt_get( f, opts=opts )
    expect_identical( t, NA )
  } 
}

# Flag (BOOLEAN) #
context("boolean flags")
opt_string <- c( '-f --flag --long-flag' )
opts <- optigrab:::str_to_opts(opt_string)
flags <- gsub( "^-+", "", opts )


for ( flag in flags  ) {
  expect_true( opt_get( name=flag, n=0, opts=opts ) )
}

# Flag Value 

for( str in opt_strings ) {
  t <- opt_get( flags, opts=optigrab:::str_to_opts( "-f value") )
  expect_equal( t, "value" )
}



# # HERE IS A TYPICAL RScirpt command
# opts <- '/opt/r/R-2.13.0-default/lib/R/bin/exec/R --slave --no-restore  
#          --file=./test.r --args --args --name fred --date 2011-05-17 -b=1 
#          --end-date=2011-05-20 -a'
# 
# opts <- optigrab:::str_to_opts( opts )
# 
# # [1] "/opt/r/R-2.13.0-default/lib/R/bin/exec/R" "--slave"                                 
# # [3] "--no-restore"                             "--file=./test.r"                         
# # [5] "--args"                                   "--args"                                  
# # [7] "--name"                                   "fred"                                    
# # [9] "--date"                                   "2011-05-17"                              
# # [11] "-b=1"                                     "--end-date=2011-05-20"                   
# # [13] "-a"                                      
# 
# # TEST: Simple
# opt_grab("--args", n=0, opts=opts) %>% expect_true
# # opt_grab("--args", n=1, opts=opts) # FAIL
# opt_grab("--name", opts=opts) %>% expect_identical('fred')
# opt_grab("--date", opts=opts) %>% expect_equal('2011-05-17')
# opt_grab("-b", opts=opts)     %>% expect_equal("1") 
# opt_grab("--end-date", opts=opts) %>% equals('2011-05-20')
# expect_error( opt_grab("-a", opts=opts) )
# 
# 
# # expect_that( opt_grab("--name", opts=opts)     , is_identical_to('fred') )
# # expect_that( opt_grab("--date", opts=opts)     , equals('2011-05-17') )
# # expect_that( opt_grab("-b", opts=opts)         , equals("1") )
# # expect_that( opt_grab("--end-date", opts=opts) , equals('2011-05-20') )
# # expect_that( opt_grab("-a", opts=opts)          , throws_error() )
# 
# 
# # TEST: MISSING VALUES
# opt_grab("--missing", opts=opts) %>% expect_equal(NA) 
# # opt_grab("--missing", opts=opts) %>% expect_equal(NULL) 
# # opt_grab("--missing", default=NA, opts=opts)   %>% expect_equal(NA) 
# # opt_grab("--missing", default=0, opts=opts)    %>% expect_equal(0) 
# 
# # TEST: logical
# 
# expect_error( opt_grab( "-a", opts=opts ) )      # Error: No argument supplied. 
# 
# opt_grab( "-a", n=0, opts=opts ) %>% expect_true
# expect_error( opt_grab( "-a", n=1, opts=opts ) )
# 
# opt_grab( "-b", opts=opts ) %>%  expect_equal("1")        
# opt_grab( "-b", n=0, opts=opts ) %>% expect_true
# # expect_that( opt_grab( "-b", n=1, opts=opts, coerce=as.logical.character ), is_true() ) #FAIL.
# opt_grab( "-b", n=1, opts=opts ) %>%  equals("1" ) 
