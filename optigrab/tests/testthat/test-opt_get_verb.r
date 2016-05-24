library(magrittr)

context( "opt_get_verb" )

# TYPICAL 
  flags <- 
    "Rscript --slave --no-restore --file=my-file --args cmd1" %>%
    str_to_opts 

  flags %>% opt_get_verb %>% expect_equal("cmd1")

  
# MULTIPLE --file
  flags <- 
    optigrab:::str_to_opts( 
      "Rscript --slave --no-restore --file=my-file --args cmd1 --date now" 
    ) 

  flags %>% opt_get_verb %>% expect_equal("cmd1")



# NOT IN FIRST POSITION
  flags  <- optigrab:::str_to_opts( "Rscript --slave --no-restore --file=my-file --args --date now cmd1" )
  date   <- opt_get( name="date", n=1, opts=flags )  
  cmd    <- opt_get_verb( flags )
  
  expect_equal( cmd, "cmd1" )

  
# NOT IN FIRST POSITION OBFUSCATED BY PARAMETER THAT TAKES NO ARGUMENTS 
  flags <- 
    "Rscript --slave --no-restore --file=my-file --args --verbose cmd1" %>%
    str_to_opts
  
  verbose <- opt_get( name="verbose", n=0, opts=flags )
  
  flags %>% opt_get_verb %>% expect_equal("cmd1")

