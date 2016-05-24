context( "opt_this_file" )
message( "THIS_FILE" )
# TYPICAL  
  flags <- str_to_opts( "Rscript --slave --no-restore --file=test-this_file.r --args sub1" )
  file  <- this_file(flags, full.path = FALSE ) 
  
  expect_equal( file, "test-this_file.r" )
  
# MULTIPLE --file
  flags <- str_to_opts( "Rscript --slave --no-restore --file=test-this_file.r --args sub1 --file other-file" )
  this_file(flags, full.path = FALSE)

  expect_equal( file, "test-this_file.r" )
