
test_that( "repeat operators work", {
    
  expect_equal( "--" %x=% 10, paste(rep("-", 20), collapse = "") )
  expect_equal( "-+" %x=% 10, paste(rep("-+", 10), collapse = ""))
	
	expect_equal( nchar("--" %x=|% 50), 50  )
	expect_equal( nchar("|-<+>-|" %x=|% 50), 50)
	expect_equal( "|-<+>-|" %x=|% 50, substring(paste(rep("|-<+>-|", 10), collapse = "" ), 0, 50) )
  
})

