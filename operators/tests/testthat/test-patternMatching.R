
test_that( "pattern matching operators work", {
  blue <- rep( "blue" , 10 )
  res  <- blue %~% "^b" 
  expect_true(all(res))
  expect_equal(length(res) , 10 )
  
  cols <- c("blue", "red" )
  expect_equal( cols %~% "^b", c(TRUE, FALSE))
  expect_equal( cols %!~% "^b", c(FALSE, TRUE))
  
  expect_equal( cols %~|% "^b" , "blue")
  expect_equal( length(cols %~|% "f") , 0  )
  
  expect_equal( cols %!~|% "^b" , "red" )
  expect_equal( cols %!~|% "f" , cols )
  expect_true( cols %~*% "e" )
  expect_false( cols %~*% "r")
  expect_false( cols %~*% "i")
  expect_false( cols %!~*% "e")
  expect_true( cols %!~*% "r" )
  expect_true( cols %!~*% "i" )
  
  expect_true( cols %~+% "e" )
  expect_true( cols %~+% "r" )
  expect_false( cols %~+% "i" )
  expect_false( cols %!~+% "e" )
  expect_false( cols %!~+% "r" )
  expect_true( cols %!~+% "i" )
  
})

test_that( "pattern removing works", {

  cols <- c( "blue", "red" )
  expect_equal( cols %-~% "e", c("blu", "rd") )
  expect_equal( cols %-~% "pp", cols )
  expect_equal( cols %-~% "blue", c( "", "red") )
  expect_equal( cols %-~|% "b", "lue" )
  expect_equal( cols %-~|% "e", c("blu", "rd" ))

  cols <- c( "blue.col", "pink", "red.stuff" )
  expect_equal( as.vector(cols %o~|% "\\..*$"  ), c(".col", ".stuff"))
  expect_equal( as.vector(cols %o~|% "\\.(.*)$"), c("col", "stuff"))

})
