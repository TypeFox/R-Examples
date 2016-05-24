
test_that( "but operators work", {
  expect_equal( rnorm %but% glm , rnorm )

	grep. <- grep %but% list( perl = FALSE, pattern = "foobar" )
	f <- formals( grep. )
	expect_equal( f$pattern , "foobar" )
	expect_equal( f$perl , FALSE)
	
	grep. <- grep %but% list( foo = "bar" )
	f <- formals( grep. )
	expect_true( "foo" %!in% names(f) )

	expect_true( formals( grep %but% "+p" ) $perl )
	expect_true( formals( grep %but% "!p" ) $perl )
	expect_true( ! formals( grep %but% "-p" ) $perl )
	
	expect_true( formals( grep %but% "pf" )$fixed )
	expect_false( formals( grep %but% "p-f" )$fixed )
	expect_true( formals( grep %but% "p!f" )$fixed )
	
	expect_warning( grep %but% "a" )
	
})

