
test.options <- function()
{
	setwd( tempdir() )
	
	fname <- "my_pkg.txt"
	unlink( fname )
	checkTrue( file.exists( fname ) == FALSE )
	
	#check initial value
	.mypkg <- new( "option", filename = "my_pkg.txt", initial.options = list( foo = 5 ) )
	checkTrue( getOptions( .mypkg, "foo" ) == 5 )
	checkTrue( getOptions( .mypkg )[[ "foo" ]] == 5 )
	saveOptions( .mypkg )
	rm( .mypkg )
	
	#check saved value
	.mypkg <<- new( "option", filename = "my_pkg.txt", initial.options = list( foo = 99999, new_var = 10 ) )
	checkTrue( getOptions( .mypkg, "foo" ) == 5 )
	checkTrue( getOptions( .mypkg, "new_var" ) == 10 )
	
	#set new filename
	new_opt_fname <- "new_opt.txt"
	unlink( fname )
	unlink( new_opt_fname )
	setOptionsFileName( .mypkg, new_opt_fname )
	saveOptions( .mypkg )
	checkTrue( file.exists( fname ) == FALSE )
	checkTrue( file.exists( new_opt_fname ) == TRUE )

	unlink( fname )
	unlink( new_opt_fname )
	
	saveOptions( .mypkg, "newest_opt.txt" )
	checkTrue( file.exists( fname ) == FALSE )
	checkTrue( file.exists( new_opt_fname ) == FALSE )
	checkTrue( file.exists( "newest_opt.txt" ) == TRUE)
	
}
