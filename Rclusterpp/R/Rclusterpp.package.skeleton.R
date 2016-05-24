# Adapted from Rclusterpp.package.skeleton.R in Rclusterpp package

Rclusterpp.package.skeleton <- function(
	name = "anRpackage", list = character(), environment = .GlobalEnv,
	path = ".", force = FALSE, namespace = TRUE, 
	code_files = character(), 
	example_code = TRUE ){
	
	env <- parent.frame(1)
	
	if( !length(list) ){
		fake <- TRUE
		assign( "Rcpp.fake.fun", function(){}, envir = env )
	} else {
		fake <- FALSE
	}
	
	# first let the traditional version do its business
	call <- match.call()
	call[[1]] <- as.name("package.skeleton")
	call[["namespace"]] <- namespace
	if( "example_code" %in% names( call ) ){
		# remove the example_code argument
		call[["example_code"]] <- NULL
	}
	if( fake ){
		call[["list"]] <- "Rcpp.fake.fun"
	}
	
	tryCatch( eval( call, envir = env ), error = function(e){
		stop( "error while calling `package.skeleton`" )
	} )
	
	message( "\nAdding Rclusterpp settings" )
	
	# now pick things up 
	root <- file.path( path, name )
	
	# Add Rcpp to the DESCRIPTION
	DESCRIPTION <- file.path( root, "DESCRIPTION" )
	if( file.exists( DESCRIPTION ) ){
		x <- cbind( read.dcf( DESCRIPTION ), 
			"Depends" = sprintf("Rcpp (>= %s), RcppEigen (>= %s), Rclusterpp (>= %s) ", 
				packageDescription("Rcpp")[["Version"]], 
				packageDescription("RcppEigen")[["Version"]], 
				packageDescription("Rclusterpp")[["Version"]]), 
			"LinkingTo" = "Rcpp, RcppEigen, Rclusterpp" )
		write.dcf( x, file = DESCRIPTION )
		message( " >> added Depends: Rcpp, RcppEigen, Rclusterpp" )
		message( " >> added LinkingTo: Rcpp, RcppEigen, Rclusterpp" )
	}
	
	# if there is a NAMESPACE, add a useDynLib
	NAMESPACE <- file.path( root, "NAMESPACE")
	if( file.exists( NAMESPACE ) ){
		lines <- readLines( NAMESPACE )
		if( ! grepl( "useDynLib", lines ) ){
			lines <- c( sprintf( "useDynLib(%s)", name), lines)
			writeLines( lines, con = NAMESPACE )
			message( " >> added useDynLib directive to NAMESPACE" )
		}
	}
	
	# lay things out in the src directory
	src <- file.path( root, "src")
	if( !file.exists( src )){
		dir.create( src )
	}
	skeleton <- system.file( "skeleton", package = "Rclusterpp" )
	
	copy <- function(path, overwrite=FALSE) {
		if (overwrite || !file.exists(path)) {
			file.copy( file.path( skeleton, basename(path) ), path)
			message( sprintf(" >> added %s file with Rclusterpp settings", basename(path)) )
		}
		path
	}

	copy( file.path(src, "Makevars.in") )
	copy( file.path(src, "Makevars.win") )
	Sys.chmod( copy( file.path(root, "cleanup") ), mode="0755")
	Sys.chmod( copy( file.path(root, "configure") ), mode="0755")
	writeLines(
		gsub("@PKG@", name, readLines( file.path(skeleton, "configure.ac") ), fixed=TRUE),
		file.path(root, "configure.ac")
	)

	if( example_code ){
		header <- readLines( file.path( skeleton, "rclusterpp_hello_world.h" ) )
		header <- gsub( "@PKG@", name, header, fixed = TRUE )
		writeLines( header , file.path( src, "rclusterpp_hello_world.h" ) )
		message( " >> added example header file using Rcpp/RcppEigen/Rclusterpp")
		
		file.copy( file.path( skeleton, "rclusterpp_hello_world.cpp" ), src )
		message( " >> added example src file using Rclusterpp classes")
		
		rcode <- readLines( file.path( skeleton, "rclusterpp_hello_world.R" ) )
		rcode <- gsub( "@PKG@", name, rcode, fixed = TRUE )
		writeLines( rcode , file.path( root, "R", "rclusterpp_hello_world.R" ) )
		message( " >> added example R file calling the C++ example")
	}
	if( fake ){
		rm( "Rcpp.fake.fun", envir = env )
		unlink( file.path( root, "R"  , "Rcpp.fake.fun.R" ) )
		unlink( file.path( root, "man", "Rcpp.fake.fun.Rd" ) )
		
	}
	
	invisible( NULL )
}


