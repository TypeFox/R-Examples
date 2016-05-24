library('testthat')
context("createSparseGrid")

test_that("check for sensible input values", {
    correct_type1 <- "GQU"
    correct_dimension <- 3
    correct_k <- 4
    correct_sym <- TRUE
    
    wrong_type1 <- "HGE"
    wrong_type2 <- function() { return( list( "nodes"=c(1), "weights"=c(1) ) ) }			# does not contain argument
    wrong_type3 <- function(n1, n2) { return( list( "nodes"=c(1), "weights"=c(1) ) ) }	    # contains too many arguments
    wrong_type4 <- function(n) { return( c(1) ) }											# does not return list
    wrong_type5 <- function(n) { return( list( "weights"=c(1) ) ) }							# return list does not contain "nodes"
    wrong_type6 <- function(n) { return( list( "nodes"=c(1) ) ) }						    # return list does not contain "weights"
    
	wrong_dimension <- 3.1
	wrong_k         <- 4.1
	wrong_sym1      <- NA
	wrong_sym2      <- 1
	
    expect_that(createSparseGrid( type=wrong_type1,   dimension=correct_dimension, k=correct_k, sym=correct_sym ), throws_error())
    expect_that(createSparseGrid( type=wrong_type2,   dimension=correct_dimension, k=correct_k, sym=correct_sym ), throws_error())
    expect_that(createSparseGrid( type=wrong_type3,   dimension=correct_dimension, k=correct_k, sym=correct_sym ), throws_error())
    expect_that(createSparseGrid( type=wrong_type4,   dimension=correct_dimension, k=correct_k, sym=correct_sym ), throws_error())
    expect_that(createSparseGrid( type=wrong_type5,   dimension=correct_dimension, k=correct_k, sym=correct_sym ), throws_error())
    expect_that(createSparseGrid( type=wrong_type6,   dimension=correct_dimension, k=correct_k, sym=correct_sym ), throws_error())
    expect_that(createSparseGrid( type=correct_type1, dimension=wrong_dimension,   k=correct_k, sym=correct_sym ), throws_error())
    expect_that(createSparseGrid( type=correct_type1, dimension=correct_dimension, k=wrong_k,   sym=correct_sym ), throws_error())
    expect_that(createSparseGrid( type=correct_type1, dimension=correct_dimension, k=correct_k, sym=wrong_sym1 ),  throws_error())
    expect_that(createSparseGrid( type=correct_type1, dimension=correct_dimension, k=correct_k, sym=wrong_sym2 ),  throws_error())
  
})

# Output values taken from www.sparse-grids.de
test_that("check output values against grids from www.sparse-grids.de", {
	
	# find files in sub-directory, and get all types, dimensions and accuracies from the filenames
	filelist 		<- dir( 'testgrids' )
	typelist        <- substr( filelist, 1, 3 )
	dimlist 		<- as.integer( gsub( "[[:upper:]][[:upper:]][[:upper:]]_d([0-9]+)_l([0-9]+).asc", "\\1", filelist ) )
	accuracylist 	<- as.integer( gsub( "[[:upper:]][[:upper:]][[:upper:]]_d([0-9]+)_l([0-9]+).asc", "\\2", filelist ) )

	for (cnt in 1:length(dimlist))
	{
		type            <- typelist[ cnt ]
		dimension		<- dimlist[ cnt ]
		accuracy		<- accuracylist[ cnt ]
		filename 		<- paste( 'testgrids/', type, "_d", dimension, "_l", accuracy, ".asc", sep='' )
		
		cat( "Testing: ", filename, '\n', sep='' )
	
		grid.read       <- readASCGrid( filename, dimension )
		expect_that( createSparseGrid( type, dimension=dimension, k=accuracy ), equals( grid.read ) )
	}
})
