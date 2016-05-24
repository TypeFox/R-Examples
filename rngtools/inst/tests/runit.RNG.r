# Unit test for getRNG
# 
# Author: Renaud Gaujoux
###############################################################################

library(stringr)

test.getRNG <- function(){

	RNGkind('default', 'default')
	on.exit( RNGrecovery() )
	
	checker <- function(x, y, ..., msg=NULL, drawRNG=TRUE){
		
		if( drawRNG ) runif(10)
		fn <- getRNG
		oldRNG <- RNGseed()
		if( !missing(x) ){
			d <- fn(x, ...)
			cl <- str_c(class(x), '(', length(x), ')')
		}else{
			d <- fn()
			cl <- 'MISSING'
		}
		newRNG <- RNGseed()
		.msg <- function(x) paste(cl, ':', x, '[', msg, ']')
		checkIdentical(oldRNG, newRNG, .msg("does not change RNG"))
		checkIdentical(d, y, .msg("result is correct") )
	}
	
	set.seed(123456)
	seed123456 <- .Random.seed
	checker(, seed123456, msg="No arguments: returns .Random.seed", drawRNG=FALSE)
	checker(123456, seed123456, msg="Single numeric argument: returns .Random.seed as it would be after setting the seed")
	checker(123456, 123456, num.ok=TRUE, msg="Single numeric argument + num.ok: returns argument unchanged")
	checker(.Random.seed, .Random.seed, msg="Integer seed argument: returns its argument unchanged")
	checker(as.numeric(.Random.seed), .Random.seed, msg="Numeric seed argument: returns its argument as an integer vector")
	checker(2:3, 2:3, msg="Integer INVALID seed vector argument: returns its argument unchanged")
	checker(c(2,3), c(2L,3L), msg="Numeric INVALID seed vector argument: returns its argument as an integer vector")
	checker(1L, 1L, msg="Single integer = Encoded RNG kind: returns it unchanged")
	checker(1000L, 1000L, msg="Invalid single integer = Encoded RNG kind: returns it unchanged")
	
}

test.setRNG <- function(){
	
	RNGkind('default', 'default')
	on.exit( RNGrecovery() )
	
	checker <- function(x, y, tset, drawRNG=TRUE){
	
		on.exit( RNGrecovery() )
		
		if( drawRNG ) runif(10)
		oldRNG <- RNGseed()
		d <- force(x)
		newRNG <- RNGseed()
		
		msg <- function(x, ...) paste(tset, ':', ...)
		checkTrue(!identical(oldRNG, newRNG), msg("changes RNG"))
		checkIdentical(getRNG(), y, msg("RNG is correctly set") )
		checkIdentical(d, oldRNG, msg("returns old RNG") )
	}
	
	set.seed(123456)
	refseed <- .Random.seed
	checker(setRNG(123456), refseed, "Single numeric: sets current RNG with seed")
	
	# setting kind with a character string
	set.seed(123)
	RNGkind('Mar')
	refseed <- .Random.seed
	RNGrecovery()
	set.seed(123)
	checker(setRNG('Mar'), refseed, "Single character: change RNG kind", drawRNG=FALSE)
	
	# setting kind with a character string
	set.seed(123)
	RNGkind('Mar', 'Ahrens')
	refseed <- .Random.seed
	RNGrecovery()
	set.seed(123)
	checker(setRNG('Mar', 'Ahrens'), refseed, "Two character strings: change RNG kind and normal kind", drawRNG=FALSE)
	RNGrecovery()
	set.seed(123)
	checker(setRNG(c('Mar', 'Ahrens')), refseed, "2-long character vector: change RNG kind and normal kind", drawRNG=FALSE)
	
	# setting kind
	set.seed(123456, kind='Mar')
	refseed <- .Random.seed
	checker(setRNG(123456, kind='Mar'), refseed, "Single numeric + kind: change RNG kind + set seed")
	
	# setting Nkind
	set.seed(123456, normal.kind='Ahrens')
	refseed <- .Random.seed
	checker(setRNG(123456, normal.kind='Ahrens'), refseed
				, "Single numeric + normal.kind: change RNG normal kind + set seed")
	
	# setting kind and Nkind
	set.seed(123456, kind='Mar', normal.kind='Ahrens')
	refseed <- .Random.seed
	checker(setRNG(123456, kind='Mar', normal.kind='Ahrens'), refseed
			, "Single numeric + kind + normal.kind: change RNG all kinds + set seed")
	
	# with seed length > 1
	refseed <- as.integer(c(201, 0, 0))
	checker(setRNG(refseed), refseed, "numeric vector: directly set seed")
	
	refseed <- .Random.seed
	checkException( setRNG(2:3), "numeric vector: throws an error if invalid value for .Random.seed")
	checkIdentical( .Random.seed, refseed, ".Random.seed is not changed in case of an error in setRNG")
    
    oldRNG <- getRNG()
    checkException(setRNG(1234L), "Error with invalid integer seed")
    checkIdentical(oldRNG, getRNG(), "RNG still valid after error")
    checkException(setRNG(123L), "Error with invalid RNG kind")
    checkIdentical(oldRNG, getRNG(), "RNG still valid after error")

    # changes in R >= 3.0.2: invalid seeds only throw warning
    if( testRversion('> 3.0.1') ){
        oldRNG <- getRNG()
        checkWarning(setRNG(1234L, check = FALSE), "\\.Random\\.seed.* is not .* valid"
                        , "Invalid integer kind: Warning only if check = FALSE")
        checkIdentical(1234L, getRNG(), "RNG has new invalid integer value")
        setRNG(oldRNG)
        checkWarning(setRNG(123L, check = FALSE), "\\.Random\\.seed.* is not .* valid"
                        , "Invalid kind: Warning only if check = FALSE")
        checkIdentical(123L, getRNG(), "RNG has new invalid RNG kind")
                                                
    }
	
}

