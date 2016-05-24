 
trait_tests <- ( function () { 
	# create an environment containing trait-
	# testing functions, for speed of access.
	# require_a checks traits using this environment as a 
	# hash table.

	test_for <- new.env(
		parent = emptyenv()
	)
	
	# useful for testing purposes
	test_for$any = function (value) TRUE

	# (mostly) builtin functions,
	# tht mostly test the class of the object
	test_for$array = is.array
	test_for$atomic = is.atomic
	test_for$call =  is.call
	test_for$character = is.character
	test_for$complex =  is.complex
	test_for$data.frame =  is.data.frame
	test_for$double = is.double
	test_for$environment = is.environment
	test_for$expression =  is.expression
	test_for$factor = is.factor
	test_for$finite =  
		function (value) {
			is.numeric(value) && is.finite(value)
		}
	test_for$'function' = is.function
	test_for$infinite =  
		function (value) {
			is.numeric(value) && is.infinite(value)
		}
	test_for$integer = is.integer
	test_for$language = is.language
	test_for$list = is.list
	test_for$logical = is.logical
	test_for$matrix = is.matrix
	test_for$na = 
		function (value) {
			is.vector(value) && 
			!is.expression(value) && is.na(value)
		}
	test_for$name = is.name
	test_for$nan = 
		function (value) {
			is.numeric(value) && is.nan(value)
		}
	test_for$null = is.null
	test_for$numeric = is.numeric
	test_for$object = 
		function (value) {
			# a decent test for objectness,
			# since the built-in is for internal use only
			!is.null(attr(value, 'class'))
		}
	test_for$pairlist = is.pairlist
	test_for$primitive = is.primitive
	test_for$raw = is.raw
	test_for$recursive = is.recursive
	test_for$s4 = isS4
	test_for$symbol = is.symbol
	test_for$true = isTRUE
	test_for$table = is.table
	test_for$vector = is.vector

	# tests I find useful
	test_for$functionable = 
		function (value) {
			(is.character(value) && length(value) == 1) ||
			is.function(value) ||
			is.symbol(value)
		}
	test_for$false = 
		function (value) {
			is.logical(value) && !is.na(value) && !value
		}
	test_for$closure = 
		function (value) {
			# is value a normal R functions?
			is.function(value) && !is.primitive(value)
		}
	test_for$whole = 
		function (value) {
			# is value a whole number, 
			# within double precision limits?

			is.numeric(value) && 
			is.finite(value) &&
			( is.integer(value) || 
				(abs(round(value) - value) < .Machine$double.eps) )
		}
	test_for$positive = 
		function (value) {
			is.numeric(value) && 
			!is.nan(value) && value > 0
		}
	test_for$nonnegative =
		function (value) {
			is.numeric(value) && 
			!is.nan(value) && value >= 0
		}
	test_for$named = 
		function (value) {
			# does a listish value have names

			( is.vector(value) || is.pairlist(value) ) &&
			!is.expression(value) &&
			(length(value) == 0 || (
				!is.null(names(value)) &&
				all(nchar(names(value)) > 0)))
		}
	test_for$boolean =
		function (value) {
			is.logical(value) && !is.na(value)
		}
	test_for$string = 
		function (value) {
			is.character(value) && length(value) == 1
		}
	test_for$listy = 
		function (value) {
			# is the value a list, vector or 
			# pairlist?

			is.vector(value) || 
			is.list(value) || 
			is.pairlist(value)
		}

	# quantifiers
	test_for$length_zero = 
		function (value) {
			length(value) == 0
		}
	test_for$length_one = 
		function (value) {
			length(value) == 1
		}
	test_for$length_two = 
		function (value) {
			length(value) == 2
		}
	test_for$length_three = 
		function (value) {
			length(value) == 3
		}

	# from my library arrow. consistent way
	# of checking function parameters/formals.
	# won't be attatched to the test_for environment

	xParams <- function (f) {
		if (is.primitive(f)) {
			head( as.list(args(f)), -1 )
		} else {
			formals(f)
		}
	}

	# check function arity. variadic is always
	# the desired arity

	test_for$nullary = 
		function (value) {
			!is.function(value) || {
				params <- xParams(value)

				"..." %in% names(params) ||
				length(params) == 0
			}
		}
	test_for$unary =  
		function (value) {
			!is.function (value) || {
				params <- xParams(value)

				"..." %in% names(params) ||
				length(params) == 1
			}
		}
	test_for$binary =
		function (value) {
			!is.function(value) || {
				params <- xParams(value)

				"..." %in% names(params) ||
				length(params) == 2
			}
		}
	test_for$ternary = 
		function (value) {
			!is.function(value) || {
				params <- xParams(value)

				"..." %in% names(params) ||
				length(params) == 3
			}
		}
	test_for$variadic =
		function (value) {
			!is.function(value) || {
				params <- xParams(value)
				"..." %in% names(params)				
			}
		}

	test_for$valid_traits <- ls(test_for)
	test_for

} )()
