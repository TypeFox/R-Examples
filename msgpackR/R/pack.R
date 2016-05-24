pack <-
function(data) {
	.pack_nil <- function() {
		return(as.raw(0xC0))
	}

	.pack_true <- function() {
		return(as.raw(0xC3))
	}

	.pack_false <- function() {
		return(as.raw(0xC2))
	}

	.pack_pfixnum <- function(value) {
		return(as.raw(0x7F) & as.raw(value))
	}

	.pack_nfixnum <- function(value) {
		return(xor(as.raw(0xFF), as.raw(abs(value)-1)))
	}

	.pack_uint8 <- function(value) {
		return(c(as.raw(0xCC), as.raw(value)))
	}

	.pack_uint16 <- function(value) {
		return(c(as.raw(0xCD), as.raw(value/2^8), as.raw(value%%2^8)))
	}

	.pack_uint32 <- function(value) {
		return(c(as.raw(0xCE), as.raw(value/2^24), as.raw((value%%2^24)/2^16), as.raw((value%%2^16)/2^8), as.raw(value%%2^8)))
	}

	.pack_uint64 <- function(value) {
		return(c(as.raw(0xCF), as.raw(value/2^56), as.raw((value%%2^56)/2^48), as.raw((value%%2^48)/2^40), as.raw((value%%2^40)/2^32), as.raw((value%%2^32)/2^24), as.raw((value%%2^24)/2^16), as.raw((value%%2^16)/2^8), as.raw(value%%2^8)))
	}

	.pack_int8 <- function(value) {
		if ( value < 0 ) {
			value <- 2^8 + value
		}
		return(c(as.raw(0xD0), as.raw(value)))
	}

	.pack_int16 <- function(value) {
		if ( value < 0 ) {
			value <- 2^16 + value
		}
		return(c(as.raw(0xD1), as.raw(value/2^8), as.raw(value%%2^8)))
	}

	.pack_int32 <- function(value) {
		if ( value < 0 ) {
			value <- 2^32 + value
		}
		return(c(as.raw(0xD2), as.raw(value/2^24), as.raw((value%%2^24)/2^16), as.raw((value%%2^16)/2^8), as.raw(value%%2^8)))
	}

	.pack_int64 <- function(value) {
		if ( value < 0 ) {
			value <- 2^64 + value
		}
		return(c(as.raw(0xD3), as.raw(value/2^56), as.raw((value%%2^56)/2^48), as.raw((value%%2^48)/2^40), as.raw((value%%2^40)/2^32), as.raw((value%%2^32)/2^24), as.raw((value%%2^24)/2^16), as.raw((value%%2^16)/2^8), as.raw(value%%2^8)))
	}

	.pack_float <- function(value) {
		# R supports not "float" but "double".
		return(.pack_double(value))
	}

	.pack_double <- function(value) {
		# the function to transform binary number into decimal number
		bit2dec <- function(bits) {
			n <- length(bits)
			result <- 0
			
			for ( i in 1:n ) {
				if ( bits[i] == 1 ) {
					result <- result + 2^(n-i)
				}
			}
			return(result)
		}
		
		# the function to transform fraction into 52 element vector 
		bit_fraction <- function(n) {
			result <- c()
			for ( i in 1:52 ) {
				b <- 0
				if ( n >= 2^(-i) ) {
					n <- n - 2^(-i)
					b <- 1
				}
				result <- c(result, b)
			}
			return(result)
		}
		# the end of the definition of two functions only used in ".pack_double()"
		
		result <- c(as.raw(0xCB))
		sign <- value<0
		exp <- floor(log2(abs(value)))
		bits <- bit_fraction(abs(value)*(2^(-exp))-1)
		exp <- exp + 1023
		
		# from 1st bit to 8th bit
		result <- c(result, as.raw(sign*0x80 + exp%/%2^4))
		# from 9th bit to 16th bit
		result <- c(result, as.raw(exp%%2^4*0x10 + bit2dec(bits[1:4])))
		# from 17th bit to 52nd bit
		result <- c(result, as.raw(bit2dec(bits[5:12])), as.raw(bit2dec(bits[13:20])), as.raw(bit2dec(bits[21:28])), as.raw(bit2dec(bits[29:36])), as.raw(bit2dec(bits[37:44])), as.raw(bit2dec(bits[45:52])))
		
		return(result)
	}

	.pack_map <- function(data) {
		# the function to return the datum and it's name
		.pack_map2 <- function(datum) {
			name <- names(datum)
			for ( d in datum ) {
				value <- d
			}
			return(c(pack(name), pack(value)))
		}
		# the end of the definition of two functions only used in ".pack_map()"
	
		result <- c()
		keys <- names(data)
		values <- as.list(data)
		n <- length(data)
		
		if ( n <= 15 ) {
			result <- c(result, (as.raw(0x8F) & as.raw(0x80+n)))
			for ( i in 1:n ) {
				result <- c(result, .pack_map2(data[i]))
			}
		}
		else if ( n <= 2^16-1 ) {
			result <- c(result, as.raw(0xDE), as.raw(n/2^8), as.raw(n%%2^8))
			for ( i in 1:n ) {
				result <- c(result, .pack_map2(data[i]))
			}
		}
		else if ( n <= 2^32-1 ) {
			result <- c(result, as.raw(0xDF), as.raw(n/2^24), as.raw((n%%2^24)/2^16), as.raw((n%%2^16)/2^8), as.raw(n%%2^8))
			for ( i in 1:n ) {
				result <- c(result, .pack_map2(data[i]))
			}
		}
		# if the length is more than 2^32-1, 
		else {
			# not implemented
		}
		return(result)
	}

	.pack_array <- function(data) {
		result <- c()
		n <- length(data)
		
		if ( n <= 15 ) {
			#FixArray
			result <- c(result, (as.raw(0x9F) & as.raw(0x90+n)))
			for ( datum in data ) {
				result <- c(result, pack(datum))
			}
		}
		else if ( n <= 2^16-1 ) {
			#array16
			result <- c(result, as.raw(0xDC), as.raw(n/2^8), as.raw(n%%2^8))
			for ( datum in data ) {
				result <- c(result, pack(datum))
			}
		}
		else if ( n <= 2^32-1 ) {
			#array32
			result <- c(result, as.raw(0xDD), as.raw(n/2^24), as.raw((n%%2^24)/2^16), as.raw((n%%2^16)/2^8), as.raw(n%%2^8))
			for ( datum in data ) {
				result <- c(result, pack(datum))
			}
		}
		# if the length is more than 2^32-1, 
		else {
			# not implemented
		}
		
		return(result)
	}

	.pack_raw <- function(data) {
		n <- length(data)
	
		if ( n <= 31 ) {
			h <- (as.raw(0xA0) | as.raw(n))
		}
		else if ( n <= 2^16-1 ) {
			h <- as.raw(0xDA)
		}
		else if ( n <= 2^32-1 ) {
			h <- as.raw(0xDB)
		}
		
		result <- c(h, data)
		return(result)
	}



	result <- c()	
	
	if ( is.vector(data) && (length(data) > 1 || !is.null(names(data))) ) {
		# if the vector has no name, pack to "array".
		if ( is.null(names(data)) ) {
			result <- c(result, .pack_array(data))
		}
		# if the vector has some name, pack to "map".
		else {
			result <- c(result, .pack_map(data))
		}
	}
	else if ( is.matrix(data) || is.data.frame(data) ) {
		result <- c(result, (as.raw(0x9F) & as.raw(0x90 + nrow(data))))
		# if the matrix has no name, pack to "array".
		if ( is.null(colnames(data)) ) {
			for ( i in 1:nrow(data) ) {
				result <- c(result, .pack_array(data[i,]))
			}
		}
		# if the matrix has some name, pack to "map".
		else {
			for ( i in 1:nrow(data) ) {
				result <- c(result, .pack_map(data[i,]))
			}
		}
	}
	# array
	else if ( is.array(data) ) {
		# not implemented
	}
	# factor
	else if ( is.factor(data) ) {
		.msgpack_factorToChar <- function(vec) {
		if ( is.factor(vec) ) {
			result <- as.character(vec)
		}
		return(result)
}
		# transfer factor into character
		mat <- .msgpack_factorToChar(data)
		result <- pack(mat)
	}
	# order
	else if ( is.ordered(data) ) {
		# not implemented
	}
	# list
	else if ( is.list(data) ) {
		result <- c(result, .pack_array(data))
	}
	# atomic data type
	else {		
		# NULL
		if ( is.null(data) ) {
			result <- c(result, .pack_nil())
		}
		# boolean
		else if ( mode(data) == mode(TRUE) && data == TRUE ) {
			result <- c(result, .pack_true())
		}
		else if ( mode(data) == mode(FALSE) && data == FALSE ) {
			result <- c(result, .pack_false())
		}
		# character
		else if ( is.character(data) ) {
			result <- c(result, .pack_raw(charToRaw(data)))
		}
		# integer
		else if ( floor(data) == data ) {
			# positive fixnum
			if ( data <= 127 && data >= 0 ) {
				result <- c(result, .pack_pfixnum(data))
			}
			# negative fixnum
			else if ( data <= -1 && data >= -32 ) {
				result <- c(result, .pack_nfixnum(data))
			}
			# int8
			else if ( data < -32 && data >= -2^7 ) {
				result <- c(result, .pack_int8(data))
			}
			# int16
			else if ( data < -2^7 && data >= -2^15 ) {
				result <- c(result, .pack_int16(data))
			}
			# int32
			else if ( data < -2^15 && data >= -2^31 ) {
				result <- c(result, .pack_int32(data))
			}
			# int64
			else if ( data < -2^31 && data >= -2^63 ) {
				result <- c(result, .pack_int64(data))
			}
			else if ( data < -2^63 ) {
				# not implemented
			}
			# uint8
			else if ( data < 2^8 ) {
				result <- c(result, .pack_uint8(data))
			}
			# uint16
			else if ( data < 2^16 ) {
				result <- c(result, .pack_uint16(data))
			}
			# uint32
			else if ( data < 2^32 ) {
				result <- c(result, .pack_uint32(data))
			}
			# uint64
			else if ( data < 2^64 ) {
				result <- c(result, .pack_uint64(data))
			}
			
			# the other cases...
			else {
				# not implemented
			}
		}
		# floating point number
		else if ( floor(data) != data ) {
			result <- .pack_double(data)
		}
	}
	
	return(result)
}
