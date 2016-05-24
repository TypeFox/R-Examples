unpack <-
function(str) {
	.unpack_bin <- function(bin) {
		e$.msgpack_index <- 1
		e$.msgpack_data <- bin
	
		return(.unpack_data())
	}

	.unpack_file <- function(filename) {
		fl <- file(filename, "rb")
		bits <- readBin(fl, raw(), file.info(filename)$size)
		close(fl)
	
		e$.msgpack_index <- 1
		e$.msgpack_data <- bits
		return(.unpack_data())
	}

	
	.unpack_pfixnum <- function() {
		num <- e$.msgpack_data[e$.msgpack_index]
		e$.msgpack_index <- e$.msgpack_index + 1
		return(as.integer(num & as.raw(0x7F)))
	}

	.unpack_fixmap <- function() {
		result <- list()
		N <- as.integer(e$.msgpack_data[e$.msgpack_index] & as.raw(0x0F))
		nms <- character(N)
		
		e$.msgpack_index <- e$.msgpack_index + 1
		
		for ( i in 1:N ) {
			nms[i] <- .unpack_data()
			result[[i]] <- .unpack_data()
		}
		result <- .unpack_checkclass(result)
		names(result) <- nms
		return(result)
	}

	.unpack_fixarray <- function() {
		result <- list()
		N <- as.integer(e$.msgpack_data[e$.msgpack_index] & as.raw(0x0F))
		
		e$.msgpack_index <- e$.msgpack_index + 1
		
		for ( i in 1:N ) {
			result[[i]] <- .unpack_data()
		}
		result <- .unpack_checkclass(result)
		
		return(result)
	}

	.unpack_fixraw <- function() {
		result <- ""
		result_byte <- list()
		N <- as.integer(e$.msgpack_data[e$.msgpack_index] & as.raw(0x1F))
		
		e$.msgpack_index <- e$.msgpack_index + 1
		
		for ( i in 1:N ) {
			result_byte[[i]] <- rawToChar(e$.msgpack_data[e$.msgpack_index])
			e$.msgpack_index <- e$.msgpack_index + 1
		}

		result_byte <- .unpack_checkclass(result_byte)
		
		for ( i in 1:length(result_byte) ) {
			result <- paste(result, result_byte[i], sep="")
		}

		return(result)
	}

	.unpack_nil <- function() {
		e$.msgpack_index <- e$.msgpack_index + 1
		return(NULL)
	}

	.unpack_false <- function() {
		e$.msgpack_index <- e$.msgpack_index + 1
		return(FALSE)
	}

	.unpack_true <- function() {
		e$.msgpack_index <- e$.msgpack_index + 1
		return(TRUE)
	}

	.unpack_float <- function() {
		# R supports not "float" but "double".
		return(.unpack_double())
	}

	.unpack_double <- function() {
		result <- 0
		e$.msgpack_index <- e$.msgpack_index + 1
	
		bits <- c()
		for ( i in 1:8 ) {
			bits[i] <- e$.msgpack_data[e$.msgpack_index]
			e$.msgpack_index <- e$.msgpack_index + 1
		}
	
		# sign
		sign <- ifelse((bits[1] & as.raw(0x80)) != as.raw(0x00), -1, 1)
		# exponent
		exp <- as.integer(bits[1] & as.raw(0x7F))*2^4 + as.integer(bits[2] & as.raw(0xF0))/2^4 - 1023
		# fraction
		frac <- (rev(c(rawToBits(bits[2] & as.raw(0x0F)))))[5:8]
		for ( i in 3:8 ) {
			frac <- c(frac, as.integer(rev(rawToBits(bits[i]))))
		}
	
		for ( i in 1:52 ) {
			result <- result + frac[i]/2^i
		}
		result <- sign*(result+1)*2^exp
		
		return(result)
	}

	.unpack_uint8 <- function() {
		result <- 0
		
		e$.msgpack_index <- e$.msgpack_index + 1
		result <- result + as.integer(e$.msgpack_data[e$.msgpack_index])
		
		e$.msgpack_index <- e$.msgpack_index + 1
		return(result)
	}

	.unpack_uint16 <- function() {
		result <- 0
		
		for ( i in seq(8,0,-8) ) {
			e$.msgpack_index <- e$.msgpack_index + 1
			result <- result + as.integer(e$.msgpack_data[e$.msgpack_index])*2^i
		}
		
		e$.msgpack_index <- e$.msgpack_index + 1
		return(result)
	}

	.unpack_uint32 <- function() {
		result <- 0
		
		for ( i in seq(24,0,-8) ) {
			e$.msgpack_index <- e$.msgpack_index + 1
			result <- result + as.integer(e$.msgpack_data[e$.msgpack_index])*2^i
		}
		
		e$.msgpack_index <- e$.msgpack_index + 1
		return(result)
	}

	.unpack_uint64 <- function() {
		result <- 0
		
		for ( i in seq(56,0,-8) ) {
			e$.msgpack_index <- e$.msgpack_index + 1
			result <- result + as.integer(e$.msgpack_data[e$.msgpack_index])*2^i
		}
		
		e$.msgpack_index <- e$.msgpack_index + 1
		return(result)
	}

	.unpack_int8 <- function() {
		result <- 0
		
		e$.msgpack_index <- e$.msgpack_index + 1
		sign <- ifelse((e$.msgpack_data[e$.msgpack_index] & as.raw(0x80)) != as.raw(0x00), -1, 1)
		
		result <- result + as.integer(e$.msgpack_data[e$.msgpack_index])
		e$.msgpack_index <- e$.msgpack_index + 1
		
		if ( sign < 0 ) {
			result <- -2^8 + result
		}
		
		return(result)
	}

	.unpack_int16 <- function() {
		result <- 0
		
		e$.msgpack_index <- e$.msgpack_index + 1
		sign <- ifelse((e$.msgpack_data[e$.msgpack_index] & as.raw(0x80)) != as.raw(0x00), -1, 1)
		
		for ( i in seq(8,0,-8) ) {
			result <- result + as.integer(e$.msgpack_data[e$.msgpack_index])*2^i
			e$.msgpack_index <- e$.msgpack_index + 1
		}
		
		if ( sign < 0 ) {
			result <- -2^16 + result
		}
		
		return(result)
	}

	.unpack_int32 <- function() {
		result <- 0
		
		e$.msgpack_index <- e$.msgpack_index + 1
		sign <- ifelse((e$.msgpack_data[e$.msgpack_index] & as.raw(0x80)) != as.raw(0x00), -1, 1)
		
		for ( i in seq(24,0,-8) ) {
			result <- result + as.integer(e$.msgpack_data[e$.msgpack_index])*2^i
			e$.msgpack_index <- e$.msgpack_index + 1
		}
		
		if ( sign < 0 ) {
			result <- -2^32 + result
		}
		
		return(result)
	}

	.unpack_int64 <- function() {
		result <- 0
		
		e$.msgpack_index <- e$.msgpack_index + 1
		sign <- ifelse((e$.msgpack_data[e$.msgpack_index] & as.raw(0x80)) != as.raw(0x00), -1, 1)
		
		for ( i in seq(56,0,-8) ) {
			result <- result + as.integer(e$.msgpack_data[e$.msgpack_index])*2^i
			e$.msgpack_index <- e$.msgpack_index + 1
		}
		
		if ( sign < 0 ) {
			result <- -2^64 + result
		}
		
		return(result)
	}

	.unpack_raw16 <- function() {
		result <- list()
		N <- 0
		
		e$.msgpack_index <- e$.msgpack_index + 1
		
		for ( i in seq(8,0,-8) ) {
			N <- N + as.integer(e$.msgpack_data[e$.msgpack_index])*2^i
			e$.msgpack_index <- e$.msgpack_index + 1
		}
		
		for ( i in 1:N ) {
			result[[i]] <- rawToChar(e$.msgpack_data[e$.msgpack_index])
			e$.msgpack_index <- e$.msgpack_index + 1
		}
		result <- .unpack_checkclass(result)
		
		return(result)
	}

	.unpack_raw32 <- function() {
		result <- list()
		N <- 0
		
		e$.msgpack_index <- e$.msgpack_index + 1
		
		for ( i in seq(24,0,-8) ) {
			N <- N + as.integer(e$.msgpack_data[e$.msgpack_index])*2^i
			e$.msgpack_index <- e$.msgpack_index + 1
		}
		
		for ( i in 1:N ) {
			result[[i]] <- rawToChar(e$.msgpack_data[e$.msgpack_index])
			e$.msgpack_index <- e$.msgpack_index + 1
		}
		result <- .unpack_checkclass(result)
		
		return(result)
	}

	.unpack_array16 <- function() {
		result <- list()
		N <- 0
		
		e$.msgpack_index <- e$.msgpack_index + 1
		
		for ( i in seq(8,0,-8) ) {
			N <- N + as.integer(e$.msgpack_data[e$.msgpack_index])*2^i
			e$.msgpack_index <- e$.msgpack_index + 1
		}
		
		for ( i in 1:N ) {
			result[[i]] <- .unpack_data()
		}
		result <- .unpack_checkclass(result)
		
		return(result)
	}

	.unpack_array32 <- function() {
		result <- list()
		N <- 0
		
		e$.msgpack_index <- e$.msgpack_index + 1
		
		for ( i in seq(24,0,-8) ) {
			N <- N + as.integer(e$.msgpack_data[e$.msgpack_index])*2^i
			e$.msgpack_index <- e$.msgpack_index + 1
		}
		
		for ( i in 1:N ) {
			result[[i]] <- .unpack_data()
		}	
		result <- .unpack_checkclass(result)
		
		return(result)
	}

	.unpack_map16 <- function() {
		result <- list()
		N <- 0
		
		e$.msgpack_index <- e$.msgpack_index + 1
		
		for ( i in seq(8,0,-8) ) {
			N <- N + as.integer(e$.msgpack_data[e$.msgpack_index])*2^i
			e$.msgpack_index <- e$.msgpack_index + 1
		}
		
		nms <- character(N)
		
		for ( i in 1:N ) {
			nms[i] <- .unpack_data()
			result[[i]] <- .unpack_data()
		}
		result <- .unpack_checkclass(result)
		names(result) <- nms
		
		return(result)
	}

	.unpack_map32 <- function() {
		result <- list()
		N <- 0
		
		e$.msgpack_index <- e$.msgpack_index + 1
		
		for ( i in seq(24,0,-8) ) {
			N <- N + as.integer(e$.msgpack_data[e$.msgpack_index])*2^i
			e$.msgpack_index <- e$.msgpack_index + 1
		}
		
		nms <- character(N)
			
		for ( i in 1:N ) {
			nms[i] <- .unpack_data()
			result[[i]] <- .unpack_data()
		}
		result <- .unpack_checkclass(result)
		names(result) <- nms
		
		return(result)
	}

	.unpack_nfixnum <- function() {
		num <- e$.msgpack_data[e$.msgpack_index]
		e$.msgpack_index <- e$.msgpack_index + 1
		return(-32 + as.integer(num & as.raw(0x1F)))
	}

	.unpack_checkclass <- function(data) {
		classes <- unique(unlist(lapply(data, class)))

		for ( i in 1:length(classes) ) {
			if ( classes[i] == "integer" ) {
				classes[i] = "numeric"
			}
		}
		classes <- unique(classes)
		
		len <- unique(unlist(lapply(data, length)))
		if( length(classes) == 1 && length(len) == 1 && len[1] == 1) {
			class(data) <- classes
		}
		
		return(data)
	}

	.unpack_data <- function() {
		if ( is.null(e$.msgpack_data) ) {
			return(NULL)
		}
		# Positive FixNum
		else if( e$.msgpack_data[e$.msgpack_index] <= as.raw(0x7F) ) {
			return(.unpack_pfixnum())
		}
		# FixMap
		else if ( e$.msgpack_data[e$.msgpack_index] <= as.raw(0x8F) ) {
			return(.unpack_fixmap())
		}
		# FixArray
		else if ( e$.msgpack_data[e$.msgpack_index] <= as.raw(0x9F) ) {
			return(.unpack_fixarray())
		}
		# FixRaw
		else if ( e$.msgpack_data[e$.msgpack_index] <= as.raw(0xBF) ) {
			return(.unpack_fixraw())
		}
		# nil
		else if ( e$.msgpack_data[e$.msgpack_index] == as.raw(0xC0) ) {
			return(.unpack_nil())
		}
		# (reserved)
		else if ( e$.msgpack_data[e$.msgpack_index] == as.raw(0xC1) ) {
			return(NULL)
		}
		# false
		else if ( e$.msgpack_data[e$.msgpack_index] == as.raw(0xC2) ) {
			return(.unpack_false())
		}
		# true
		else if ( e$.msgpack_data[e$.msgpack_index] == as.raw(0xC3) ) {
			return(.unpack_true())
		}
		# (reserved)
		else if ( e$.msgpack_data[e$.msgpack_index] <= as.raw(0xC9) ) {
			return(NULL)
		}
		# float
		else if ( e$.msgpack_data[e$.msgpack_index] == as.raw(0xCA) ) {
			return(.unpack_float())
		}
		# double
		else if ( e$.msgpack_data[e$.msgpack_index] == as.raw(0xCB) ) {
			return(.unpack_double())
		}
		# uint8
		else if ( e$.msgpack_data[e$.msgpack_index] == as.raw(0xCC) ) {
			return(.unpack_uint8())
		}
		# uint16
		else if ( e$.msgpack_data[e$.msgpack_index] == as.raw(0xCD) ) {
			return(.unpack_uint16())
		}
		# uint32
		else if ( e$.msgpack_data[e$.msgpack_index] == as.raw(0xCE) ) {
			return(.unpack_uint32())
		}
		# uint64
		else if ( e$.msgpack_data[e$.msgpack_index] == as.raw(0xCF) ) {
			return(.unpack_uint64())
		}
		# int8
		else if ( e$.msgpack_data[e$.msgpack_index] == as.raw(0xD0) ) {
			return(.unpack_int8())
		}
		# int16
		else if ( e$.msgpack_data[e$.msgpack_index] == as.raw(0xD1) ) {
			return(.unpack_int16())
		}
		# int32
		else if ( e$.msgpack_data[e$.msgpack_index] == as.raw(0xD2) ) {
			return(.unpack_int32())
		}
		# int64
		else if ( e$.msgpack_data[e$.msgpack_index] == as.raw(0xD3) ) {
			return(.unpack_int64())
		}
		# (reserved)
		else if ( e$.msgpack_data[e$.msgpack_index] <= as.raw(0xD9) ) {
			return(NULL)
		}
		# raw16
		else if ( e$.msgpack_data[e$.msgpack_index] == as.raw(0xDA) ) {
			return(.unpack_raw16())
		}
		# raw32
		else if ( e$.msgpack_data[e$.msgpack_index] == as.raw(0xDB) ) {
			return(.unpack_raw32())
		}
		# array16
		else if ( e$.msgpack_data[e$.msgpack_index] == as.raw(0xDC) ) {
			return(.unpack_array16())
		}
		# array32
		else if ( e$.msgpack_data[e$.msgpack_index] == as.raw(0xDD) ) {
			return(.unpack_array32())
		}
		# map16
		else if ( e$.msgpack_data[e$.msgpack_index] == as.raw(0xDE) ) {
			return(.unpack_map16())
		}
		# map32
		else if ( e$.msgpack_data[e$.msgpack_index] == as.raw(0xDF) ) {
			return(.unpack_map32())
		}
		# Negative FixNum
		else {
			return(.unpack_nfixnum())
		}
	}

	e <- new.env()
	e$.msgpack_index <- 1
	e$.msgpack_data <- list()
	# the case when str is filename
	if ( mode(str) == "character" ) {
		return(.unpack_file(str))
	}
	# the case when str is raw array
	else if ( mode(str) == "raw" ) {
		return(.unpack_bin(str))
	}
	else {

	}
}
