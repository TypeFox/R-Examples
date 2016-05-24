# Unit test for stat functions.
#
# Copyright (C) 2011 Renaud Gaujoux
# 
# This file is part of RcppOctave.
#
# RcppOctave is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# RcppOctave is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with RcppOctave.  If not, see <http://www.gnu.org/licenses/>.

check.conversion <- function(t, msg=NULL){
	lapply(t, function(x){
		d <- paste(if( is.null(dim(x)) ) length(x) else dim(x), collapse=' x ')
		checkIdentical( o_identity(x), x
					, paste(msg, " - Conversions to and from ", class(x)
					, "<", d, "> [", paste(head(x), collapse=', '), "] is OK", sep=''))
	})
}

check.dimnames_dropped <- function(x, msg){
	x <- add_names(x)
	ox <- o_identity(x)
	checkIdentical(ox, setNames(x, NULL), paste(msg, " - Values of named vector converted properly"))
	checkIdentical(names(ox), NULL, paste(msg, " - Names of named vector are dropped"))
}

add_names <- function(x) setNames(x, letters[1:length(x)])

# Unit test for logical conversions
test.bool <- function(){
	
	bool <- c(TRUE, FALSE)
	
	set.seed(1234)
	t <- list(TRUE
			, FALSE
			, c(TRUE, FALSE, TRUE, TRUE, FALSE)
			, matrix(sample(bool, 50, rep=TRUE), 10, 5)
			, logical(0)
	)
	
	check.conversion(t)
	
	# check dimnames
	check.dimnames_dropped(bool, "boolean")
}

# Unit test for string conversions
test.char <- function(){
	
	char <- c('absjdkslj', 'sewae', '!@#$%^&(*)')
	
	set.seed(1234)
	t <- list( char[1]
			, ""
			, char
			, character(0)
	)
	
	check.conversion(t)
	
	# check dimnames
	check.dimnames_dropped(char, "character")
	
	checkException(o_identity(matrix(sample(char, 50, rep=TRUE), 10, 5)), "Character matrices are not supported")
}

# Unit test for double conversions
test.double <- function(){
	
	set.seed(1234)
	t <- list(0.0
			, 1.0
			, numeric(0)
			, 94.23
			, -123.323
			, runif(10)
			, matrix(runif(50), 10, 5)			
	)
	
	check.conversion(t)
	
	# check dimnames
	check.dimnames_dropped(runif(10), "double")
}

# Unit test for integer conversions
test.int <- function(){
	
	set.seed(1234)
	t <- list(0L
			, 1L
			, integer(0)
			, 100L
			, -3L
			, 1:10
			, matrix(1:50, 10, 5)
	)
		
	check.conversion(t)
	
	# check dimnames
	check.dimnames_dropped(25:40, "integer")
}

# Unit test for list conversions
test.list <- function(){
	
	set.seed(1234)
		
	t <- list( 
			list()
			, list('z', runif(10))
			, list(1L, 2L, 3L)
			, list(1, 2, 3)
			, list(1:2, 2)
			, list('x', 'y', 'z')
			, list(1, 2, 3, runif(5), letters[10:17])
			, list(matrix(1:10, 2,5), matrix(1:10, 5, 2), matrix(c(TRUE, FALSE), 2,5))
	)
	
	check.conversion(t, "Unamed lists")
	t <- t[-1]
	
	tn <- lapply(t, function(x){ names(x) <- letters[1:length(x)]; x} )	
	check.conversion(tn, "Named lists")
	
	tn.last <- lapply(t, function(x){ names(x) <- letters[1:length(x)]; names(x)[length(x)] <- ""; x} )	
	check.conversion(tn.last, "Lists with last name missing")
	
	tn.first <- lapply(t, function(x){ names(x) <- letters[1:length(x)]; names(x)[1] <- ""; x} )	
	check.conversion(tn.first, "Lists with first name missing")
	
	tn.middle <- lapply(t, function(x){ names(x) <- letters[1:length(x)]; names(x)[2] <- ""; x} )
	check.conversion(tn.middle, "Lists with middle name missing")
		
	checkException(o_identity(list(a=1, b=5, 2, 3)), "More than 1 missing name throws an error")
}