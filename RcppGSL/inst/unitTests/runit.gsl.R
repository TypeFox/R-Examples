#!/usr/bin/r -t
#                        Emacs make this -*- mode: R; tab-width: 4 -*-
#
# Copyright (C) 2010 - 2013  Romain Francois and Dirk Eddelbuettel
#
# This file is part of RcppGSL.
#
# RcppGSL is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# RcppGSL is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with RcppGSL.  If not, see <http://www.gnu.org/licenses/>.

.setUp <- function() {
    if (!exists("pathRcppTests"))
        pathRcppTests <- system.file("unitTests", package="RcppGSL")
    sourceCpp(file.path(pathRcppTests, "cpp/gsl.cpp"))
}

test.gsl.vector.wrappers <- function(){
	fx <- test_gsl_vector_wrapper
	res <- fx()
	checkEquals(res,
                list("gsl_vector" = numeric(10),
                     "gsl_vector_float" = numeric(10),
                     "gsl_vector_int" = integer(10),
                     ##"gsl_vector_long" = numeric(10),
                     "gsl_vector_char" = raw(10),
                     "gsl_vector_complex" = complex(10),
                     "gsl_vector_complex_float" = complex(10),
                     "gsl_vector_complex_long_double" = complex(10),
                     "gsl_vector_long_double" = numeric(10),
                     "gsl_vector_short" = integer(10),
                     "gsl_vector_uchar" = raw(10),
                     "gsl_vector_uint" = integer(10),
                     "gsl_vector_ushort" = integer(10)
                     ##,"gsl_vector_ulong" = numeric(10)
                     ),
                msg = "wrap( gsl_vector )" )
}

test.gsl.vector <- function(){
    fx <- test_gsl_vector
    res <- fx()
    checkEquals(res,
                list("gsl_vector" = numeric(10),
                     "gsl_vector_float" = numeric(10),
                     "gsl_vector_int" = integer(10),
                     ##"gsl_vector_long" = numeric(10),
                     "gsl_vector_char" = raw(10),
                     "gsl_vector_complex" = complex(10),
                     "gsl_vector_complex_float" = complex(10),
                     "gsl_vector_complex_long_double" = complex(10),
                     "gsl_vector_long_double" = numeric(10),
                     "gsl_vector_short" = integer(10),
                     "gsl_vector_uchar" = raw(10),
                     "gsl_vector_uint" = integer(10),
                     "gsl_vector_ushort" = integer(10)
                     ##,"gsl_vector_ulong" = numeric(10)
                     ),
                msg = "wrap( gsl_vector )" )
}

test.gsl.matrix <- function(){
	helper <- function(what){
		as.what <- get( paste( "as.", deparse(substitute(what)), sep = "" ) )
		x <- what(10)
		x[1] <- as.what(1)
		x[7] <- as.what(1)
		dim( x )  <- c(5,2)
		x
	}
    fx <- test_gsl_matrix
    res <- fx()
	checkEquals(res,
                list("gsl_matrix"                     = helper( numeric ),
                     "gsl_matrix_float"               = helper( numeric ),
                     "gsl_matrix_int"                 = helper( integer ),
                     ##"gsl_matrix_long"                = helper( numeric ),
                     "gsl_matrix_char"                = helper( raw ),
                     "gsl_matrix_complex"             = helper( complex ),
                     "gsl_matrix_complex_float"       = helper( complex ),
                     "gsl_matrix_complex_long_double" = helper( complex ),
                     "gsl_matrix_long_double"         = helper( numeric ),
                     "gsl_matrix_short"               = helper( integer ),
                     "gsl_matrix_uchar"               = helper( raw ),
                     "gsl_matrix_uint"                = helper( integer ),
                     "gsl_matrix_ushort"              = helper( integer )
                     ##,"gsl_matrix_ulong"               = helper( numeric )
                     ),
                msg = "wrap( gsl_matrix )" )

}

test.gsl.vector.view <- function(){
    fx <- test_gsl_vector_view
    res <- fx()
	checkEquals(res,
                list( even = 2.0 * 0:4, odd = 2.0 * 0:4 + 1.0 ),
                msg = "wrap( gsl.vector.view )" )

    fx <- test_gsl_vector_view_wrapper
    res <- fx()
    checkEquals( res,
                list( even = 2.0 * 0:4, odd = 2.0 * 0:4 + 1.0 ),
                msg = "wrap( gsl.vector.view.wrapper )" )
}

test.gsl.matrix.view <- function(){
    fx <- test_gsl_matrix_view
    res <- fx()
	checkEquals( res$full[3:4, 3:4], res$view, msg = "wrap(gsl.matrix.view)" )

    fx <- test_gsl_matrix_view_wrapper
    res <- fx()
	checkEquals( res$full[3:4, 3:4], res$view, msg = "wrap(gsl.matrix.view.wrapper)" )
}

test.gsl.vector.input.SEXP <- function(){
	x <- rnorm( 10 )
    fx <- test_gsl_vector_input
    res <- fx(x)
	checkEquals( res, sum(x), msg = "RcppGSL::vector<double>(SEXP)" )
}

test.gsl.matrix.input.SEXP <- function(){
	x <- matrix( rnorm(20), nc = 4 )
    fx <- test_gsl_matrix_input
    res <- fx( x)
	checkEquals( res, sum(x[,1]), msg = "RcppGSL::matrix<double>(SEXP)" )
}

test.gsl.RcppGSL.vector <- function(){
    fx <- test_gsl_vector_conv
    res <- fx()
	checkEquals( res, 0:9, msg = "RcppGSL::vector<int> -> IntegerVector" )
}

test.gsl.RcppGSL.vector.indexing <- function(){
    fx <- test_gsl_vector_indexing
    res <- fx( seq(0.5, 10.5) )
	checkEquals( res, seq( 1.5, 11.5 ) )
}

test.gsl.RcppGSL.vector.iterating <- function(){
	x   <-  seq(0.5, 10.5)
    fx <- test_gsl_vector_iterating
    res <- fx(x)
	checkEquals( res, sum(x) )
}

test.gsl.RcppGSL.vector.iterator.transform <- function() {
	x   <-  seq(0.5, 10.5)
    fx <- test_gsl_vector_iterator_transform
    res <- fx(x)
	checkEquals(res, sqrt(x))
}

test.gsl.RcppGSL.matrix.indexing <- function(){
	m   <- matrix( 1:16+.5, nr = 4 )
    fx <- test_gsl_matrix_indexing
    res <- fx(m)
	checkEquals( res, m+1 )
}

test.gsl.RcppGSL.vector.view.iterating <- function(){
	x   <-  seq(1.5, 10.5)
    fx <- test_gsl_vector_view_iterating
    res <- fx(x)
	checkEquals( res, sum( x[ seq(1, length(x), by = 2 ) ] ) )
}

test.gsl.RcppGSL.matrix.view.indexing <- function(){
    fx <- test_gsl_matrix_view_indexing
    res <- fx()
	checkEquals( res, 110.0 )
}

