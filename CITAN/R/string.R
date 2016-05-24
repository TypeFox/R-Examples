## This file is part of the CITAN package for R
##
## Copyright 2011-2015 Marek Gagolewski
##
##
## CITAN is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## CITAN is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
## GNU Lesser General Public License for more details.
##
## You should have received a copy of the GNU Lesser General Public License
## along with CITAN. If not, see <http://www.gnu.org/licenses/>.


# Converts a numeric vector to its character representation or returns \code{"NULL"}.
#
# @title Get nullable numeric values for use in an SQL query
# @param value numeric vector to be processed.
# @return The function returns a character vector containing a
# textual representations of each numeric value or \code{"NULL"}.
# @export
# @seealso \code{\link{sqlSwitchOrNULL}}, \code{\link{sqlStringOrNULL}}
sqlNumericOrNULL <- function(value)
{
	value <- as.character(as.numeric(value));
	ifelse(is.na(value) | is.null(value), "NULL", value);
}



# Returns trimmed and single-quote-escaped (see \code{\link{sqlEscapeTrim}})
# character strings or \code{"NULL"}s on empty input.
#
# @title Get nullable character string values for use in an SQL query
# @param value character vector to be processed.
# @return The function returns a character vector containing a
# processed input values and/or \code{"NULL"}s.
# @export
# @seealso \code{\link{sqlNumericOrNULL}}, \code{\link{sqlSwitchOrNULL}}, \code{\link{sqlEscapeTrim}}
sqlStringOrNULL <- function(value)
{
	value <- sqlEscapeTrim(value);
	ifelse(is.na(value) | is.null(value) | nchar(value)==0,
		"NULL", paste("'", value, "'", sep=""));
}


# Given a character vector, the function
# tries to match each of its values up in the array \code{search}.
# On success, corresponding value of \code{replace} will is returned.
# Otherwise, it outputs \code{"NULL"}.
#
# @title Switch-like nullable construct for use in an SQL query
# @param value character vector to be processed.
# @param search  character vector of length \code{n}.
# @param replace character vector of length \code{n}.
# @return The function returns a character vector containing
# values from \code{replace} that correspond to elements of \code{search} and/or \code{"NULL"}s.
# @export
# @seealso \code{\link{sqlNumericOrNULL}}, \code{\link{sqlStringOrNULL}}, \code{\link{sqlTrim}}
sqlSwitchOrNULL <- function(value, search, replace)
{
	stopifnot(length(search) == length(replace) && length(search) > 0);

	value <- sqlEscapeTrim(value);

	h <- hash(keys=search, values=replace);
	out <- sapply(value, function(x)
	{
		if (is.na(x) || nchar(x) == 0) return("NULL");

		ret <- h[[x]];
		ifelse(is.null(ret), "NULL", ret);
	});

	clear(h);
	rm(h);

	return(out);
}



# Escapes character strings for use in an SQL query.
#
# The SQL standard specifies that single-quotes in strings should be
# escaped by putting two single quotes in a row.
# This function repeats the quotes using \code{\link{stri_replace_all_fixed}}.
#
# @title Escape character strings for use in an SQL query
# @param str a character vector where matches are sought, or an object which can be coerced by \code{as.character} to a character vector.
# @return See 'Value' for \code{\link{stri_replace_all_fixed}}.
# @export
# @seealso \code{\link{sqlTrim}}, \code{\link{sqlEscapeTrim}}
sqlEscape <- function(str, useBytes=FALSE)
{
	stri_replace_all_fixed(str, "'", "''")
}


# Trims white-spaces on both sides of a string.
#
# This function uses \code{\link{stri_replace_all_fixed}}.
#
# @title Trim white-spaces on both sides of character strings
# @param str a character vector where matches are sought, or an object which can be coerced by \code{as.character} to a character vector.
# @return See 'Value' for \code{\link{stri_trim_both}}.
# @export
# @seealso \code{\link{stri_replace_all_fixed}}, \code{\link{sqlEscape}}, \code{\link{sqlEscapeTrim}}
sqlTrim <- function(str)
{
	stri_trim_both(str)
}


# Escapes given character strings for use in an SQL query and trims white-spaces on both sides.
#
# The SQL standard specifies that single-quotes in strings should be
# escaped by putting two single quotes in a row.
# This function repeats the quotes using \code{\link{gsub}}.
#
# @title Escape character strings for use in an SQL query and trim white-spaces on both sides
# @param str a character vector where matches are sought, or an object which can be coerced by \code{as.character} to a character vector.
# @return character vector
# @export
# @seealso \code{\link{sqlEscape}}, \code{\link{sqlTrim}}
sqlEscapeTrim <- function(str)
{
	sqlEscape(sqlTrim(str))
}
