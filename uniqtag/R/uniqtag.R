#' Abbreviate strings to short, unique identifiers.
#' @docType package
#' @name uniqtag-package
#' @author Shaun Jackman \email{sjackman@@gmail.com}
NULL

#' Return the k-mers of a string.
#'
#' Return the k-mers (substrings of size \code{k}) of the string \code{x}, or
#' return the string \code{x} itself if it is shorter than k.
#' @describeIn kmers_of Return the k-mers of the string \code{x}.
#' @param k the size of the substrings, an integer
#' @param x a character string
#' @return kmers_of: a character vector of the k-mers of \code{x}
#' @export
kmers_of <- function(x, k)
	if (nchar(x) < k) x else
		substring(x, 1:(nchar(x) - k + 1), k:nchar(x))

#' @describeIn kmers_of Return the k-mers of the strings \code{xs}.
#' @param xs a character vector
#' @return vkmers_of: a list of character vectors of the k-mers of \code{xs}
#' @export
vkmers_of <- function(xs, k)
	Vectorize(kmers_of, SIMPLIFY = FALSE)(xs, k)

#' Cumulative count of strings.
#'
#' Return an integer vector counting the number of occurrences of each string up to that position in the vector.
#' @param xs a character vector
#' @return an integer vector of the cumulative string counts
#' @examples
#' cumcount(abbreviate(state.name, 3, strict = TRUE))
#' @export
cumcount <- function(xs) {
	counts <- new.env(parent = emptyenv())
	setNames(vapply(xs, function(x)
		counts[[x]] <- 1L + mget(x, counts, ifnotfound = 0L)[[1]],
		integer(1)), xs)
}

#' Make character strings unique.
#'
#' Apppend sequence numbers to duplicate elements to make all elements of a character vector unique.
#' @param xs a character vector
#' @param sep a character string used to separate a duplicate string from its sequence number
#' @describeIn make_unique Append a sequence number to duplicated elements, including the first occurence.
#' @seealso make.unique
#' @examples
#' abcb <- c("a", "b", "c", "b")
#' make_unique(abcb)
#' make_unique_duplicates(abcb)
#' make_unique_all(abcb)
#' make_unique_all_or_none(abcb)
#' make_unique_all_or_none(c("a", "b", "c"))
#' x <- make_unique(abbreviate(state.name, 3, strict = TRUE))
#' x[grep("-", x)]
#' @export
make_unique <- function(xs, sep = '-') {
	i <- xs %in% xs[duplicated(xs)]
	xs[i] <- make_unique_all(xs[i], sep)
	xs
}

#' @describeIn make_unique Append a sequence number to duplicated elements, except the first occurence.
#'
#' This function behaves similarly to make.unique
#' @export
make_unique_duplicates <- function(xs, sep = '-') {
	i <- duplicated(xs)
	xs[i] <- make_unique_all(xs[i], sep)
	xs
}

#' @describeIn make_unique Append a sequence number to every element.
#' @export
make_unique_all <- function(xs, sep = "-") {
	xs[] <- paste(xs, cumcount(xs), sep = sep)
	xs
}

#' @describeIn make_unique Append a sequence number to every element or no elements.
#'
#' Return \code{xs} unchanged if the elements of the character vector \code{xs} are already unique.
#' Otherwise append a sequence number to every element.
#' @export
make_unique_all_or_none <- function(xs, sep = '-')
	if (anyDuplicated(xs)) make_unique_all(xs, sep) else xs

#' Abbreviate strings to short, unique identifiers.
#'
#' Abbreviate strings to unique substrings of \code{k} characters.
#'
#' For each string in a set of strings, determine a unique tag that is a substring of fixed size \code{k} unique to that string, if it has one. If no such unique substring exists, the least frequent substring is used. If multiple unique substrings exist, the lexicographically smallest substring is used. This lexicographically smallest substring of size \code{k} is called the UniqTag of that string.
#'
#' The lexicographically smallest substring depend on the locale's sort order.
#' You may wish to first call \code{Sys.setlocale("LC_COLLATE", "C")}
#'
#' @examples
#' Sys.setlocale("LC_COLLATE", "C")
#' states <- sub(" ", "", state.name)
#' uniqtags <- uniqtag(states)
#' uniqtags4 <- uniqtag(states, k = 4)
#' uniqtags3 <- uniqtag(states, k = 3)
#' uniqtags3x <- uniqtag(states, k = 3, uniq = make_unique)
#' table(nchar(states))
#' table(nchar(uniqtags))
#' table(nchar(uniqtags4))
#' table(nchar(uniqtags3))
#' table(nchar(uniqtags3x))
#' uniqtags3[grep("-", uniqtags3x)]
#' @param xs a character vector
#' @param k the size of the identifier, an integer
#' @param uniq a function to make the abbreviations unique, such as make_unique, make_unique_duplicates, make_unique_all_or_none, make_unique_all, make.unique, or to disable this function, identity or NULL
#' @param sep a character string used to separate a duplicate string from its sequence number
#' @return a character vector of the UniqTags of the strings \code{x}
#' @seealso abbreviate, locales, make.unique
#' @export
uniqtag <- function(xs, k = 9, uniq = make_unique_all_or_none, sep = '-') {
	if (is.null(uniq)) {
		uniq <- identity
		sep <- NA
	}
	counts <- table(unlist(lapply(vkmers_of(xs, k), unique)))
	counts_kmers <- setNames(
		paste0(format(counts, justify = "right"), names(counts)),
		names(counts))
	tags <- vapply(xs, function(x)
		names(counts_kmers)[match(min(counts_kmers[kmers_of(x, k)]), counts_kmers)],
		character(1))
	if (is.na(sep)) uniq(tags) else uniq(tags, sep)
}
