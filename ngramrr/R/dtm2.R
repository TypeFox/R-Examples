#' @rdname dtmwrappers
#' @name dtmwrappers
#' @title Wrappers to DocumentTermMatrix and DocumentTermMatrix to use n-gram tokenizaion
#' @description
#' Wrappers to \code{DocumentTermMatrix} and \code{DocumentTermMatrix} to use n-gram tokenization provided by \code{ngramrr}.
#' @param x character vector, \code{Source} or \code{Corpus} to be converted
#' @param char logical, using character n-gram. char = FALSE denotes word n-gram.
#' @param ngmin integer, minimun order of n-gram
#' @param ngmax integer, maximun order of n-gram
#' @param rmEOL logical, remove ngrams wih EOL character
#' @param ... Additional options for \code{DocumentTermMatrix} or \code{DocumentTermMatrix}
#' @return \code{DocumentTermMatrix} or \code{DocumentTermMatrix}
#' @seealso \code{ngramrr}, \code{\link[tm]{DocumentTermMatrix}}, \code{\link[tm]{TermDocumentMatrix}}
#' @examples
#' nirvana <- c("hello hello hello how low", "hello hello hello how low",
#' "hello hello hello how low", "hello hello hello",
#' "with the lights out", "it's less dangerous", "here we are now", "entertain us",
#' "i feel stupid", "and contagious", "here we are now", "entertain us",
#' "a mulatto", "an albino", "a mosquito", "my libido", "yeah", "hey yay")
#' dtm2(nirvana, ngmax = 3, removePunctuation = TRUE)
NULL

#' @rdname dtmwrappers
#' @export
dtm2 <- function(x, char = FALSE, ngmin = 1, ngmax = 2, rmEOL = TRUE, ...) {
    return(DocumentTermMatrix(conv_corpus(x), control = list(tokenize = function(z) ngramrr(z, char = char, ngmin = ngmin, ngmax = ngmax, rmEOL = rmEOL), ...)))
}

#' @rdname dtmwrappers
#' @export
tdm2 <- function(x, char = FALSE, ngmin = 1, ngmax = 2, rmEOL = TRUE, ...) {
    return(TermDocumentMatrix(conv_corpus(x), control = list(tokenize = function(z) ngramrr(z, char = char, ngmin = ngmin, ngmax = ngmax, rmEOL = rmEOL), ...)))
}

## helper function, convert everything to corpus
conv_corpus <- function(x) {
    if ("character" %in% class(x)) {
        return(Corpus(VectorSource(x)))
    } else if ("Source" %in% class(x)) {
        return(Corpus(x))
    } else if ("Corpus" %in% class(x)) {
        return(x)
    } else {
        stop("x must be character vector, Source or Corpus.")
    }
}
