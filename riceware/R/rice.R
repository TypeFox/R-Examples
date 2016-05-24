##' Checks the validity of the token.
##'
##' The token is a 5-digit number representing the outcomes of 5 rolls
##' of a 6-faced dice. The function returns \code{TRUE} if the
##' \code{token} is correctly formatted and \code{FALSE} otherwise.
##' @title Check the Token
##' @param token a character vector of length 1 representing the outcomes of 5 rolls of a dice.
##' @return \code{TRUE} if correctly formatted, \code{FALSE} otherwise.
##' @author Francois Michonneau
##' @export
check_token <- function(token) {
    (length(token) == 1L) && grepl("^[1-6]{5}$", as.character(token))
}

##' Generates the tokens.
##'
##' This function generates as many tokens as needed to create the
##' passphrase. Currently, two methods can be used to generate the
##' random tokens.
##' \enumerate{
##'
##'  \item The simplest and the default (\code{pseudorandom}) uses the
##' function \code{sample} to simulate the dice rolls. Numbers
##' generated this way are not truly random but are a decent
##' approximation.
##'
##'  \item The other option (\code{random}) uses the \code{random} package
##' that gets truly random numbers by converting atmospheric noise
##' into numbers. The main issue is that someone could monitor your
##' network and intercept the numbers that are being used. If you are
##' concerned about this, use a physical dice.
##' }
##'
##' Note that if you want to use the \code{random} method, you will
##' need an internet connection. The service that provides these
##' random numbers (\url{http://www.random.org}) has daily quotas,
##' don't go too crazy (if you need to, you can purchase additional
##' bits, see \url{http://www.random.org/quota}).
##' @title Generate the tokens
##' @param n_words The number of tokens to generate.
##' @param method The method used to draw the random numbers. See below for more details.
##' @return A character vector representing the generated tokens.
##' @author Francois Michonneau
##' @seealso \url{http://www.random.org} the website that generates
##' the true random numbers, and the random package from Dirk
##' Eddelbuettel.
##' @export
generate_token <- function(n_words, method = c("pseudorandom", "random")) {
    method <- match.arg(method)
    if (identical(method,  "pseudorandom")) {
        rolls <- replicate(sample(1:6, 5), n = n_words)
        tok <- apply(rolls, 2, paste0, collapse = "")
    } else if (identical(method, "random")) {
        n_words <- as.integer(n_words)
        rolls <- random::randomNumbers(n = n_words * 5L, min = 1, max = 6, col = 5)
        tok <- apply(rolls, 1, paste0, collapse = "")
    }
    tok
}

##' Retrieves the word corresponding to a given token.
##'
##' Given a token and a list of words, this function returns the word
##' matching the supplied token.
##' @title Match the token to a word
##' @param token A correctly formatted token
##' @param wordlist A \code{data.frame} with at least 2 columns: one
##' named \sQuote{token} that contains the tokens, one named
##' \sQuote{word} that contains the words corresponding to the tokens.
##' @param title_case If \code{TRUE}, the first letter of each word in
##' the passphrase will be capitalized
##' @return A character vector of length 1 representing the word
##' corresponding to the token.
##' @author Francois Michonneau
##' @export
match_token <- function(token, wordlist = riceware::wordlist_en,
                        title_case = TRUE) {
    if (check_token(token)) {
        wrd <- wordlist[wordlist$token %in% token, "word", drop = TRUE]
        if (title_case) {
            wrd <- paste0(toupper(substr(wrd, 1L, 1L)), substring(wrd, 2L))
        }
    } else {
        stop("invalid token")
    }
    wrd
}

##' Generates a passphrase.
##'
##' Given a wordlist and a number of words, this function generates a
##' passphrase. You can control the wordlist you choose and whether
##' the passphrase uses title case by providing additional arguments
##' that will be passed to \code{\link{match_token}}.
##' @title Generates a passphrase
##' @param tokens a vector of character representing the tokens to be
##' used to generate the passphrase. By default, 7 are randomly
##' generated using \code{generate_token(7)}.
##' @param verbose if \code{TRUE} the passphrase is displayed as a message
##' @param ... additional parameters to be passed to \code{\link{match_token}}
##' @return a character string representing the passphrase
##' @author Francois Michonneau
##' @seealso \code{\link{match_token}}, \code{\link{generate_token}}
##' @examples
##'   generate_passphrase(tokens = generate_token(7, "pseudorandom"),
##'                       verbose = FALSE)
##' @export
generate_passphrase <- function(tokens = generate_token(7), verbose = TRUE, ...) {
    pass <- sapply(tokens, match_token, ...)
    if (verbose) message("Your passphrase is: ", paste(pass, collapse = " "))
    paste0(pass, collapse = "")
}
