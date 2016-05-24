## Copyright (c) 2016, James P. Howard, II <jh@jameshoward.us>
##
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are
## met:
##
##     Redistributions of source code must retain the above copyright
##     notice, this list of conditions and the following disclaimer.
##
##     Redistributions in binary form must reproduce the above copyright
##     notice, this list of conditions and the following disclaimer in
##     the documentation and/or other materials provided with the
##     distribution.
##
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
## "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
## LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
## A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
## HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
## SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
## LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
## DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
## THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
## OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#' @title Cologne Phonetic Name Coding
#'
#' @description
#' The Cologne phonetic name coding procedure.
#'
#' @param word string or vector of strings to encode
#' @param maxCodeLen   maximum length of the resulting encodings, in characters
#'
#' @details
#'
#' The variable \code{word} is the name to be encoded.  The variable
#' \code{maxCodeLen} is the limit on how long the returned name code
#' should be.  The default is 4.
#'
#' @return the Cologne encoded character vector
#'
#' @references
#'
#' Hans Joachim Postel. "Die Koelner Phonetik. Ein Verfahren zur
#' Identifizierung von Personennamen auf der Grundlage der
#' Gestaltanalyse."  \emph{IBM-Nachrichten} 19. Jahrgang, 1969,
#' p. 925-931.
#'
#' @family phonics
#'
#' @examples
#' lein("William")
#' lein(c("Peter", "Peady"))
#' lein("Stevenson", maxCodeLen = 8)
#'
#' @export
cologne <- function(word, maxCodeLen = NULL) {

    ## First, remove any nonalphabetical characters and uppercase it
    word <- gsub("[^[:alpha:]]*", "", word)
    word <- toupper(word)

    ## Remove umlauts and eszett
    word <- gsub("\u00C4", "A", word)
    word <- gsub("\u00DC", "U", word)
    word <- gsub("\u00D6", "O", word)
    word <- gsub("\u00DF", "S", word)

    ## Work through the rules...but backwards, mostly, here's 8s
    word <- gsub("([CKQ])X", "\\18", word)
    word <- gsub("[DT]([CSZ])", "8\\1", word)
    word <- gsub("([SZ])C", "\\18", word)
    word <- gsub("^C([^AHKLOQRUX])", "8\\1", word)
    word <- gsub("C([^AHKOQUX])", "8\\1", word)
    word <- gsub("[SZ]", "8", word)

    ## Rule #7
    word <- gsub("R", "7", word)

    ## Rule #6
    word <- gsub("[MN]", "6", word)

    ## Rule #5
    word <- gsub("L", "5", word)

    ## Rule #48
    word <- gsub("X", "48", word)

    ## Rule #4
    word <- gsub("[CGKQ]", "4", word)

    ## Rule #3
    ## And we can strip the H since it will not be coded
    word <- gsub("PH|[FVW]", "3", word)

    ## Rule #2
    word <- gsub("[DT]", "2", word)

    ## Rule #1
    word <- gsub("[BP]", "1", word)

    ## Rule #H
    word <- gsub("H", "", word)

    ## Rule #0
    word <- gsub("[AEIJOUY]", "0", word)

    ## Remove duplicate consecutive characters
    word <- gsub("([0-9])\\1+", "\\1\\2", word)

    ## Remove all 0s, except first
    first <- substr(word, 1, 1)
    word <- substr(word, 2, nchar(word))
    word <- gsub("0", "", word)
    word <- paste(first, word, sep = "")

    return(word)
}
