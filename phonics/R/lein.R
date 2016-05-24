## Copyright (c) 2015, James P. Howard, II <jh@jameshoward.us>
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

#' @title Lein Name Coding
#'
#' @description
#' The Lein name coding procedure.
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
#' @return the Lein encoded character vector
#'
#' @references
#'
#' Billy T. Lynch and William L. Arends. "Selection of surname coding
#' procedure for the SRS record linkage system." United States
#' Department of Agriculture, Sample Survey Research Branch, Research
#' Division, Washington, 1977.
#'
#' @family phonics
#'
#' @examples
#' lein("William")
#' lein(c("Peter", "Peady"))
#' lein("Stevenson", maxCodeLen = 8)
#'
#' @export
lein <- function(word, maxCodeLen = 4) {

    ## First, remove any nonalphabetical characters and uppercase it
    word <- gsub("[^[:alpha:]]*", "", word)
    word <- toupper(word)

    ## First character of key = first character of name
    first <- substr(word, 1, 1)
    word <- substr(word, 2, nchar(word))

    ## Delete vowels and Y, W, and H
    word <- gsub("A|E|I|O|U|Y|W|H", "", word)

    ## Append word except for first character to first
    word <- paste(first, word, sep = "")
    
    ## Remove duplicate consecutive characters
    word <- gsub("([A-Z])\\1+", "\\1", word)
    
    ## Grab the last n - 1 characters again, for recoding
    word <- substr(word, 2, nchar(word))

    ## D, T -> 1
    word <- gsub("D|T", "1", word)

    ## M, N -> 2
    word <- gsub("M|N", "2", word)

    ## L, R -> 3
    word <- gsub("L|R", "3", word)

    ## B, F, P, V -> 4
    word <- gsub("B|F|P|V", "4", word)

    ## C, J, K, G, Q, S, X, Z -> 5
    word <- gsub("C|J|K|G|Q|S|X|Z", "5", word)

    ## Zero-pad and truncate to requested length
    word <- gsub("$", paste(rep(0, maxCodeLen), collapse = ""), word)
    word <- substr(word, 1, maxCodeLen)

    return(word)
}
