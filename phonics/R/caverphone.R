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

#' @title Caverphone
#'
#' @description
#' The Caverphone family of phonetic algorithms
#'
#' @param word string or vector of strings to encode
#' @param maxCodeLen   maximum length of the resulting encodings, in characters
#' @param modified     if \code{TRUE}, use the Caverphone 2 algorithm
#'
#' @details
#'
#' The variable \code{maxCodeLen} is the limit on how long the returned
#' Caverphone code should be.  The default is 6, unless \code{modified}
#' is set to \code{TRUE}, then the default is 10.
#'
#' The variable \code{modified} directs \code{caverphone} to use the
#' Caverphone2 method, instead of the original.
#'
#' @return the Caverphone encoded character vector
#'
#' @references
#'
#' David Hood, "Caverphone: Phonetic matching algorithm," Technical
#' Paper CTP060902, University of Otago, New Zealand, 2002.
#'
#' David Hood, "Caverphone Revisited," Technical Paper CTP150804
#' University of Otago, New Zealand, 2004.
#'
#' @family phonics
#'
#' @examples
#' caverphone("William")
#' caverphone(c("Peter", "Peady"), modified = TRUE)
#' caverphone("Stevenson", maxCodeLen = 4)
#'
#' @export
caverphone <- function(word, maxCodeLen = NULL, modified = FALSE) {
    ## From here on, this is a line-for-line translation of the Apache
    ## Commons Caverphone and Caverphone2 implementations, which both
    ## used regular expressions for substantially all of the work.

    ## Set the maxCodeLen if not set
    if(is.null(maxCodeLen))
        if(modified == TRUE)
            maxCodeLen <- 10
        else
            maxCodeLen <- 6
    
    ## First, remove any nonalphabetical characters and lowercase it
    word <- gsub("[^[:alpha:]]*", "", word)
    word <- tolower(word)
    
    if(modified == TRUE)
        word <- caverphone_modified(word)
    else
        word <- caverphone_original(word)

    ## Pad the wording with maxCodeLen 1s and truncate
    word <- gsub("$", paste(rep(1, maxCodeLen), collapse = ""), word)
    word <- substr(word, 1, maxCodeLen)
    return(word)
}

caverphone_original <- function(word) {

    ## Starting and ending special cases
    word <- gsub("^cough", "cou2f", word)
    word <- gsub("^rough", "rou2f", word)
    word <- gsub("^tough", "tou2f", word)
    word <- gsub("^enough", "enou2f", word)
    word <- gsub("^gn", "2n", word)
    word <- gsub("mb$", "m2", word)

    ## Core encodings
    word <- gsub("cq", "2q", word)
    word <- gsub("ci", "si", word)
    word <- gsub("ce", "se", word)
    word <- gsub("cy", "sy", word)
    word <- gsub("tch", "2ch", word)
    word <- gsub("c", "k", word)
    word <- gsub("q", "k", word)
    word <- gsub("x", "k", word)
    word <- gsub("v", "f", word)
    word <- gsub("dg", "2g", word)
    word <- gsub("tio", "sio", word)
    word <- gsub("tia", "sia", word)
    word <- gsub("d", "t", word)
    word <- gsub("ph", "fh", word)
    word <- gsub("b", "p", word)
    word <- gsub("sh", "s2", word)
    word <- gsub("z", "s", word)
    word <- gsub("^[aeiou]", "A", word)
    word <- gsub("[aeiou]", "3", word)
    word <- gsub("3gh3", "3kh3", word)
    word <- gsub("gh", "22", word)
    word <- gsub("g", "k", word)
    word <- gsub("ss*", "S", word)
    word <- gsub("tt*", "T", word)
    word <- gsub("pp*", "P", word)
    word <- gsub("kk*", "K", word)
    word <- gsub("ff*", "F", word)
    word <- gsub("mm*", "M", word)
    word <- gsub("nn*", "N", word)
    word <- gsub("w3", "W3", word)
    word <- gsub("wy", "Wy", word)
    word <- gsub("wh3", "Wh3", word)
    word <- gsub("why", "Why", word)
    word <- gsub("w", "2", word)
    word <- gsub("^h", "A", word)
    word <- gsub("h", "2", word)
    word <- gsub("r3", "R3", word)
    word <- gsub("ry", "Ry", word)
    word <- gsub("r", "2", word)
    word <- gsub("l3", "L3", word)
    word <- gsub("ly", "Ly", word)
    word <- gsub("l", "2", word)
    word <- gsub("j", "y", word)
    word <- gsub("y3", "Y3", word)
    word <- gsub("y", "2", word)
    word <- gsub("[23]", "", word)

    return(word)
}

caverphone_modified <- function(word) {

    ## Starting and ending special cases
    word <- gsub("e$", "", word)
    word <- gsub("^cough", "cou2f", word)
    word <- gsub("^rough", "rou2f", word)
    word <- gsub("^tough", "tou2f", word)
    word <- gsub("^enough", "enou2f", word)
    word <- gsub("^trough", "trou2f", word)
    word <- gsub("^gn", "2n", word)
    word <- gsub("mb$", "m2", word)

    ## Core encodings
    word <- gsub("cq", "2q", word)
    word <- gsub("ci", "si", word)
    word <- gsub("ce", "se", word)
    word <- gsub("cy", "sy", word)
    word <- gsub("tch", "2ch", word)
    word <- gsub("c", "k", word)
    word <- gsub("q", "k", word)
    word <- gsub("x", "k", word)
    word <- gsub("v", "f", word)
    word <- gsub("dg", "2g", word)
    word <- gsub("tio", "sio", word)
    word <- gsub("tia", "sia", word)
    word <- gsub("d", "t", word)
    word <- gsub("ph", "fh", word)
    word <- gsub("b", "p", word)
    word <- gsub("sh", "s2", word)
    word <- gsub("z", "s", word)
    word <- gsub("^[aeiou]", "A", word)
    word <- gsub("[aeiou]", "3", word)
    word <- gsub("j", "y", word)
    word <- gsub("^y3", "Y3", word)
    word <- gsub("^y", "a", word)
    word <- gsub("y", "3", word)
    word <- gsub("3gh3", "3kh3", word)
    word <- gsub("gh", "22", word)
    word <- gsub("g", "k", word)
    word <- gsub("ss*", "S", word)
    word <- gsub("tt*", "T", word)
    word <- gsub("pp*", "P", word)
    word <- gsub("kk*", "K", word)
    word <- gsub("ff*", "F", word)
    word <- gsub("mm*", "M", word)
    word <- gsub("nn*", "N", word)
    word <- gsub("w3", "W3", word)
    word <- gsub("wh3", "Wh3", word)
    word <- gsub("w$", "3", word)
    word <- gsub("w", "2", word)
    word <- gsub("^h", "A", word)
    word <- gsub("h", "2", word)
    word <- gsub("r3", "R3", word)
    word <- gsub("r$", "3", word)
    word <- gsub("r", "2", word)
    word <- gsub("l3", "L3", word)
    word <- gsub("l$", "3", word)
    word <- gsub("l", "2", word)
    word <- gsub("2", "", word)
    word <- gsub("3$", "A", word)
    word <- gsub("3", "", word)

    return(word)
}
