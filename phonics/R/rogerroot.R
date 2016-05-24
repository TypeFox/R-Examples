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

#' @title Roger Root Name Coding Procedure
#'
#' @description
#' Provides the Roger Root name coding system
#'
#' @param word string or vector of strings to encode
#' @param maxCodeLen  maximum length of the resulting encodings, in characters
#'
#' @details
#'
#' The \code{rogerroot} function phentically encodes the given string
#' using the Roger Root algorithm.  The variable \code{word} is a string
#' or vector of strings to encode.
#'
#' The variable \code{maxCodeLen} is the limit on how long the returned
#' code should be.  The default is 5.
#'
#' @return the Roger Root encoded character vector
#'
#' @references
#'
#' Robert L. Taft, \emph{Name search techniques}, Bureau of Systems
#' Development, Albany, New York, 1970.
#'
#' @family phonics
#'
#' @examples
#' rogerroot("William")
#' rogerroot(c("Peter", "Peady"))
#' rogerroot("Stevenson")
#'
#' @importFrom utils read.csv
#'
#' @export
rogerroot <- function(word, maxCodeLen = 5) {

    ## First, remove any nonalphabetical characters and uppercase it
    word <- gsub("[^[:alpha:]]*", "", word)
    word <- toupper(word)

    ## First letter table...these are write-once tables...
    letterTable<-"letter,code\n^A,1\n^B,09\n^CE,00\n^CH,06\n^CI,00\n^CY,00\n^C,07\n^DG,07\n^D,01\n^E,1\n^F,08\n^GF,08\n^GM,03\n^GN,02\n^G,07\n^H,2\n^I,1\n^J,3\n^KN,02\n^K,07\n^L,05\n^M,03\n^N,02\n^O,1\n^PF,08\n^PH,08\n^PN,02\n^P,09\n^Q,07\n^R,04\n^SCH,06\n^SH,06\n^S,00\n^TSCH,06\n^TSH,06\n^TS,00\n^T,01\n^U,1\n^V,08\n^WR,04\n^W,4\n^X,07\n^Y,5\n^Z,00\n"
    letters <- read.csv(colClasses=c("character", "character"), text = letterTable)
    for(i in 1:nrow(letters))
        word <- gsub(letters$letter[i], letters$code[i], word)

    ## Basic letter table
    letterTable<-"letter,code\nB,9\nCE,0\nCH,6\nCI,0\nCY,0\nC,7\nDG,7\nD,1\nF,8\nG,7\nJ,6\nK,7\nL,5\nM,3\nN,2\nPH,8\nP,8\nQ,7\nR,4\nSCH,6\nSH,6\nS,0\nTSCH,6\nTSH,6\nTS,0\nT,1\nV,8\nX,7\nZ,0"
    letters <- read.csv(colClasses=c("character", "character"), text = letterTable)
    for(i in 1:nrow(letters))
        word <- gsub(letters$letter[i], letters$code[i], word)

    ## Remove duplicate consecutive characters
    word <- gsub("([1-9])\\1+", "\\1", word)
    word <- gsub(".([0])\\1+", "\\1", word)

    ## Remove non-numeric characters
    word <- gsub("[^0-9]", "", word)

    ## Truncate to requested length
    zeros <- paste(rep(0, maxCodeLen), sep = "", collapse = "")
    word <- gsub("$", zeros, word)
    word <- substr(word, 1, maxCodeLen)

    return(word)
}
