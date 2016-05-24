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

#' @title New York State Identification and Intelligence System
#'
#' @description
#' The NYSIIS phonetic algorithm
#'
#' @param word string or vector of strings to encode
#' @param maxCodeLen   maximum length of the resulting encodings, in characters
#' @param modified     if \code{TRUE}, use the modified NYSIIS algorithm
#'
#' @details The \code{nysiis} function phentically encodes the given
#' string using the New York State Identification and Intelligence
#' System (NYSIIS) algorithm. The algorithm is based on the
#' implementation provided by Wikipedia and is implemented in pure R
#' using regular expressions.
#'
#' The variable \code{maxCodeLen} is the limit on how long the returned
#' NYSIIS code should be.  The default is 6.
#'
#' The variable \code{modified} directs \code{nysiis} to use the
#' modified method instead of the original.
#'
#' @return the NYSIIS encoded character vector
#'
#' @references
#'
#' Robert L. Taft, \emph{Name search techniques}, Bureau of Systems
#' Development, Albany, New York, 1970.
#'
#' @family phonics
#'
#' @examples
#' nysiis("Robert")
#' nysiis("rupert")
#' nysiis(c("Alabama", "Alaska"), modified = TRUE)
#' nysiis("mississippi", 4)
#'
#' @export
nysiis <- function(word, maxCodeLen = 6, modified = FALSE) {
    ## Both NYSIIS and the modified NYSIIS are based on the
    ## implementation described at
    ## http://www.dropby.com/NYSIISTextStrings.html

    if(modified == TRUE)
        return(nysiis_modified(word, maxCodeLen))
    else
        return(nysiis_original(word, maxCodeLen))
}

nysiis_original <- function(word, maxCodeLen = 6) {

    ## First, remove any nonalphabetical characters and capitalize it
    word <- gsub("[^[:alpha:]]*", "", word)
    word <- toupper(word)

    ## Translate first characters of name: MAC to MCC, KN to N, K to C, PH,
    ## PF to FF, SCH to SSS
    word <- gsub("^MAC", "MCC", word)
    word <- gsub("KN", "NN", word)
    word <- gsub("K", "C", word)
    word <- gsub("^PF", "FF", word)
    word <- gsub("PH", "FF", word)
    word <- gsub("SCH", "SSS", word)

    ## Translate last characters of name: EE to Y, IE to Y, DT, RT, RD,
    ## NT, ND to D
    word <- gsub("EE$", "Y", word)
    word <- gsub("IE$", "Y", word)
    word <- gsub("DT$", "D", word)
    word <- gsub("RT$", "D", word)
    word <- gsub("RD$", "D", word)
    word <- gsub("NT$", "D", word)
    word <- gsub("ND$", "D", word)

    ## First character of key = first character of name.
    first <- substr(word, 1, 1)
    word <- substr(word, 2, nchar(word))

    ## EV to AF else A, E, I, O, U to A
    word <- gsub("EV", "AF", word)
    word <- gsub("E|I|O|U", "A", word)

    ## Q to G, Z to S, M to N
    word <- gsub("Q", "G", word)
    word <- gsub("Z", "S", word)
    word <- gsub("M", "N", word)

    ## KN to N else K to C
    ## SCH to SSS, PH to FF
    ## Rules are implemented as part of opening block

    ## H to If previous or next is non-vowel, previous.
    word <- gsub("([^AEIOU])H", "\\1", word)
    word <- gsub("(.)H[^AEIOU]", "\\1", word)

    ## W to If previous is vowel, A
    word <- gsub("([AEIOU])W", "A", word)

    ## If last character is S, remove it
    word <- gsub("S$", "", word)

    ## If last characters are AY, replace with Y
    word <- gsub("AY$", "Y", word)

    ## Remove duplicate consecutive characters
    word <- gsub("([A-Z])\\1+", "\\1", word)

    ## If last character is A, remove it
    word <- gsub("A$", "", word)

    ## Append word except for first character to first
    word <- paste(first, word, sep = "")

    ## Truncate to requested length
    word <- substr(word, 1, maxCodeLen)

    return(word)
}

nysiis_modified <- function(word, maxCodeLen = 6) {

    ## First, remove any nonalphabetical characters and capitalize it
    word <- gsub("[^[:alpha:]]*", "", word)
    word <- toupper(word)

    ## Translate first characters of name: MAC to MC, PF to FF
    word <- gsub("^MAC", "MC", word)
    word <- gsub("^PF", "FF", word)

    ## First character of key = first character of name.
    first <- substr(word, 1, 1)

    ## Remove a trailing S
    word <- gsub("S$|Z$", "", word)

    ## Translate last characters of name
    word <- gsub("IX$", "IC", word)
    word <- gsub("EX$", "EC", word)
    word <- gsub("YE$|EE$|IE$", "Y", word)
    word <- gsub("DT$|RT$|RD$|NT$|ND$", "D", word)

    ## transcode 'EV' to 'EF' if not at start of name
    word <- gsub("(.+)EV", "\1EF", word)

    ## EV to AF else A, E, I, O, U to A
    word <- gsub("E|I|O|U", "A", word)

    ## W to If previous is vowel, A
    word <- gsub("([AEIOU])W", "A", word)

    ## transcode 'GHT' to 'GT'
    word <- gsub("GHT", "GT", word)

    ## transcode 'DG' to 'G'
    word <- gsub("DG", "G", word)

    ## Q to G,  M to N, SH to S, SCH to S, YW to Y, Y to A,
    word <- gsub("M", "N", word)
    word <- gsub("Q", "G", word)
    word <- gsub("(.+)SH", "\\1S", word)
    word <- gsub("(.+)SCH", "\\1S", word)
    word <- gsub("YW", "Y", word)

    ## if not first or last character, change 'Y' to 'A'
    last <- substring(word, nchar(word), nchar(word))
    word <- substring(word, 1, nchar(word) - 1)
    word <- gsub("Y", "A", word)
    word <- paste(word, last, sep = "")

    ## WR to R, Z to S
    word <- gsub("WR", "R", word)
    word <- gsub("Z", "S", word)

    ## Remove duplicate consecutive characters
    word <- gsub("([A-Z])\\1+", "\\1", word)

    ## If last character is A, remove it
    word <- gsub("A$", "", word)

    ## Append word except for first character to first
    word <- substr(word, 2, nchar(word))
    word <- paste(first, word, sep = "")

    ## transcode 'PH' to 'F'
    word <- gsub("PH", "F", word)

    ## change 'KN' to 'N', else 'K' to 'C'
    word <- gsub("KN", "N", word)
    word <- gsub("K", "C", word)

    ## H to If previous or next is non-vowel, previous.
    word <- gsub("([^AEIOU])H", "\\1", word)
    word <- gsub("(.)H[^AEIOU]", "\\1", word)

    ## transcode terminal 'AY' to 'Y'
    word <- gsub("AY$", "Y", word)

    ## Remove duplicate consecutive characters
    word <- gsub("([A-Z])\\1+", "\\1", word)

    ## Truncate to requested length
    word <- substr(word, 1, maxCodeLen)

    return(word)
}
