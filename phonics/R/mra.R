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

#' @rdname mra
#' @title Match Rating Approach Encoder
#'
#' @description
#' The Western Airlines matching rating approach name encoder
#'
#' @param word string or vector of strings to encode
#' @param x MRA-encoded character vector
#' @param y MRA-encoded character vector
#'
#' @details
#'
#' The variable \code{word} is the name to be encoded.  The variable
#' \code{maxCodeLen} is \emph{not} supported in this algorithm encoder
#' because the algorithm itself is dependent upon its six-character
#' length.  The variables \code{x} and \code{y} are MRA-encoded and are
#' compared to each other using the MRA comparison specification.
#'
#' @return The \code{mra_encode} function returns match rating approach
#' encoded character vector.  The \code{mra_compare} returns a boolean
#' vector which is \code{TRUE} if \code{x} and \code{y} pass the MRA
#' comparison test.
#'
#' @references
#'
#' G.B. Moore, J.L. Kuhns, J.L. Treffzs, and C.A. Montgomery,
#' \emph{Accessing Individual Records from Personal Data Files Using
#' Nonunique Identifiers,} US National Institute of Standards and
#' Technology, SP-500-2 (1977), p. 17.
#'
#' @family phonics
#'
#' @examples
#' mra_encode("William")
#' mra_encode(c("Peter", "Peady"))
#' mra_encode("Stevenson")

#' @rdname mra
#' @name mra_encode
#' @export
mra_encode <- function(word) {

    ## First, remove any nonalphabetical characters and uppercase it
    word <- gsub("[^[:alpha:]]*", "", word)
    word <- toupper(word)

    ## First character of key = first character of name
    first <- substr(word, 1, 1)
    word <- substr(word, 2, nchar(word))

    ## Delete vowels not at the start of the word
    word <- gsub("[AEIOU]", "", word)
    word <- paste(first, word, sep = "")

    ## Remove duplicate consecutive characters
    word <- gsub("([A-Z])\\1+", "\\1", word)

    ## If longer than 6 characters, take first and last 3...and we have
    ## to vectorize it
    for(i in 1:length(word)) {
        if((l = nchar(word[i])) > 6) {
            first <- substr(word[i], 1, 3)
            last <- substr(word[i], l - 2, l)
            word[i] <- paste(first, last, sep = "");
        }
    }

    return(word)
}

#' @rdname mra
#' @name mra_compare
#' @export
mra_compare <- function(x, y) {
    mra <- data.frame(x = x, y = y, sim = 0, min = 100, stringsAsFactors = FALSE)

    ## Obtain the minimum rating value by calculating the length sum of
    ## the encoded strings and using table A (from Wikipedia).  We start
    ## by setting the minimum to be the sum and move from there.
    mra$lensum <- nchar(mra$x) + nchar(mra$y)
    mra$min[mra$lensum == 12] <- 2
    mra$min[mra$lensum > 7 && mra$lensum <= 11] <- 3
    mra$min[mra$lensum > 4 && mra$lensum <= 7] <- 4
    mra$min[mra$lensum <= 4] <- 5

    ## If the length difference between the encoded strings is 3 or
    ## greater, then no similarity comparison is done.  For us, we
    ## continue the similarity comparison out of laziness and ensure the
    ## minimum is impossibly high to meet.
    mra$min[abs(nchar(mra$x) - nchar(mra$y)) >= 3] <- 100

    ## Start the comparison.
    x <- strsplit(mra$x, split = "")
    y <- strsplit(mra$y, split = "")
    rows <- nrow(mra)
    for(i in 1:rows) {
        ## Process the encoded strings from left to right and remove any
        ## identical characters found from both strings respectively.
        j <- 1
        while(j < min(length(x[[i]]), length(y[[i]]))) {
            if(x[[i]][j] == y[[i]][j]) {
                x[[i]] <- x[[i]][-j]
                y[[i]] <- y[[i]][-j]
            } else
                j <- j + 1
        }

        ## Process the unmatched characters from right to left and
        ## remove any identical characters found from both names
        ## respectively.
        x[[i]] <- rev(x[[i]])
        y[[i]] <- rev(y[[i]])
        j <- 1
        while(j < min(length(x[[i]]), length(y[[i]]))) {
            if(x[[i]][j] == y[[i]][j]) {
                x[[i]] <- x[[i]][-j]
                y[[i]] <- y[[i]][-j]
            } else
                j <- j + 1
        }
        ## Subtract the number of unmatched characters from 6 in the
        ## longer string. This is the similarity rating.
        len <- min(length(x[[i]]), length(y[[i]]))
        mra$sim[i] <- 6 - len
    }

    ## If the similarity is greater than or equal to the minimum
    ## required, it is a successful match.
    mra$match <- (mra$sim >= mra$min)
    return(mra$match)
}
