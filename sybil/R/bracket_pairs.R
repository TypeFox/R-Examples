#  bracket_pairs.R
#  FBA and friends with R.
#
#  Copyright (C) 2010-2014 Gabriel Gelius-Dietrich, Dpt. for Bioinformatics,
#  Institute for Informatics, Heinrich-Heine-University, Duesseldorf, Germany.
#  All right reserved.
#  Email: geliudie@uni-duesseldorf.de
#
#  This file is part of sybil.
#
#  Sybil is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  Sybil is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with sybil.  If not, see <http://www.gnu.org/licenses/>.


################################################
# Function: .bracket_pairs
#
#
#


# ---------------------------------------------------------------------------- #

# This function finds corresponding pairs of brackets in a logical rule.
# It returns a matrix with three columns. Each row containing one bracket
# pair, the first value is the position of the opening bracket in the
# string, the second value is the position of the closing bracket and the
# third element is the length of the text between the brackets (excluding
# the brackets).
# Argument rule is a character vector, each element being one letter of
# the original text string.

.bracket_pairs <- function(rule) {

    if (!any(is.na(rule))) {

        # stack for opening brackets
        #s_open <- stapel()
        st <- "s_open"
        stinit(st)

        pairc <- 0   # counter for pairs

        # positions of opening brackets
        nop   <- length(grep("(", rule, fixed = TRUE))

        # result matrix, nop is the number of bracket pairs
        bpair <- matrix(0, nrow = nop, ncol = 3)


        # push each the position of each '(' onto the stack
        for (i in seq(along = rule)) {
            if (rule[i] == "(") {
                #push(s_open, i)
                stpush(st, i)
            }
            # if we find a ')' the '(' on top of the stack is the
            # corrseponding bracket.
            else if (rule[i] == ")") {
                pairc <- pairc + 1
                #ob <- as.integer(pop(s_open))
                ob <- as.integer(stpop(st))
                bpair[pairc,] <- c(ob, i, i-ob-1)
            }
        }

        colnames(bpair) <- c("open", "close", "length")

        # sort the pairs according to their length
        bpair <- bpair[order(bpair[, "length"]), , drop = FALSE]


    }
    else {
        bpair <- NA
    }

    return(bpair)

}
