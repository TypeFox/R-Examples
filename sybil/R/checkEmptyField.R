#  checkEmptyField.R
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
# Function: .checkEmptyField
#
#
#

.checkEmptyField <- function(vec, name) {

    nalines <- which(is.na(vec) | nchar(vec, allowNA = TRUE) == 0)
    if (length(nalines) > 0) {
        msg <- sprintf(ngettext(length(nalines),
                                paste(" for", sQuote(name), "in %d line: %s"),
                                paste("s for", sQuote(name), "in %d lines: %s")),
                       length(nalines), paste(nalines+1, collapse = ", "))
        msg <- paste("empty field", msg, "\n", sep = "")
        out <- list(nalines = nalines, msg = msg)
    }
    else {
        out <- NULL
    }

    return(out)

}
