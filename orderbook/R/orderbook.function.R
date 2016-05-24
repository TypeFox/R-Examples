################################################################################
##
##
## orderbook.function.R: Returns an object of class limitob
##
##
## orderbook is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 2 of the License, or
## (at your option) any later version.
##
## orderbook is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with orderbook.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

## Returns an orderbook object. For input it takes a data frame, and names for
## price, size, type, time, id, as well as what ASK and BID are denoted as.

orderbook <- function(file  = NULL)
{

        ## Create an empty current order book data frame

        current.ob <- data.frame(price = numeric(0), size =
                                 numeric(0), type = character(0), time
                                 = numeric(0), id = character(0))

        ## Check to see that the file is valid and can be opened

        obfile <- file(file, open = "r")

        stopifnot(isOpen(obfile, "r"))


        ## Read header file

         x <- scan(file, nlines = 1, sep = ",", what = "character",
              quiet = TRUE)

        ## Close file

        close(obfile)

        stopifnot(identical(x[1], "type"), identical(x[2], "time"),
                  identical(x[3], "id"), identical(x[4], "price"),
                  identical(x[5], "size"), identical(x[6], "type"))

        ## See if its a trader file or not

        if(identical(x[7], "status"))
            trader = TRUE
        else
            trader = FALSE


        ## Return a new order book object.

        invisible(new("orderbook",
                      current.ob   = current.ob,
                      file         = file,
                      trader       = trader
                      ))
}

