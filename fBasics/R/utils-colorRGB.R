
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Library General Public License for more details.
#
# You should have received A copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################
# FUNCTION:                 DESCRIPTION:
#  .asRGB                    Converts any R color to RGB (red/green/blue)
#  .chcode                   Changes from one to another number system
#  .hex.to.dec               Converts heximal numbers do decimal numbers
#  .dec.to.hex               Converts decimal numbers do heximal numbers
################################################################################


.asRGB <-
function(col = rainbowPalette(64), alpha = FALSE)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Converts any R color to RGB (red/green/blue)

    # Arguments:
    #   col - vector of any of the three kind of R colors, i.e., either a
    #       color name (an element of colors()), a hexadecimal string of
    #       the form "#rrggbb", or an integer i meaning palette()[i].
    #   alpha - a logical value indicating whether alpha channel values
    #       should be returned.

    # FUNCTION:

    # Color Conversion:
    result <- col2rgb(col = col, alpha = alpha)

    # Return Value:
    t(result)
}


# ------------------------------------------------------------------------------


.chcode <-
function(b, base.in = 2, base.out = 10, digits = "0123456789ABCDEF")
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Changes from one to another number system

    # Arguments:
    #   b - number specified in the input base
    #   b.in - input base
    #   b.out - output base
    #   digits - digits string

    # Value:
    #   returns the input in the form represented by the output base

    # Author:
    #   Peter Wolf Universitaet Bielefeld
    #   from: http://tolstoy.newcastle.edu.au/R/help/05/04/2085.html

    # FUNCTION:

    # Change Number System:
    digits = substring(digits,1:nchar(digits),1:nchar(digits))
    if (length(base.in) == 1)
        base.in = rep(base.in, max(nchar(b) - 1))
    if (is.numeric(b))
        b = as.character(as.integer(b))
    b.num = lapply(strsplit(b, ""),
    function(x) {match(x, digits)-1} )
    result = lapply(b.num,
    function(x) {cumprod(rev(c(base.in,1))[1:length(x)]) %*% rev(x)} )
    number = unlist(result)
    # DW Print Output Suppressed
    # cat("decimal representation:",number,"\n")
    if (length(base.out) == 1) {
        base.out<-rep(base.out, 1+ceiling(log(max(number), base = base.out)))
    }
    n.base = length(base.out)
    result = NULL
    for(i in n.base:1){
        result = rbind(number %% base.out[i], result)
        number = floor(number/base.out[i])
    }
    result[]<-digits[result+1]
    ans = apply(result, 2, paste, collapse = "")

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


.hex.to.dec <-
function(b)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Converts heximal numbers do decimal numbers

    # Arguments:
    #   b - a heximal number

    # Value:
    #   returns a heximal numbers as decimal numbers

    # FUNCTION:

    # Hex to Bin:
    ans = as.numeric(.chcode(b, base.in = 16, base.out = 10))

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


.dec.to.hex <-
function(b)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Converts decimal numbers do heximal numbers

    # Arguments:
    #   x - a decimal number

    # Value:
    #   returns a decimal numbers as heximal numbers

    # FUNCTION:

    # Decimal to Hex:
    ans = .chcode(b, base.in = 10, base.out = 16)

    # Return Value:
    ans
}


################################################################################


