
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

# Copyrights (C) for this R-port:
#   1999 - 2012 Diethelm Wuertz, Zurich, <wuertz@itp.phys.ethz.ch>
#   2009 - 2012 Rmetrics Association, Zurich, www.rmetrics.org


################################################################################
# FUNCTION:               DESCRIPTION:
#  read.lynx               Uses the lynx browser to read from a web page
################################################################################


read.lynx <- 
    function(url, intern = TRUE, bin = NULL, pipe = FALSE, ...) 
{
    # A function implemented by Diethelm Wuertz and Yohan Chalabi

    # Description:
    #   Uses the lynx browser to read from a web page
    
    # Example:
    #   read.lynx("www.rmetrics.org", bin = "C:/cygwin/bin/lynx.exe")
    
    # FUNCTION:

    # Paste Arguments:
    ldots <- list(...)
    args <-
        if (length(ldots))
            paste(paste("-", names(ldots) , "=", sep = ""), ldots, sep = "", 
                collapse = " ")
        else
            ""

    # Compose Command:
    cmd <-
        if (is.null(bin))
            paste("lynx", args, "-dump -dont_wrap_pre", shQuote(url))
        else
            paste(bin, args, "-dump -dont_wrap_pre", shQuote(url))

    # Download URL:
    ans <-
        if (pipe)
            pipe(cmd)
        else
            system(cmd, intern = intern)

    # Return Value:
    ans
}


################################################################################
