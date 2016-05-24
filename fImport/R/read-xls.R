
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
# You should have received a copy of the GNU Library General 
# Public License along with this library; if not, write to the 
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, 
# MA  02111-1307  USA

# Copyrights (C) for this R-port:
#   1999 - 2012 Diethelm Wuertz, Zurich, <wuertz@itp.phys.ethz.ch>
#   2009 - 2012 Rmetrics Association, Zurich, www.rmetrics.org

################################################################################
# FUNCTION:             DESCRIPTION:
#  read.xls              Downloads an xls file and converts it in a csv file
################################################################################


read.xls <- 
function(url, sheet=1, lines=-1, verbose=FALSE, encoding="unknown")
{
    # Description:
    #   Downloads an xls and converts it in a csv file
    
    # Arguments:
    #   url - a character string specifying the URL
    #   sheet - an integer denoting whih sheet should be extracted
    #   lines - a negative integer with the lines to be skipped
    #   verbose - a logical decides about verbose mode
    #   encoding - a character string with the type of encodeing
    
    # Example:
    #   read.xls("http://www.jrsainfo.org/jabg/state_data2/Tribal_Data00.xls")
    
    # FUNCTION:
    
    # Download:
    con <- .xls2csv(url, sheet, verbose)
    ans = readLines(con, n = lines, encoding = encoding)
    ans = gsub('"', "", ans)
    file.remove(summary(con)$description)
    
    # Verbose ?
    if(verbose) print(ans)

    # Return Value:
    invisible(ans)
}


################################################################################
 
