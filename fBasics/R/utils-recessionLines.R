
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
#  .recessionLines           Adds vertical recession lines to a plot
################################################################################


.recessionLines <-
function() 
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Adds vertical recession lines to a plot
    
    # Arguments:
    #   none
       
    # References:
    #   http://www.nber.org/cycles/cyclesmain.html
    
    # FUNCTION:
    
    # Add Recession Liones:
    abline(v = as.POSIXct(c(
        "1920-01-01", "1921-08-01",
        "1923-05-01", "1924-07-01",
        "1926-10-01", "1927-11-01",
        "1929-08-01", "1933-03-01",
        "1945-01-01", "1945-10-01",
        "1948-11-01", "1949-10-01",
        "1953-07-01", "1954-05-01",
        "1957-08-01", "1958-04-01",
        "1960-04-01", "1961-02-01",
        "1969-12-01", "1970-11-01",
        "1973-11-01", "1975-03-01",
        "1990-07-01", "1991-03-01",
        "2001-03-01", "2001-11-01", 
        "2007-12-07")), col = "grey")
}


################################################################################

