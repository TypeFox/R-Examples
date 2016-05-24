
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


################################################################################
# FUNCTION:             DESCRIPTION:
#  ks2Test               Performs a two sample Kolmogorov-Smirnov test
################################################################################


ks2Test <- 
function(x, y, title = NULL, description = NULL)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Two-sample Kolmogorov-Smirnov test.
    
    # Arguments:
    #   x - a numeric vector of data values.
    #   description - a brief description of the project of type 
    #       character.
    #   title - a character string which allows for a project title.
    
    # Note:
    #   A function partly copied from "stats"
    
    # FUNCTION:
    
    # Call:
    call = match.call()
    
    # Test:
    test = list()
    
    # Data Set Name:
    DNAME = paste(deparse(substitute(x)), "and", deparse(substitute(y)))
    test$data.name = DNAME
    
    # Convert Type:
    x = as.vector(x)
    y = as.vector(y)
    
    # Compute Test: 
    two.sided = ks.test(x = x, y = y, alternative = "two.sided")
    exact     = ks.test(x = x, y = y, exact = TRUE, alternative = "two.sided")
    less      = ks.test(x = x, y = y, alternative = "less")
    greater   = ks.test(x = x, y = y, alternative = "greater")

    # P Value:
    PVAL = c(
        two.sided$p.value, 
        exact$p.value, 
        less$p.value, 
        greater$p.value)
    names(PVAL) = c(
        "Alternative       Two-Sided", 
        "Alternative Exact Two-Sided",
        "Alternative            Less", 
        "Alternative         Greater")
    test$p.value = PVAL
    
    # Statistic:
    STATISTIC = c(
        two.sided$statistic, 
        less$statistic, 
        greater$statistic)
    names(STATISTIC) = c(
        "D | Two Sided", 
        "   D^- | Less", 
        "D^+ | Greater")
    test$statistic = STATISTIC
    
    # Add:
    if (is.null(title)) title = "Kolmogorov-Smirnov Two Sample Test"
    if (is.null(description)) description = date()
    
    # Return Value:
    new("fHTEST",     
        call = call,
        data = list(x = x, y = y), 
        test = test,
        title = as.character(title), 
        description = as.character(description) ) 
}


################################################################################

