
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
#  varianceTest          Performs variance tests on two samples
#  .varfTest             F test for differences in variances
#  .bartlett2Test        Bartlett's test for differences in variances
#  .fligner2Test         Fligner-Killeen test for differences in variances
################################################################################


varianceTest <- 
function(x, y, method = c("varf", "bartlett", "fligner"), 
    title = NULL, description = NULL) 
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Performs Variance Tests
       
    # FUNCTION:
    
    # Test:
    method = match.arg(method)
    if (method == "varf") {
        ans = .varfTest(x, y, title = title, description = description) 
    }
    if (method == "bartlett") {
        ans = .bartlett2Test(x, y, title = title, description = description) 
    }  
    if (method == "fligner") {
        ans = .fligner2Test(x, y, title = title, description = description) 
    } 
        
    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


.varfTest <- 
function(x, y, title = NULL, description = NULL)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Performs an F test to compare the variances of two samples 
    #   from normal populations. 

    # Arguments:
    #   x, y - a numeric vector of data values.
    #   description - a brief description of the porject of type 
    #       character.
    #   title - a character string which allows for a project title.
    
    # Notes:
    #   A modified copy originally from R's ctest package Version 1.8.1

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
    
    # Estimate - Hypothesized Equal Variance:
    ratio = 1
    DF.x = length(x) - 1
    DF.y = length(y) - 1
    VAR.x = var(x)
    VAR.y = var(y)
    ESTIMATE = VAR.x / VAR.y / ratio
    names(ESTIMATE) = "Ratio of Variances"
    test$estimate = ESTIMATE
   
    # Parameter:
    PARAMETER = c(ratio, DF.x, DF.y)
    names(PARAMETER) = c(
        "Hypothesized Ratio", 
        "Numerator   df", 
        "Denumerator df")
    test$parameter = PARAMETER
    
    # Statistic:
    STATISTIC = ESTIMATE / ratio
    names(STATISTIC) = "F"
    test$statistic = STATISTIC
    
    # P Value:
    p = pf(STATISTIC, DF.x, DF.y)
    PVAL = c(
        two.sided = 2 * min(p, 1 - p), 
        less = p, 
        greater = 1 - p)
    names(PVAL) = c(
        "Alternative Two-Sided", 
        "Alternative      Less",
        "Alternative   Greater")
    test$p.value = PVAL
    
    # Confidence Interval:
    conf.level = 0.95
    B = (1 - conf.level) / 2
    two.sided = c(ESTIMATE/qf(1-B, DF.x, DF.y), ESTIMATE/qf(B, DF.x, DF.y))
    less = c(0, ESTIMATE/qf(1-conf.level, DF.x, DF.y)) 
    greater = c(ESTIMATE/qf(conf.level, DF.x, DF.y), Inf) 
    CONF.INT = cbind(
        a = two.sided, 
        b = less, 
        c = greater)
    dimnames(CONF.INT)[[2]] = c(
        "Two-Sided", 
        "     Less", 
        "  Greater")
    test$conf.int = CONF.INT
   
    # Add:
    if(is.null(title)) title = "F Test of Variances"
    if(is.null(description)) description = date()
        
    # Return Value:
    new("fHTEST",     
        call = call,
        data = list(x = x, y = y), 
        test = test,
        title = as.character(title), 
        description = as.character(description) )
}


# ------------------------------------------------------------------------------


.bartlett2Test <- 
function(x, y, title = NULL, description = NULL)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Bartlett's test for differences in variances
    
    # Arguments:
    #   x - a numeric vector of data values.
    
    # Note:
    #   # A function linked to "stats"

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
    
    # Settings:
    x = list(x = x, y = y)
    k = length(x)
    n = sapply(x, "length") - 1
    v = sapply(x, "var")
    n.total = sum(n)
    v.total = sum(n * v) / n.total
    
    # Statistic:
    STATISTIC = ((n.total * log(v.total) - sum(n * log(v))) /
        (1 + (sum(1 / n) - 1 / n.total) / (3 * (k - 1))))
    names(STATISTIC) = "Bartlett's Chi-squared"
    test$statistic = STATISTIC
    
    # P Value:
    PVAL = 1 - pchisq(STATISTIC, 1)    
    names(PVAL) = ""
    test$p.value = PVAL  
    
    # Add:
    if(is.null(title)) title = "Bartlett Test for Homogeneity of Variances"
    if(is.null(description)) description = date()  
    
    # Return Value:
    new("fHTEST",     
        call = call,
        data = list(x = x, y = y), 
        test = test,
        title = as.character(title), 
        description = as.character(description) )
}


# ------------------------------------------------------------------------------
  

.fligner2Test <- 
function(x, y, title = NULL, description = NULL)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Fligner-Killeen's rank based test for homogeneity of variances
 
    # Arguments:
    #   x - a numeric vector of data values.
    
    # Note:
    #   A function linked to "stats"

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
    
    # Settings:
    x = list(x = x, y = y)
    k = length(x)
    l = sapply(x, "length")
    g = factor(rep(1 : k, l))
    x = unlist(x)
   
    # Statistic:
    n = length(x)
    x = unlist(tapply(x, g, function(u) u - median(u)))
    a = qnorm((1 + rank(abs(x)) / (n + 1)) / 2)
    STATISTIC = sum(tapply(a, g, "sum")^2 / tapply(a, g, "length"))
    STATISTIC = (STATISTIC - n * mean(a)^2) / var(a)
    names(STATISTIC) = "FK:med chi-squared"
    test$statistic = STATISTIC

    # P Value:
    PVAL = 1 - pchisq(STATISTIC, 1)
    names(PVAL) = ""
    test$p.value = PVAL
    
    # Add:
    if(is.null(title)) title = "Fligner-Killeen Test for Homogeneity of Variances"
    if(is.null(description)) description = date()  
    
    # Return Value:
    new("fHTEST",     
        call = call,
        data = list(x = x, y = y), 
        test = test,
        title = as.character(title), 
        description = as.character(description) )
}


################################################################################

