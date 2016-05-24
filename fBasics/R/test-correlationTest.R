
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
#  correlationTest       Performs correlation tests on two samples
#  pearsonTest           Pearson product moment correlation coefficient
#  kendallTest           Kendall's tau correlation test
#  spearmanTest          Spearman's rho correlation test
################################################################################


correlationTest <- 
function(x, y, method = c("pearson", "kendall", "spearman"), 
    title = NULL, description = NULL) 
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Performs Correlation Tests
        
    # FUNCTION:
    
    # Test:
    method = match.arg(method)
    if (method[1] == "pearson") {
        ans = pearsonTest(x, y, title = title, description = description) 
    }
    if (method[1] == "kendall") {
        ans = kendallTest(x, y, title = title, description = description) 
    }  
    if (method[1] == "spearman") {
       ans = spearmanTest(x, y, title = title, description = description)
    }
        
    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


pearsonTest <- 
function(x, y, title = NULL, description = NULL)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   A test for association between paired samples
    
    # Arguments:
    #   x - a numeric vector of data values.
    #   description - a brief description of the porject of type 
    #       character.
    #   title - a character string which allows for a project title.
    
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
    stopifnot(length(x) == length(y))
    
    # Test:
    two.sided = cor.test(x = x, y = y, alternative = "two.sided",  
        method = "pearson")
    less = cor.test(x = x, y = y, alternative = "less",  
        method = "pearson")
    greater = cor.test(x = x, y = y, alternative = "greater",  
        method = "pearson")  
        
    # Sample Estimates:   
    ESTIMATE = two.sided$estimate
    names(ESTIMATE) = "Correlation"
    test$estimate = ESTIMATE
    
    # Parameter
    DF = two.sided$parameter
    names(DF) = "Degrees of Freedom"
    test$parameter = DF
    
    # P Values:
    PVAL = c(
        two.sided$p.value, 
        less$p.value, 
        greater$p.value)
    names(PVAL) = c(
        "Alternative Two-Sided", 
        "Alternative      Less", 
        "Alternative   Greater")
    test$p.value = PVAL     
    
    # Confidences Levels:
    if (!is.null(two.sided$conf.int)) {
        CONF.INT = cbind(
            a = two.sided$conf.int, 
            b = less$conf.int, 
            c = greater$conf.int)
        dimnames(CONF.INT)[[2]] = c(
            "Two-Sided", 
            "     Less", 
            "  Greater")
        test$conf.int = CONF.INT  
    }  
    
    # Statistic:    
    STATISTIC = two.sided$statistic
    names(STATISTIC) = "t"
    test$statistic = STATISTIC
        
    # Add:
    if (is.null(title)) title = "Pearson's Correlation Test"
    if (is.null(description)) description = date()
    
    # Return Value:
    new("fHTEST",     
        call = call,
        data = list(x = x, y = y), 
        test = test,
        title = as.character(title), 
        description = as.character(description) ) 
}


# ------------------------------------------------------------------------------


kendallTest <- 
function(x, y, title = NULL, description = NULL)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   A test for association between paired samples
    
    # Arguments:
    #   x - a numeric vector of data values.
    #   description - a brief description of the porject of type 
    #       character.
    #   title - a character string which allows for a project title.
    
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
    stopifnot(length(x) == length(y))
    
    # Test:
    two.sided = cor.test(x = x, y = y, alternative = "two.sided",  
        method = "kendall")
    less = cor.test(x = x, y = y, alternative = "less",  
        method = "kendall")
    greater = cor.test(x = x, y = y, alternative = "greater",  
        method = "kendall") 
        
    # Exact Test:
    if (class(version) != "Sversion") {
        two.sided.exact = cor.test(x = x, y = y, exact = TRUE, 
            alternative = "two.sided",  method = "kendall")
        less.exact = cor.test(x = x, y = y, exact = TRUE, 
            alternative = "less", method = "kendall")
        greater.exact = cor.test(x = x, y = y, exact = TRUE, 
            alternative = "greater",  method = "kendall")
    } else {
        two.sided.exact = list()
        two.sided.exact$p.value = two.sided.exact$statistic = NA
        less.exact = list()
        less.exact$p.value = less.exact$statistic = NA
        greater.exact = list()
        greater.exact$p.value = greater.exact$statistic = NA 
    }
            
    # Sample Estimates:
    ESTIMATE = two.sided$estimate
    names(ESTIMATE) = "tau"
    test$estimate = ESTIMATE
    
    # P Values:
    PVAL = c(
        two.sided$p.value, 
        two.sided.exact$p.value,
        less$p.value, 
        less.exact$p.value, 
        greater$p.value, 
        greater.exact$p.value)
    if (is.na(two.sided.exact$p.value)) {
        names(PVAL) = c(
            "Alternative Two-Sided", 
            "Alternative Two-Sided | Exact",
            "Alternative      Less", 
            "Alternative      Less | Exact", 
            "Alternative   Greater", 
            "Alternative   Greater | Exact")
    } else {
        names(PVAL) = c(
            "Alternative         Two-Sided", 
            "Alternative Two-Sided | Exact",
            "Alternative              Less", 
            "Alternative      Less | Exact", 
            "Alternative           Greater", 
            "Alternative   Greater | Exact")
    }
    test$p.value = PVAL     
    
    # Statistic:
    # STATISTIC = c(
    #   two.sided$statistic, two.sided.exact$statistic,
    #   less$statistic, less.exact$statistic,
    #   greater$statistic, greater.exact$statistic)
    STATISTIC = c(
        two.sided$statistic, 
        two.sided.exact$statistic)
    # names(STATISTIC) = c(
    #   "z | Two-Sided", "T | Two-Sided | Exact",
    #   "z | Less", "T | Less | Exact",
    #   "z | Greater", "T | Greater | Exact")
    names(STATISTIC) = c(
        "z", 
        "T | Exact")
    test$statistic = STATISTIC
        
    # Add:
    if (is.null(title)) title = "Kendall's tau Correlation Test"
    if (is.null(description)) description = date()
    
    # Return Value:
    new("fHTEST",     
        call = call,
        data = list(x = x, y = y), 
        test = test,
        title = as.character(title), 
        description = as.character(description) ) 
}


# ------------------------------------------------------------------------------


spearmanTest <- 
function(x, y, title = NULL, description = NULL)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   A test for association between paired samples
    
    # Arguments:
    #   x - a numeric vector of data values.
    #   description - a brief description of the porject of type 
    #       character.
    #   title - a character string which allows for a project title.
    
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
    stopifnot(length(x) == length(y))
    
    # Test:
    two.sided = cor.test(x = x, y = y, alternative = "two.sided",  
        method = "spearman")
    less = cor.test(x = x, y = y, alternative = "less",  
        method = "spearman")
    greater = cor.test(x = x, y = y, alternative = "greater",  
        method = "spearman")
        
    # Sample Estimates:
    ESTIMATE = two.sided$estimate
    names(ESTIMATE) = "rho"
    test$estimate = ESTIMATE
    
    # P Values:
    PVAL = c(
        two.sided$p.value, 
        less$p.value, 
        greater$p.value)
    names(PVAL) = c(
        "Alternative Two-Sided", 
        "Alternative      Less", 
        "Alternative   Greater")
    test$p.value = PVAL   
    
    # Statistic:  
    STATISTIC = two.sided$statistic
    names(STATISTIC) = "S"
    test$statistic = STATISTIC
      
    # Add:
    if (is.null(title)) title = "Spearman's rho Correlation Test"
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

