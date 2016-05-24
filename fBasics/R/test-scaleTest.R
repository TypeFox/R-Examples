
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
#  scaleTest             Performs scale tests on two samples
#  .ansariTest           Ansari-Bradley test for differences in scale
#  .moodTest             Mood test for differences in scale
################################################################################


scaleTest <- 
function(x, y, method = c("ansari", "mood"), 
    title = NULL, description = NULL) 
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Performs Scale Tests
    
    # FUNCTION:
    
    # Test:
    method = match.arg(method)
    if (method == "ansari") {
        ans = .ansariTest(x, y, title = title, description = description) 
    }
    if (method == "mood") {
        ans = .moodTest(x, y, title = title, description = description) 
    }  
        
    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


.ansariTest <- 
function(x, y, title = NULL, description = NULL)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Ansari-Bradley's test for differences in scale
    
    # Arguments:
    #   x - a numeric vector of data values.
    #   description - a brief description of the porject of type 
    #       character.
    #   title - a character string which allows for a project title.
  
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
    
    # Test:
    two.sided = .ansari2Test(x = x, y = y, alternative = "two.sided",
        exact = FALSE, conf.int = TRUE, conf.level = 0.95)
    less = .ansari2Test(x = x, y = y, alternative = "less",
        exact = FALSE, conf.int = TRUE, conf.level = 0.95)
    greater = .ansari2Test(x = x, y = y, alternative = "greater",
        exact = FALSE, conf.int = TRUE, conf.level = 0.95)
        
    two.sided.exact = .ansari2Test(x = x, y = y, alternative = "two.sided",
        exact = TRUE, conf.int = TRUE, conf.level = 0.95)
    less.exact = .ansari2Test(x = x, y = y, alternative = "less",
        exact = TRUE, conf.int = TRUE, conf.level = 0.95)
    greater.exact = .ansari2Test(x = x, y = y, alternative = "greater",
        exact = TRUE, conf.int = TRUE, conf.level = 0.95)

    # Statistic:
    STATISTIC = c(two.sided$statistic)
    names(STATISTIC) = "AB"
    test$statistic = STATISTIC
    
    # P Values:
    PVAL = c(
        two.sided$p.value, 
        two.sided.exact$p.value,
        less$p.value, 
        less.exact$p.value,
        greater$p.value, 
        greater.exact$p.value)
    names(PVAL) = c(
        "Alternative Two-Sided        ", 
        "Alternative Two-Sided | Exact",
        "Alternative      Less        ", 
        "Alternative      Less | Exact",
        "Alternative   Greater        ", 
        "Alternative   Greater | Exact")
    test$p.values = PVAL
    
    # Confidence Levels:
    CONF.INT = cbind(
        a = two.sided$conf.int, 
        b = two.sided.exact$conf.int,
        c = less$conf.int, 
        d = less.exact$conf.int,
        e = greater$conf.int, 
        f = greater.exact$conf.int)
    # For Splus compatibility use named a CONF.INT
    # and dimnames instead of colnames!
    dimnames(CONF.INT)[[2]] = c(
        "Two-Sided | Asymptotic ", 
        "Two-Sided |      Exact ", 
        "Less      | Asymptotic ", 
        "Less      |      Exact ",
        "Greater   | Asymptotic ", 
        "Greater   |      Exact ")
    test$conf.int = CONF.INT
    
    # Add:
    if(is.null(title)) title = "Ansari-Bradley Test for Scale"
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


.ansari2Test <- 
function(x, y, alternative = c("two.sided", "less", "greater"),
    exact = TRUE, conf.int = FALSE, conf.level = 0.95, ...)
{
    # A function implemented by Diethelm Wuertz
    
    # Arguments:
    #   x - numeric vector of data values. 
    #   y - numeric vector of data values. 
    #   alternative - indicates the alternative hypothesis and must 
    #       be one of "two.sided", "greater" or "less". You can specify 
    #       just the initial letter. 
    #   exact - a logical indicating whether an exact p-value should 
    #       be computed. 
    #   conf.int - a logical,indicating whether a confidence interval 
    #       should be computed. 
    #   conf.level - confidence level of the interval. 

    # FUNCTION:
    
    # Return Value:
    ansari.test(x, y, alternative, exact, conf.int, conf.level, ...)     
}


# ------------------------------------------------------------------------------


.moodTest <- 
function(x, y, title = NULL, description = NULL)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Performs Mood's two-sample test for a difference in 
    #   scale parameters. 
    
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

    # Check Data:
    m = length(x)
    n = length(y)
    s = m + n
    if (s < 3) stop("not enough observations")
  
    # Statistic:
    r = rank(c(x, y)) 
    # From R:
    # z = ((sum((r[1:length(x)] - (s + 1) / 2)^2) - m * (s^2 - 1) / 12)
    #      / sqrt(m * n * (s + 1) * (s + 2) * (s - 2) / 180) )
    # To run also under S-Plus use ...      
    a = sum( (r[1:length(x)]-0.5*(s+1))^2 ) - m*(s*s-1)/12
    b = sqrt(as.numeric(m) * n * (s + 1) * (s + 2) * (s - 2) / 180) 
    STATISTIC = a/b
    names(STATISTIC) = "Z"
    test$statistic = STATISTIC
    
    # P Values:
    p = pnorm(STATISTIC)
    PVAL = c(2 * min(p, 1 - p), p, 1 - p)
    names(PVAL) = c(
        "Alternative Two-Sided", 
        "Alternative      Less",
        "Alternative   Greater")
    test$p.value = PVAL

    # Add:
    if(is.null(title)) title = "Mood Two-Sample Test of Scale"
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

