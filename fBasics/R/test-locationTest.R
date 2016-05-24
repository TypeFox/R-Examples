
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
# FUNCTION:             LOCATION TESTS:
#  locationTest          Performs locations tests on two samples
#  .tTest                Unpaired t test for differences in mean
#  .kw2Test              Kruskal-Wallis test for differences in locations
################################################################################


locationTest <-
function(x, y, method = c("t", "kw2"),
    title = NULL, description = NULL)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Correlation Tests

    # FUNCTION:

    # Test:
    method = match.arg(method)
    if (method == "t") {
        ans = .tTest(x, y, title = title, description = description)
    }
    if (method == "kw2") {
        ans = .kw2Test(x, y, title = title, description = description)
    }

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


.tTest <-
function(x, y, title = NULL, description = NULL)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #    Tests if two population means are equal.

    # Arguments:
    #   x, y - two numeric vector of data values or time series objects
    #   description - a brief description of the porject of type character.
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

    # Asymptotic Test:
    two.sided = t.test(x = x, y = y, alternative = "two.sided",
        mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95)
    less = t.test(x = x, y = y, alternative = "less",
        mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95)
    greater = t.test(x = x, y = y, alternative = "greater",
        mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95)

    # Assume Equal Variances:
    two.sided.equal = t.test(x = x, y = y, alternative = "two.sided",
        mu = 0, paired = FALSE, var.equal = TRUE, conf.level = 0.95)
    less.equal = t.test(x = x, y = y, alternative = "less",
        mu = 0, paired = FALSE, var.equal = TRUE, conf.level = 0.95)
    greater.equal = t.test(x = x, y = y, alternative = "greater",
        mu = 0, paired = FALSE, var.equal = TRUE, conf.level = 0.95)

    # Sample Estimates:
    PARAMETER = c(length(x), length(y), 0)
    names(PARAMETER) = c(
        "x Observations",
        "y Observations",
        "mu")
    test$parameter = PARAMETER

    # Sample Estimates:
    ESTIMATE = c(two.sided$estimate, var(x), var(y))
    names(ESTIMATE) = c("Mean of x", "Mean of y", "Var  of x", "Var  of y")
    test$estimate = ESTIMATE

    # P Values:
    PVAL = c(
        two.sided$p.value,
        less$p.value,
        greater$p.value,
        two.sided.equal$p.value,
        less.equal$p.value,
        greater.equal$p.value)
    names(PVAL) = c(
        "Alternative Two-Sided",
        "Alternative      Less",
        "Alternative   Greater",
        "Alternative Two-Sided | Equal Var",
        "Alternative      Less | Equal Var",
        "Alternative   Greater | Equal Var")
    test$p.value = PVAL

    # Statistic:
    STATISTIC = c(
        two.sided$statistic,
        two.sided.equal$statistic)
    names(STATISTIC) = c(
        "            T",
        "T | Equal Var")
    test$statistic = STATISTIC

    # Confidence Intervals:
    CONF.INT = cbind(
        a = two.sided$conf.int,
        b = less$conf.int,
        c = greater$conf.int,
        d = two.sided.equal$conf.int,
        e = less.equal$conf.int,
        f = greater.equal$conf.int)
    # For Splus compatibility use named a CONF.INT
    # and dimnames instead of colnames!
    dimnames(CONF.INT)[[2]] = c(
        "Two-Sided",
        "     Less",
        "  Greater",
        "Two-Sided | Equal Var",
        "     Less | Equal Var",
        "  Greater | Equal Var")
    test$conf.int = CONF.INT

    # Add:
    if (is.null(title)) title = "t Test"
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


.kw2Test <-
function(x, y, title = NULL, description = NULL)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Performs a Kruskal-Wallis rank sum test of the null that
    #   the location parameters of the distribution of x are the
    #   same in each group (sample). The alternative is that they
    #   differ in at least one.

    # Arguments:
    #   x, y - two numeric vector of data values or time series objects
    #   description - a brief description of the porject of type character.
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

    # Sample Estimates:
    ESTIMATE = c(mean(x), mean(y), var(x), var(y))
    names(ESTIMATE) = c("Mean of x", "Mean of y", "Var  of x", "Var  of y")
    test$estimate = ESTIMATE

    # Parameter:
    PARAMETER = c(length(x), length(y))
    names(PARAMETER) = c(
        "x Observations",
        "y Observations")
    test$parameter = PARAMETER

    # Operate on Lists:
    x = list(x = x, y = y)
    if (length(x) < 2) stop("x must be a list with at least 2 elements")
    k = length(x)
    l = sapply(x, "length")
    g = factor(rep(1 : k, l))
    x = unlist(x)

    # Test:
    n = length(x)
    if (n < 2) stop("not enough observations")
    r = rank(x)
    TIES = table(x)

    # Statistic:
    STATISTIC = sum(tapply(r, g, "sum")^2 / tapply(r, g, "length"))
    STATISTIC = ((12 * STATISTIC / (n * (n + 1)) - 3 * (n + 1)) /
        (1 - sum(TIES^3 - TIES) / (n^3 - n)))
    names(STATISTIC) = "KW chi-squared"
    test$statistic = STATISTIC

    # P Value:
    PVAL = 1 - pchisq(STATISTIC, 1)
    names(PVAL) = ""
    test$p.value = PVAL

    # Add:
    if(is.null(title)) title = "Kruskal-Wallis Two Sample Test"
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

