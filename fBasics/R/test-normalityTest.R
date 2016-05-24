
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
#  ksnormTest            One sample Kolmogorov-Smirnov normality test
#  shapiroTest           Shapiro-Wilk normality test
#  jarqueberaTest        Jarque-Bera normality test
#  dagoTest              D'Agostino normality test
#   .skewness.test        ... internal function called by dagoTest
#   .kurtosis.test        ... internal function called by dagoTest
#   .omnibus.test         ... internal function called by dagoTest
# FUNCTION:             DESCRIPTION:
#  normalTest            Normality tests S-Plus compatible
# FUNCTION:             FROM NORTEST PACKAGE:
#  adTest                Anderson-Darling normality test
#  cvmTest               Cramer-von Mises normality test
#  lillieTest            Lilliefors (Kolmogorov-Smirnov) normality test
#  pchiTest              Pearson chi-square normality test
#  sfTest                Shapiro-Francia normality test
# FUNCTION ADDON:       AUGMENTED FINITE SAMPLE JB TEST:
#  jbTest                Performs finite sample adjusted JB LM and ALM test
#  .jb.test               S3 version type finite sample adjusted JB test
################################################################################


ksnormTest <-
function(x, title = NULL, description = NULL)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Performs a one-sample Kolmogorov-Smirnov normality test

    # Arguments:
    #   x - a numeric vector or an univariate 'timeSeries' object.
    #   description - a brief description of the porject of type
    #       character.
    #   title - a character string which allows for a project title.

    # FUNCTION:

    # Convert Type:
    if (class(x) == "fREG") x = residuals(x)
    x = as.vector(x)

    # Call:
    call = match.call()

    # Test:
    test = ks.test(x, "pnorm", alternative = "two.sided")
    less = ks.test(x, "pnorm", alternative = "less")
    greater = ks.test(x, "pnorm", alternative = "greater")
    PVAL = c(test$p.value, less$p.value, greater$p.value)
    names(PVAL) = c(
        "Alternative Two-Sided",
        "Alternative      Less",
        "Alternative   Greater")
    test$metod = "Kolmogorov - Smirnow One Sample Normality Test"
    test$p.value = PVAL
    class(test) = "list"

    # Add:
    if (is.null(title)) title = test$method
    if (is.null(description)) description = description()

    # Return Value:
    new("fHTEST",
        call = call,
        data = list(x = x),
        test = test,
        title = as.character(title),
        description = as.character(description) )
}

# ------------------------------------------------------------------------------


shapiroTest <-
function(x, title = NULL, description = NULL)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Performs the Shapiro-Wilk test for normality.

    # Arguments:
    #   x - a numeric vector or an univariate 'timeSeries' object.
    #   description - a brief description of the porject of type
    #       character.
    #   title - a character string which allows for a project title.

    # Note:
    #   A function linked to "stats"

    # FUNCTION:

    # Convert Type:
    if (class(x) == "fREG") x = residuals(x)
    x = as.vector(x)

    # Call:
    call = match.call()

    # Works both in R and SPlus:
    test = shapiro.test(x)
    names(test$p.value) = ""
    class(test) = "list"

    # Add:
    if (is.null(title)) title = "Shapiro - Wilk Normality Test"
    if (is.null(description)) description = description()

    # Return Value:
    new("fHTEST",
        call = call,
        data = list(x = x),
        test = test,
        title = as.character(title),
        description = as.character(description) )
}


# ------------------------------------------------------------------------------


jarqueberaTest <-
function(x, title = NULL, description = NULL)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #  Jarque-Bera Test for Normality

    # Arguments:
    #   x - a numeric vector or an univariate 'timeSeries' object.
    #   description - a brief description of the porject of type
    #       character.
    #   title - a character string which allows for a project title.

    # Authors:
    #   from A. Trapletti's tseries Package.

    # FUNCTION:

    # Data Set Name:
    DNAME = deparse(substitute(x))

    # Convert Type:
    if (class(x) == "fREG") x = residuals(x)
    x = as.vector(x)

    # Call:
    call = match.call()

    # Test:

    n = length(x)
    m1 = sum(x)/n
    m2 = sum((x-m1)^2)/n
    m3 = sum((x-m1)^3)/n
    m4 = sum((x-m1)^4)/n
    b1 = (m3/m2^(3/2))^2
    b2 = (m4/m2^2)

    # Statistic:
    STATISTIC = n*b1/6+n*(b2-3)^2/24
    names(STATISTIC) = "X-squared"

    # P Value:
    PVAL = 1 - pchisq(STATISTIC,df = 2)
    names(PVAL) = "Asymptotic p Value"

    # Method:
    METHOD = "Jarque Bera Test"

    # Result:
    test = list(
        statistic = STATISTIC,
        p.value = PVAL,
        method = METHOD,
        data.name = DNAME)

    # Add:
    if (is.null(title)) title = "Jarque - Bera Normalality Test"
    if (is.null(description)) description = description()

    # Return Value:
    new("fHTEST",
        call = call,
        data = list(x = x),
        test = test,
        title = as.character(title),
        description = as.character(description) )
}


# ------------------------------------------------------------------------------


.skewness.test <-
function(x)
{
    # Internal Function for D'Agostino Normality Test:

    # Note:
    #   D'Agostino Test
    #   http://adela.karlin.mff.cuni.cz/~klaster/vyuka/
    #   Materià¥‡ly pro cvicenö€Œ, kterà¥‡ byla v labu, jsou zde: cv01.txt,
    #   cv02.txt, cv03.txt, cv05.txt, cv06.txt, data maths, police a
    #   vysky a makro dagost.r. Vö€œber nejakö€œch prö€Œkladu ze cvicenö€Œ je
    #   tady.
    #   Program R lze zdarma (GNU General Public Licence) stà¥‡hnout z
    #   www.r-project.org. Alespon k letmà¤¼mu nahlà¤¼dnutö€Œ doporucuji tà¤¼ž
    #   minimanuà¥‡l An Introduction to R, kterö€œ roste tamtà¤¼ž. Dalšö€Œ
    #   materià¥‡ly vcetne dvou zacà¥‡tecnickö€œch prö€Œrucek najdete na
    #   strà¥‡nkà¥‡ch Dr. Kulicha.

    # FUNCTION:

    DNAME = deparse(substitute(x))
    if (exists("complete.cases")) {
        test = complete.cases(x)
    } else {
        test = !is.na(x)
    }
    x = x[test]
    n = length(x)
    if (n < 8) stop("Sample size must be at least 8")
    meanX = mean(x)
    s =  sqrt(mean((x-meanX)**2))
    a3 = mean((x-meanX)**3)/s**3
    SD3 = sqrt(6*(n-2)/((n+1)*(n+3)))
    U3 = a3/SD3
    b  = (3*(n**2+27*n-70)*(n+1)*(n+3))/((n-2)*(n+5)*(n+7)*(n+9))
    W2 = sqrt(2*(b-1))-1
    delta = 1/sqrt(log(sqrt(W2)))
    a = sqrt(2/(W2-1))
    Z3 = delta*log((U3/a)+sqrt((U3/a)**2+1))
    pZ3 = 2*(1-pnorm(abs(Z3),0,1))
    names(Z3) = "Z3"

    # Result:
    RVAL = list(
        statistic = Z3,
        p.value = pZ3,
        method = "D'Agostino Skewness Normality Test",
        data.name = DNAME)

    # Return Value:
    class(RVAL) = "htest"
    RVAL
}


# ------------------------------------------------------------------------------


.kurtosis.test <-
function(x)
{
    # Internal Function for D'Agostino Normality Test:

    # FUNCTION:

    DNAME = deparse(substitute(x))

    if (exists("complete.cases")) {
        test = complete.cases(x)
    } else {
        test = !is.na(x)
    }
    x = x[test]
    n = length(x)
    if (n < 20) stop("Sample size must be at least 20")
    meanX = mean(x)
    s =  sqrt(mean((x-meanX)**2))
    a4 = mean((x-meanX)**4)/s**4
    SD4 = sqrt(24*(n-2)*(n-3)*n/((n+1)**2*(n+3)*(n+5)))
    U4 = (a4-3+6/(n+1))/SD4
    B = (6*(n*n-5*n+2)/((n+7)*(n+9)))*sqrt((6*(n+3)*(n+5))/(n*(n-2)*(n-3)))
    A = 6+(8/B)*((2/B)+sqrt(1+4/(B**2)))
    jm = sqrt(2/(9*A))
    pos = ((1-2/A)/(1+U4*sqrt(2/(A-4))))**(1/3)
    Z4 = (1-2/(9*A)-pos)/jm
    pZ4 = 2*(1-pnorm(abs(Z4),0,1))
    names(Z4) = "Z4"

    # Result:
    RVAL = list(
        statistic = Z4,
        p.value = pZ4,
        method = "D'Agostino Kurtosis Normality Test",
        data.name = DNAME)

    # Return Value:
    class(RVAL) = "htest"
    RVAL
}


# ------------------------------------------------------------------------------


.omnibus.test <-
function(x)
{
    # Internal Function for D'Agostino Normality Test:

    # FUNCTION:

    DNAME = deparse(substitute(x))
    if (exists("complete.cases")) {
        test = complete.cases(x)
    } else {
        test = !is.na(x)
    }
    x = x[test]
    n = length(x)
    if (n < 20) stop("sample size must be at least 20")
    meanX = mean(x)
    s =  sqrt(mean((x-meanX)**2))
    a3 = mean((x-meanX)**3)/s**3
    a4 = mean((x-meanX)**4)/s**4
    SD3 = sqrt(6*(n-2)/((n+1)*(n+3)))
    SD4 = sqrt(24*(n-2)*(n-3)*n/((n+1)**2*(n+3)*(n+5)))
    U3 = a3/SD3
    U4 = (a4-3+6/(n+1))/SD4
    b  = (3*(n**2+27*n-70)*(n+1)*(n+3))/((n-2)*(n+5)*(n+7)*(n+9))
    W2 = sqrt(2*(b-1))-1
    delta = 1/sqrt(log(sqrt(W2)))
    a = sqrt(2/(W2-1))
    Z3 = delta*log((U3/a)+sqrt((U3/a)**2+1))
    B = (6*(n*n-5*n+2)/((n+7)*(n+9)))*sqrt((6*(n+3)*(n+5))/(n*(n-2)*(n-3)))
    A = 6+(8/B)*((2/B)+sqrt(1+4/(B**2)))
    jm = sqrt(2/(9*A))
    pos = ((1-2/A)/(1+U4*sqrt(2/(A-4))))**(1/3)
    Z4 = (1-2/(9*A)-pos)/jm
    omni = Z3**2+Z4**2
    pomni = 1-pchisq(omni,2)
    names(omni) = "Chi2"

    # Result:
    RVAL = list(
        statistic = omni,
        method = "D'Agostino Omnibus Normality Test",
        p.value = pomni,
        data.name = DNAME)

    # Return Value:
    class(RVAL) = "htest"
    RVAL
}


# ------------------------------------------------------------------------------


dagoTest =
function(x, title = NULL, description = NULL)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Performs the D'Agostino normality test

    # Arguments:
    #   x - a numeric vector or an univariate 'timeSeries' object.
    #   description - a brief description of the porject of type
    #       character.
    #   title - a character string which allows for a project title.

    # Source:
    #   This function was inspired by ...
    #   http://adela.karlin.mff.cuni.cz/~klaster/vyuka/

    # FUNCTION:

    # Data Set Name:
    DNAME = deparse(substitute(x))

    # Convert Type:
    if (class(x) == "fREG") x = residuals(x)
    x = as.vector(x)

    # Call:
    call = match.call()

    # Test:
    ans = NA
    test = .omnibus.test(x)
    skew = .skewness.test(x)
    kurt = .kurtosis.test(x)
    test$data.name = DNAME
    PVAL = c(test$p.value, skew$p.value, kurt$p.value)
    names(PVAL) = c(
        "Omnibus  Test",
        "Skewness Test",
        "Kurtosis Test")
    test$p.value = PVAL
    STATISTIC = c(test$statistic, skew$statistic, kurt$statistic)
    names(STATISTIC) = c(
        "Chi2 | Omnibus",
        "Z3  | Skewness",
        "Z4  | Kurtosis")
    test$statistic = STATISTIC
    class(test) = "list"

    # Add:
    if (is.null(title)) title = "D'Agostino Normality Test"
    if (is.null(description)) description = description()

    # Return Value:
    new("fHTEST",
        call = call,
        data = list(x = x),
        test = test,
        title = as.character(title),
        description = as.character(description) )
}


################################################################################


normalTest <-
function(x, method = c("sw", "jb"), na.rm = FALSE)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Shapiro-Wilk and Jarque-Bera Test

    # Notes:
    #   This function is for S-Plus compatibility

    # FUNCTION:

    # Convert Type:
    if (class(x) == "fREG") x = residuals(x)
    x = as.vector(x)
    if (na.rm) x = x[!is.na(x)]

    # Method:
    #   Don't use: method = match.arg(method)
    method = method[1]

    # Test:
    if (method == "sw") {
        ans = shapiroTest(x)
    } else if (method == "jb") {
       ans = jarqueberaTest(x)
    }

    # Additional Tests:
    if (method == "ks") {
        ans = ksnormTest(x)
    }
    if (method == "da") {
        ans = dagoTest(x)
    }
    if (method == "ad") {
        ans = adTest(x)
    }

    # Return Value:
    ans
}


################################################################################
# FROM NORTEST PACKAGE:


adTest <-
function(x, title = NULL, description = NULL)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Performs the Anderson-Darling normality test

    # Arguments:
    #   x - a numeric vector of data values.
    #   description - a brief description of the porject of type
    #       character.
    #   title - a character string which allows for a project title.

    # Source:
    #   Package: nortest
    #   Title: Tests for Normality
    #   Version: 1.0
    #   Author: Juergen Gross
    #   Description: 5 omnibus tests for the composite hypothesis of normality
    #   Maintainer: Juergen Gross <gross@statistik.uni-dortmund.de>
    #   License: GPL version 2 or newer

    # Thanks:
    #   to Spencer Grave for contributions.

    # FUNCTION:

    # Data Set Name:
    DNAME = deparse(substitute(x))

    # Convert Type:
    if (class(x) == "fREG") x = residuals(x)
    x = as.vector(x)

    # Call:
    call = match.call()

    # Test:
    x = sort(x)
    n = length(x)
    if (n < 8) stop("sample size must be greater than 7")
    var.x <- var(x)

    if(var.x > 0){
        p = pnorm((x - mean(x))/sqrt(var.x))
        h = (2 * seq(1:n) - 1) * (log(p) + log(1 - rev(p)))
        ### DW modified:
        h = h[is.finite(h)]
        n = length(h)
        ### Continue:
        A = -n - mean(h)
        AA = (1 + 0.75/n + 2.25/n^2) * A

        if (AA < 0.2) {
            PVAL = 1 - exp(-13.436 + 101.14 * AA - 223.73 * AA^2)
        } else if (AA < 0.34) {
            PVAL = 1 - exp(-8.318 + 42.796 * AA - 59.938 * AA^2)
        } else if (AA < 0.6) {
            PVAL = exp(0.9177 - 4.279 * AA - 1.38 * AA^2)
        } else {
            PVAL = exp(1.2937 - 5.709 * AA + 0.0186 * AA^2)
        }
        # Result:
        if (PVAL > 1) {
            PVAL = 1 # was NA, suggested by Spencer Graves to modify
            W = NA
        }
    } else {
        A <- Inf
        PVAL <- 0
    }
    names(PVAL) = ""

    test = list(
        statistic = c(A = A),
        p.value = PVAL,
        method = "Anderson - Darling Normality Test",
        data.name = DNAME)

    # Add:
    if (is.null(title)) title = "Anderson - Darling Normality Test"
    if (is.null(description)) description = description()

    # Return Value:
    ans = new("fHTEST",
        call = call,
        data = list(x = x),
        test = test,
        title = as.character(title),
        description = as.character(description))
}


# ------------------------------------------------------------------------------


cvmTest <-
function(x, title = NULL, description = NULL)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Performs the Cramer-von Mises normality test

    # Arguments:
    #   x - a numeric vector of data values.
    #   description - a brief description of the porject of type
    #       character.
    #   title - a character string which allows for a project title.

    # Notes:
    #   A copy from contributed R-package 'nortest'
    #   Source:
    #       Package: nortest
    #       Title: Tests for Normality
    #       Version: 1.0
    #       Author: Juergen Gross
    #       Description: 5 omnibus tests for the composite hypothesis of normality
    #       Maintainer: Juergen Gross <gross@statistik.uni-dortmund.de>
    #       License: GPL version 2 or newer

    # FUNCTION:

    # Data Set Name:
    DNAME = deparse(substitute(x))

    # Convert Type:
    if (class(x) == "fREG") x = residuals(x)
    x = as.vector(x)

    # Call:
    call = match.call()

    # Test:
    x = sort(x)
    n = length(x)
    if (n < 8) stop("sample size must be greater than 7")
    p = pnorm((x - mean(x))/sqrt(var(x)))
    W = (1/(12 * n) + sum((p - (2 * seq(1:n) - 1)/(2 * n))^2))
    WW = (1 + 0.5/n) * W

    # P Value:
    if (WW < 0.0275) {
        PVAL = 1 - exp(-13.953 + 775.5 * WW - 12542.61 * WW^2)
    } else if (WW < 0.051) {
        PVAL = 1 - exp(-5.903 + 179.546 * WW - 1515.29 * WW^2)
    } else if (WW < 0.092) {
        PVAL = exp(0.886 - 31.62 * WW + 10.897 * WW^2)
    } else {
        PVAL = exp(1.111 - 34.242 * WW + 12.832 * WW^2)
    }

    # Result:
    if (PVAL > 1) {
        PVAL = 1
        W = NA
    }
    names(PVAL) = ""

    test = list(
        statistic = c(W = W),
        p.value = PVAL,
        method = "Cramer - von Mises Test",
        data.name = DNAME)

    # Add:
    if (is.null(title)) title = "Cramer - von Mises Normality Test"
    if (is.null(description)) description = description()

    # Return Value:
    new("fHTEST",
        call = call,
        data = list(x = x),
        test = test,
        title = as.character(title),
        description = as.character(description) )
}


# ------------------------------------------------------------------------------


lillieTest <-
function(x, title = NULL, description = NULL)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Performs the Lilliefors (Kolmogorov-Smirnov) normality test

    # Arguments:
    #   x - a numeric vector of data values.
    #   description - a brief description of the porject of type
    #       character.
    #   title - a character string which allows for a project title.

    # Notes:
    #   A copy from contributed R-package 'nortest'
    #   Source:
    #       Package: nortest
    #       Title: Tests for Normality
    #       Version: 1.0
    #       Author: Juergen Gross
    #       Description: 5 omnibus tests for the composite hypothesis of normality
    #       Maintainer: Juergen Gross <gross@statistik.uni-dortmund.de>
    #       License: GPL version 2 or newer

    # FUNCTION:

    # Data Set Name:
    DNAME = deparse(substitute(x))

    # Convert Type:
    if (class(x) == "fREG") x = residuals(x)
    x = as.vector(x)

    # Call:
    call = match.call()

    # Test:
    x = sort(x)
    n = length(x)
    if (n < 5) stop("sample size must be greater than 4")
    p = pnorm((x - mean(x))/sqrt(var(x)))
    Dplus = max(seq(1:n)/n - p)
    Dminus = max(p - (seq(1:n) - 1)/n)
    K = max(Dplus, Dminus)
    if (n <= 100) {
        Kd = K
        nd = n
    } else {
        Kd = K * ((n/100)^0.49)
        nd = 100 }

    # P Value:
    PVAL = exp(-7.01256 * Kd^2 * (nd + 2.78019) + 2.99587 * Kd *
        sqrt(nd + 2.78019) - 0.122119 + 0.974598/sqrt(nd) + 1.67997/nd)
    if (PVAL > 0.1) {
        KK = (sqrt(n) - 0.01 + 0.85/sqrt(n)) * K
        if (KK <= 0.302) {
            PVAL = 1
        } else if (KK <= 0.5) {
            PVAL = 2.76773 - 19.828315 * KK + 80.709644 *
                KK^2 - 138.55152 * KK^3 + 81.218052 * KK^4
        } else if (KK <= 0.9) {
            PVAL = -4.901232 + 40.662806 * KK - 97.490286 *
                KK^2 + 94.029866 * KK^3 - 32.355711 * KK^4
        } else if (KK <= 1.31) {
            PVAL = 6.198765 - 19.558097 * KK + 23.186922 *
                KK^2 - 12.234627 * KK^3 + 2.423045 * KK^4
        } else {
            PVAL = 0
        }
    }
    names(PVAL) = ""

    # Result:
    test = list(
        statistic = c(D = K),
        p.value = PVAL,
        method = "Lilliefors Test", data.name = DNAME)

    # Add:
    if (is.null(title)) title = "Lilliefors (KS) Normality Test"
    if (is.null(description)) description = description()

    # Return Value:
    new("fHTEST",
        call = call,
        data = list(x = x),
        test = test,
        title = as.character(title),
        description = as.character(description) )
}


# ------------------------------------------------------------------------------


pchiTest <-
function(x, title = NULL, description = NULL)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Performs the Pearson chi-square normality test

    # Arguments:
    #   x - a numeric vector of data values.
    #   description - a brief description of the porject of type
    #       character.
    #   title - a character string which allows for a project title.

    # Details:
    #   n.classes - an integer value, the classes are build is such
    #       a way that they are equiprobable under the hypothesis
    #       of normality. We use Moore's default value.
    #   adjust - a logical flag,  if TRUE (default), the p-value
    #       is computed from a chi-square distribution with
    #       n.classes-3 degrees of freedom, otherwise from a
    #       chi-square distribution with n.classes-1 degrees of
    #       freedom. We print the results for both

    # Notes:
    #   A copy from contributed R-package 'nortest'
    #   Source:
    #       Package: nortest
    #       Title: Tests for Normality
    #       Version: 1.0
    #       Author: Juergen Gross
    #       Description: 5 omnibus tests for the composite hypothesis of normality
    #       Maintainer: Juergen Gross <gross@statistik.uni-dortmund.de>
    #       License: GPL version 2 or newer

    # FUNCTION:

    # Data Set Name:
    DNAME = deparse(substitute(x))

    # Convert Type:
    if (class(x) == "fREG") x = residuals(x)
    x = as.vector(x)

    # Call:
    call = match.call()

    # Moore's Default Value:
    n.classes = ceiling(2 * (length(x)^(2/5)))
    names(n.classes) = "Number of Classes"

    # Test:
    n = length(x)
    dfd = c(2, 0)
    names(dfd) = c("TRUE", "FALSE") # adjust
    num = floor(1 + n.classes * pnorm(x, mean(x), sqrt(var(x))))
    count = tabulate(num, n.classes)
    prob = rep(1/n.classes, n.classes)
    xpec = n * prob
    h = ((count - xpec)^2)/xpec
    P = sum(h)

    # For Splus compatibility:
    # pvalue = pchisq(P, n.classes - dfd - 1, lower.tail = FALSE)
    PVAL = c(1 - pchisq(P, n.classes - dfd[1] - 1),
        1 - pchisq(P, n.classes - dfd[2] - 1))
    names(PVAL) = c("Adhusted", "Not adjusted")

    # Return Value:
    test = list(
        statistic = c(P = P),
        p.value = PVAL,
        method = "Pearson Chi-Square Normality Test",
        data.name = DNAME,
        parameter = n.classes)

    # Add:
    if (is.null(title)) title = "Pearson Chi-Square Normality Test"
    if (is.null(description)) description = description()

    # Return Value:
    new("fHTEST",
        call = call,
        data = list(x = x),
        test = test,
        title = as.character(title),
        description = as.character(description) )
}


# ------------------------------------------------------------------------------


sfTest <-
function(x, title = NULL, description = NULL)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Performs the Shapiro-Francia normality test

    # Arguments:
    #   x - a numeric vector of data values.
    #   description - a brief description of the porject of type
    #       character.
    #   title - a character string which allows for a project title.

    # Notes:
    #   A copy from contributed R-package 'nortest'
    #   Source:
    #       Package: nortest
    #       Title: Tests for Normality
    #       Version: 1.0
    #       Author: Juergen Gross
    #       Description: 5 omnibus tests for the composite hypothesis of normality
    #       Maintainer: Juergen Gross <gross@statistik.uni-dortmund.de>
    #       License: GPL version 2 or newer

    # FUNCTION:

    # Data Set Name:
    DNAME = deparse(substitute(x))

    # Convert Type:
    if (class(x) == "fREG") x = residuals(x)
    x = as.vector(x)

    # Call:
    call = match.call()

    # Test:
    x = sort(x)
    n = length(x)
    if ((n < 5 || n > 5000)) stop("sample size must be between 5 and 5000")
    y = qnorm(ppoints(n, a = 3/8))
    W = cor(x, y)^2
    u = log(n)
    v = log(u)
    mu = -1.2725 + 1.0521 * (v - u)
    sig = 1.0308 - 0.26758 * (v + 2/u)
    z = (log(1 - W) - mu)/sig

    # For Splus compatibility:
    # pval = pnorm(z, lower.tail = FALSE)
    PVAL = 1 - pnorm(z)
    names(PVAL) = ""

    # Test:
    test = list(
        statistic = c(W = W),
        p.value = PVAL,
        method = "Shapiro-Francia Normality Test", data.name = DNAME)

    # Add:
    if (is.null(title)) title = "Shapiro - Francia Normality Test"
    if (is.null(description)) description = description()

    # Return Value:
    new("fHTEST",
        call = call,
        data = list(x = x),
        test = test,
        title = as.character(title),
        description = as.character(description))
}


################################################################################
# AUGMENTED FINITE SAMPLE JB TEST:


.jb.test <-
function(x)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   S3 version type finite sample adjusted JB test

    # Arguments:

    # Notes:
    #   S3 Version type of test.
    #   See also the Jarque-Bera test in Adrian Trapletti's
    #   contributed "tseries" R package.

    # FUNCTION:

    # Data Set Name:
    DNAME = deparse(substitute(x))

    # Check:
    if (NCOL(x) > 1) stop("x is not a vector or univariate time series")
    if (any(is.na(x))) stop("NAs in x")

    # LM Coefficients:
    n = length(x)
    c1lm = 6/n
    c2lm = 3
    c3lm = 24/n

    # ALM Coefficients:
    c1alm = 6*(n-2) / ((n+1)*(n+3))
    c2alm = 3*(n-1) / (n+1)
    c3alm = 24*n*(n-2)*(n-3) / ((n+1)^2*(n+3)*(n+5))

    # Statistic:
    m1 = sum(x)/n
    m2 = sum((x - m1)^2)/n
    m3 = sum((x - m1)^3)/n
    m4 = sum((x - m1)^4)/n
    b1 = (m3/m2^(3/2))
    b2 = (m4/(m2*m2))
    STATISTIC = b1*b1/c1lm + ((b2-c2lm)^2)/c3lm
    ALM = b1*b1/c1alm + ((b2-c2alm)^2)/c3alm
    names(STATISTIC) = "LM"

    # The rest goes as parameters ...
    pvalLM = 1 - .pjb(STATISTIC, N = n, type = "LM")
    pvalALM = 1 - .pjb(ALM, N = n, type = "ALM")
    PARAMETER = c(ALM, n, pvalLM, pvalALM)
    names(PARAMETER) = c( "ALM", "Sample Size", "LM p-value", "ALM p-value" )

    # P Values:
    PVAL = 1 - pchisq(STATISTIC[1], df = 2)
    names(PVAL) = ""

    # Method String:
    METHOD = "Jarque Bera Test"

    # Result:
    Digits = 3
    RVAL = list(
        statistic = round(STATISTIC, digits = Digits),
        parameter = round(PARAMETER, digits = Digits),
        p.value = round(PVAL, digits = Digits),
        method = METHOD,
        data.name = DNAME)

    # Return Value:
    class(RVAL) = "htest"
    RVAL
}


# ------------------------------------------------------------------------------


jbTest <-
function(x, title = NULL, description = NULL)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Performs finite sample adjusted Jarque-Bera LM and ALM test

    # Arguments:
    #   x - a univariate timeSeries object or numeric vector

    # Notes:
    #   S3 Version type of test.

    # FUNCTION:

    # Data Set Name:
    DNAME = deparse(substitute(x))

    # Call:
    call = match.call()

    # Convert Type:
    if (class(x) == "fREG") x = residuals(x)
    x = as.vector(x)

    # Test:
    test = .jb.test(x)
    class(test) = "list"
    parameter = test$parameter[2]
    statistic = c(test$statistic, test$parameter[1])
    p.value = c(test$parameter[3], test$parameter[4],
        Asymptotic = test$p.value[[1]])
    test$parameter = parameter
    test$statistic = statistic
    test$p.value = p.value

    # Add Title:
    if (is.null(title)) title = "Jarque - Bera Normality Test"
    if (is.null(description)) description = description()

    # Return Value:
    new("fHTEST",
        call = call,
        data = list(x = x),
        test = test,
        title = as.character(title),
        description = as.character(description) )
}


################################################################################
