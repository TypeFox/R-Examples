
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA

# Copyrights (C)
# for this R-port:
#   1999 - 2008, Diethelm Wuertz, Rmetrics Foundation, GPL
#   Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
#   info@rmetrics.org
#   www.rmetrics.org
# for the code accessed (or partly included) from other R-ports:
#   see R's copyright and license files
# for the code accessed (or partly included) from contributed R-ports
# and other sources
#   see Rmetrics's copyright file


################################################################################
# FUNCTION:                 DESCRIPTION:
#  summary                   Summary method for an object of class 'fGARCH'
################################################################################


setMethod(f = "summary", signature(object = "fGARCH"),
          function(object)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Summary method for an object of class "fGARCH"

    # Arguments:
    #   object - an object of class 'fGARCH'

    # FUNCTION:

    # same output as show method
    show(object)

    # Lagged Series:
    .tslagGarch = function (x, k = 1) {
        ans = NULL
        for (i in k) ans = cbind(ans, .tslag1Garch(x, i))
        indexes = (1:length(ans[, 1]))[!is.na(apply(ans, 1, sum))]
        ans = ans[indexes, ]
        if (length(k) == 1) ans = as.vector(ans)
        ans }
    .tslag1Garch = function (x, k) {
        c(rep(NA, times = k), x[1:(length(x) - k)]) }

    # Statistical Tests:
    cat("\nStandardised Residuals Tests:\n")
    r.s = object@residuals/object@sigma.t
    ans = NULL
    # Normality Tests:
    jbtest = jarqueberaTest(r.s)@test
    ans = rbind(ans, c(jbtest[1], jbtest[2]))
    if (length(r.s) < 5000) {
        swtest = shapiro.test(r.s)
        if (swtest[2] < 2.6e-16) swtest[2] = 0
        ans = rbind(ans, c(swtest[1], swtest[2]))
    } else {
        ans = rbind(ans, c(NA, NA))
    }
    # Ljung-Box Tests:
    box10 = Box.test(r.s, lag = 10, type = "Ljung-Box")
    box15 = Box.test(r.s, lag = 15, type = "Ljung-Box")
    box20 = Box.test(r.s, lag = 20, type = "Ljung-Box")
    ans = rbind(ans, c(box10[1], box10[3]))
    ans = rbind(ans, c(box15[1], box15[3]))
    ans = rbind(ans, c(box20[1], box20[3]))
    box10 = Box.test(r.s^2, lag = 10, type = "Ljung-Box")
    box15 = Box.test(r.s^2, lag = 15, type = "Ljung-Box")
    box20 = Box.test(r.s^2, lag = 20, type = "Ljung-Box")
    ans = rbind(ans, c(box10[1], box10[3]))
    ans = rbind(ans, c(box15[1], box15[3]))
    ans = rbind(ans, c(box20[1], box20[3]))
    # Ljung-Box Tests - tslag required
    lag.n = 12
    x.s = as.matrix(r.s)^2
    n = nrow(x.s)
    tmp.x = .tslagGarch(x.s[, 1], 1:lag.n)
    tmp.y = x.s[(lag.n + 1):n, 1]
    fit = lm(tmp.y ~ tmp.x)
    stat = (n-lag.n) * summary.lm(fit)$r.squared
    ans = rbind(ans, c(stat, p.value = 1 - pchisq(stat, lag.n)) )
    # Add Names:
    rownames(ans) = c(
        " Jarque-Bera Test   R    Chi^2 ",
        " Shapiro-Wilk Test  R    W     ",
        " Ljung-Box Test     R    Q(10) ",
        " Ljung-Box Test     R    Q(15) ",
        " Ljung-Box Test     R    Q(20) ",
        " Ljung-Box Test     R^2  Q(10) ",
        " Ljung-Box Test     R^2  Q(15) ",
        " Ljung-Box Test     R^2  Q(20) ",
        " LM Arch Test       R    TR^2  ")
    colnames(ans) = c("Statistic", "p-Value")
    print(ans)

    # Information Criterion Statistics:
    cat("\nInformation Criterion Statistics:\n")
    print(object@fit$ics)

    # Return Value:
    cat("\n")
    invisible()
})


################################################################################

