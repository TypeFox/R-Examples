
brunner.munzel.test<-function (x, y, alternative = c("two.sided", "greater", "less"), 
    alpha = 0.05) 
{
    alternative <- match.arg(alternative)
    DNAME = paste(deparse(substitute(x)), "and", deparse(substitute(y)))
    x <- na.omit(x)
    y <- na.omit(y)
    n1 = length(x)
    n2 = length(y)
    r1 = rank(x)
    r2 = rank(y)
    r = rank(c(x, y))
    m1 = mean(r[1:n1])
    m2 = mean(r[n1 + 1:n2])
    pst = (m2 - (n2 + 1)/2)/n1
    v1 = sum((r[1:n1] - r1 - m1 + (n1 + 1)/2)^2)/(n1 - 1)
    v2 = sum((r[n1 + 1:n2] - r2 - m2 + (n2 + 1)/2)^2)/(n2 - 1)
    statistic = n1 * n2 * (m2 - m1)/(n1 + n2)/sqrt(n1 * v1 + 
        n2 * v2)
    dfbm = ((n1 * v1 + n2 * v2)^2)/(((n1 * v1)^2)/(n1 - 1) + 
        ((n2 * v2)^2)/(n2 - 1))
    if ((alternative == "greater") | (alternative == "g")) {
        p.value = pt(statistic, dfbm)
    }
    else if ((alternative == "less") | (alternative == "l")) {
        p.value = 1-pt(statistic, dfbm)
    }
    else {
        alternative = "two.sided"
        p.value = 2 * min(pt(abs(statistic), dfbm), (1 - pt(abs(statistic), 
            dfbm)))
    }
    conf.int = c(pst - qt(1 - alpha/2, dfbm) * sqrt(v1/(n1 * 
        n2^2) + v2/(n2 * n1^2)), pst + qt(1 - alpha/2, dfbm) * 
        sqrt(v1/(n1 * n2^2) + v2/(n2 * n1^2)))
    estimate = pst
    ESTIMATE = pst
    names(ESTIMATE) = "P(X<Y)+.5*P(X=Y)"
    STATISTIC = statistic
    names(STATISTIC) = "Brunner-Munzel Test Statistic"
    PARAMETER = dfbm
    names(PARAMETER) = "df"
    CONF.INT = conf.int
    names(CONF.INT) = c("lower", "upper")
    attr(CONF.INT, "conf.level") = (1 - alpha)
    METHOD = "Brunner-Munzel Test"
    structure(list(estimate = ESTIMATE, conf.int = CONF.INT, 
        statistic = STATISTIC, parameter = PARAMETER, p.value = p.value, 
        method = METHOD, data.name = DNAME), class = "htest")
}