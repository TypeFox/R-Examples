qbinom(0.975,200,0.5)
qbinom(0.025,200,0.5)
pbinom(85:86,200,0.5)
1 - pbinom(114:115,200,0.5)

# define a function to calculate power for given sample size.
power = function(size,null=0.50,alt=0.55){ 
    leftCritical  <- -1 + qbinom(0.025, round(size),null)
    rightCritical <-  1 + qbinom(0.975, round(size),null)
    leftPower <- pbinom( leftCritical, round(size),alt)
    rightPower <- 1 - pbinom( rightCritical - 1, round(size),alt)

    result <- leftPower + rightPower

    # the rest of this allows for more verbose output
    attr(result,"power") = leftPower + rightPower
    attr(result,"null") = null
    attr(result,"alt") = alt
    attr(result,"leftCritical") = leftCritical
    attr(result,"rightCritical") = rightCritical
    attr(result,"leftPower") = leftPower
    attr(result,"rightPower") = rightPower

    return (result)
}

# just the power values:
as.numeric(power(200))
as.numeric(power(400))

# more verbose output:
power(200)

power(400)

# find sample size with 90% power
uniroot(function(size){ as.numeric(power(size)) -0.90} ,c(400,5000))
