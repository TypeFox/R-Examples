# Function for power calculations
#
# Sample size for multiple groups
sample.size <- function(my.frame, sig.level, power) {
    sample.calculations <- sapply(1:length(my.frame$variable), FUN = function(x) {
        pwr::pwr.f2.test(u = my.frame$df[x], f2 = my.frame$f2[x], sig.level = sig.level, power = power)
    }
    )
    required.sample <- sapply(1:length(my.frame$variable), FUN = function(x) {
        c(sample.calculations[,x]$v + max(my.frame$lev)) #changed to max
    }
    )
    return(ceiling(required.sample))
}
