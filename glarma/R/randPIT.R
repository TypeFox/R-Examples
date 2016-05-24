### Function to randomize PITs
normRandPIT <- function(object){
    temp <- glarmaPredProb(object)
    rt <- qnorm(runif(length(temp$lower), temp$lower, temp$upper))
    rtMid <- qnorm((temp$lower + temp$upper)/2)
    list(rt = rt, rtMid = rtMid)
}

