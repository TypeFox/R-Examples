# calculates prob of shared birthday
birthdayprob <- function(n) { 
        last <- 366 - n
        1 - ( prod(seq(last,365)) / 365^n )
}

birthdayprob(10)
cbind(20:25,sapply(20:25,birthdayprob))
