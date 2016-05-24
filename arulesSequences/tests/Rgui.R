
### ceeboo 2012

if (.Platform$OS  == "windows" && 
    .Platform$GUI == "RTerm") {

    library("arulesSequences")

    ## use example data from paper
    data(zaki)
    ## mine frequent sequences
    s1 <- cspade(zaki, parameter = list(support = 0.4)) 

    ## fake
    .Platform$GUI <- "Rgui"

    s2 <- cspade(zaki, parameter = list(support = 0.4))

    stopifnot(identical(s1, s2))
}

###
