nominalTest <- function(tab, limit.exp = 5){

    c <- suppressWarnings(chisq.test(tab))
    f <- fisher.test(tab, simulate.p.value = TRUE, B = 10000)
    
    p <- c$p.value
    test <- expression(Chi^2*"-test")
    if (min(c$expected) <= limit.exp){
        p <- f$p.value
        test <- "Fisher's exact test"}    
    return(list("p" = p, "test" = test))
}
