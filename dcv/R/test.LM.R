test.LM <-
function(x){
    linmod <- lm(x[,1] ~ x[,2])
    summary(linmod)
}

