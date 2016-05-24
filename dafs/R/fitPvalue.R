fitPvalue = function(fit){
    Fstats = summary(fit)$fstatistic
    return(1-pf(Fstats[1],Fstats[2],Fstats[3]))
}
