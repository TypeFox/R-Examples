altfavstats <- function(x) {
    cat(paste("  mean:", format(mean(x),4),"\n"))
    cat(paste(" edian:", format(median(x),4),"\n"))
    cat(paste("    sd:", format(sd(x),4),"\n"))
}
altfavstats((1:20)^2)
