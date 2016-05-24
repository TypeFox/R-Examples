gx.summary.groups <-
function (group, x, xname = deparse(substitute(x)), log = FALSE) 
{
    groupname <- deparse(substitute(group))
    group.stats <- tapply(x, group, gx.summary, log = log, iftell = FALSE)
    if (log) 
        cat("  Data log10 transformed: SD, CV% and SE in log10 units\n")
    cat("  Summary Stats for", xname, "subset by", groupname, 
        "\n\t", "N NAs - Min Q1 M Q2 Max - MAD IQR_SD - Mean SD CV% - SE 95% CI on Mean\n")
    nstats <- length(group.stats)
    for (i in 1:nstats) {
        gi <- names(group.stats[i])
        ii <- unlist(group.stats[i], use.names = FALSE)
        cat("\n  ", gi, ":\t ", ii[1], " ", ii[2], " - ", ii[3], 
            " ", ii[4], " ", ii[5], " ", ii[6], " ", ii[7], " - ", 
            ii[8], " ", ii[9], " - ", ii[10], " ", ii[11], " ", 
            ii[12], "% - ", ii[13], " ", ii[14], " <-> ", ii[15], 
            "\n", sep = "")
    }
    cat("\n")
    invisible()
}
