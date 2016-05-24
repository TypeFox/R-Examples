"plot.relimplmbootMI" <- 
function (x, ..., lev = max(x@level), names.abbrev = 4, ylim=NULL, main=NULL, cex.title=1.5) 
{
    # function shows barplots with error bars indicating confidence interval for chosen level
    # if chosen level not available, confidence interval for largest available level is produced
    # this function is created for S3-type referencing
plot.relimplmbooteval(x, ..., lev = lev, names.abbrev = names.abbrev, ylim=ylim, main=main, cex.title=cex.title)
}
