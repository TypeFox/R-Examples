

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Output consistency of unify Ranges:
# Consecutive ranges must not ovelap.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

urg <- unifyRanges(enex)
clevels <- levels(urg@ev$gtf$seqid)
n <- length(clevels)

for(i in 1:n)
{
    xtr <- urg@ev$gtf[urg@ev$gtf$seqid == clevels[i], ]
    xb <- c(xtr$begin, 0)
    xe <- c(0,xtr$end)
    diff <- xb - xe
    if(table(diff > 0)[1] > 1)
        stop("[unifyRanges] inconsistent output!")
}

