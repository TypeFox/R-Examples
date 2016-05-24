cramer <- function (tbl) 
{
#
# Compute Cramer's V from a two-way table.
#
if (!is.matrix (tbl)) stop ("Two-way table needed")
k <- min (nrow(tbl), ncol(tbl))
# Turn off warnings for the moment; revert to current setting when done
current.warning.level <- options ()$warn
on.exit (options(warn = current.warning.level))
options(warn = -1) # no warnings at all!
#
chisq <- chisq.test (tbl)$statistic
n <- sum (tbl)
return (unname (sqrt (chisq / (n * (k - 1)))))
}

