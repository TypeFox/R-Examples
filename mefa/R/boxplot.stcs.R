`boxplot.stcs` <-
function(x, stat=1:4, all = TRUE, show=TRUE, ylab, xlab, ...)
{
boxplot.mefa(mefa(x), stat, all, show, ylab, xlab, ...)
}

