`plot.stcs` <-
function(x, stat=1:4, type=c("hist", "rank"), trafo=c("none", "log",
"ratio"), show=TRUE, ylab, xlab, ...)
{
plot.mefa(mefa(x), stat, type, trafo, show, ylab, xlab, ...)
}
