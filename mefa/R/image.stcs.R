`image.stcs` <-
function(x, segm=NULL, trafo=c("none", "log", "bins", "prab"), 
probs = seq(0, 1, 0.05), ordering=TRUE, reverse=TRUE, names = FALSE,
show=TRUE, ylab, xlab, ...)
{
image.mefa(mefa(x), segm, trafo, probs, ordering, reverse, show, names, ylab, xlab, ...)
}
