## =============================
## Plotting individual sequences
## =============================

seqiplot <- function(seqdata, group=NULL, title=NULL, ...) {
	seqplot(seqdata, group=group, type="i", title=title, ...)
}
seqIplot <- function(seqdata, group=NULL, title=NULL, ...) {
	seqplot(seqdata, group=group, type="I", title=title, ...)
}
