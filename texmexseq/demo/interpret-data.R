# in this demo, we'll load in the example data and examine it from
# a more scientific point of view. maybe do the load-from-tables demo first.

# load in the example data.
otu.table.filename <- system.file('demo', 'otus.dat', package='texmexseq')
pairs.table.filename <- system.file('demo', 'pairs', package='texmexseq')
quads.table.filename <- system.file('demo', 'quads', package='texmexseq')
my.experiment <- Experiment(otu.table.filename, pairs.table.filename, quads.table.filename)

# let's compare the two quads
inoc1 <- my.experiment$quads$inoculum1
inoc2 <- my.experiment$quads$inoculum2

# it looks like, in inoculum 1, the control unit went through community
# composition changes that were mostly unrelated to the effects the
# experimental unit went through: the points in the middle of the plot look
# mostly uncorrelated.
PlotQuad(inoc1)

# in inoculum 2, however, the control and experimental unit underwent similar
# changes, suggesting that the treatment for inoculum 2 had a weak effeck: it
# didn't do much beyond what the control did.
PlotQuad(inoc2)

# maybe we're interested in OTUs that bloomed in inoculum 2. to pick some
# OTUs that bloomed in the treatment but not the control, maybe I would pick
otus <- inoc2$control$dz < 1 & inoc2$treatment$dz > 1

# unfortunately, lots of the dz are NA, and R gives NA when trying to compare
# NA > 1. so let's replace all NA entries in otus with FALSE
otus <- replace(otus, is.na(otus), FALSE)

# now we can see which OTUs those are
PlotQuad(inoc2, highlight=otus)

# most of them are on the top axis, meaning that their dz is +inf, which means
# that they went from 0 counts in the before timepoint to more than 0 in the
# after timepoint. let's look at which OTUs these are and what their relative
# abundances were in the before sample:
head(inoc2$treatment$sample0$ra[otus], n=20)

# and in the after sample
head(inoc2$treatment$sample1$ra[otus], n=20)

# you can ask which rows in the OTU table those OTUs were in
which(otus)

# or, more directly, ask which OTUs those are
rownames(my.experiment$otu.table)[otus]