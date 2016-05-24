# in this demo, we'll load an OTU table, pairs table, and quads table
# that are included with the package. this approach fits multiple samples
# and packages them up with the least amount of typing.

# load in the example data. in the source, it's located in inst/extdata
# we use the system.file command here to locate the package however it's
# been installed on your machine.
otu.table.filename <- system.file('demo', 'otus.dat', package='texmexseq')
pairs.table.filename <- system.file('demo', 'pairs', package='texmexseq')
quads.table.filename <- system.file('demo', 'quads', package='texmexseq')

# load in an experiment, consisting of 8 samples in 4 pairs in 2 quads
my.experiment <- Experiment(otu.table.filename, pairs.table.filename, quads.table.filename)

# it's all done!
# now we can plot the distribution of counts in an individual sample
hist(my.experiment$quads$inoculum1$control$sample0$n)

# visualize the quality of the fit of the Poisson lognormal distribution
# to counts in this sample
# fit0 = the fit for the first sample in the quad
FitPpPlot(my.experiment$quads$inoculum1$control$fit0)

# i'm interested in this inoculum1 quad, so i'll just define...
my.quad <- my.experiment$quads$inoculum1

# or, we could show how the counts for every OTU compare in two samples
# (i.e., a pair) before any Poisson lognormal fitting is applied
PairPlot(my.quad$control)

# how do the rescaled reads for OTUs in a quad behave?
QuadPlot(my.quad)

# or, the same quad, but now showing the change in F (cdf) rather than
# the change in rescaled reads
QuadPlot(my.quad, dF=TRUE)

# hmm.. that might be a little tricky to see, since there are a lot of
# points right on top of one another. let's add a little noise to get
# a better sense of the point density
QuadPlot(my.quad, dF=TRUE, jitter.amount=0.1)
