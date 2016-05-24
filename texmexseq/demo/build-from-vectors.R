# in this demo, we'll build up a quad and plot it. this is a little slower
# than the approach used in the other demo (which uses the Experiment command),
# but it will show how the objects are related to one another.

# i'll start by loading in the example OTU data. you could imagine that
# you somehow loaded in your own count data and were interested in the
# results.
otu.table.filename <- system.file('demo', 'otus.dat', package='texmexseq')
otu.table <- ReadOtuTable(otu.table.filename)
inoculum1.control.before <- otu.table$inoculum1.control.before 
inoculum1.control.after <- otu.table$inoculum1.control.after 
inoculum1.treatment.before <- otu.table$inoculum1.treatment.before 
inoculum1.treatment.after <- otu.table$inoculum1.treatment.after 

# now we have four vectors that show the count data from a set of
# four samples. there are two timepoints from the control unit and
# two timepoints from the experiment (aka "treatment") unit.
# first, we'll make these count vectors into Sample objects
control.before.sample <- Sample(inoculum1.control.before, name="control before")
control.after.sample <- Sample(inoculum1.control.after, name="control after")
treatment.before.sample <- Sample(inoculum1.treatment.before, name="treatment before")
treatment.after.sample <- Sample(inoculum1.treatment.after, name="treatment after")

# now we'll join up the samples into pairs
control.pair <- SamplePair(control.before.sample, control.after.sample, name="control")
treatment.pair <- SamplePair(treatment.before.sample, treatment.after.sample, name="treatment")

# and then link those two pairs into a quad
my.quad <- SampleQuad(control.pair, treatment.pair, name="my quad")

# from here we can plot it or whatever else, which is shown in the
# other demo
