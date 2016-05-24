
####  this demo shows a short workflow that could lead to an elaborate analysis.

####  first we create a file of metagenome IDs.
####  a user would have their own, in a real application.

writeLines (colnames (xx1), "my_ids.txt")
readLines ("my_ids.txt")

####  annotation data is easy to retrieve.

zz <- biomRequest (file="my_ids.txt", group_level="level2", evalue=1)

####  a log (or other) transformation is often a good idea.

zz0 <- transform (zz, t_Log)

####  metadata is easy to inspect, as follows.

columns (zz0, "host_common_name|samp_store_temp|material")

####  metadata of interest can then be incorporated into a plot.

princomp (zz0, map=c(col="host_common_name", pch="samp_store_temp"), labels="$$pubmed_id", cex=2)

####  next, a grouped distance calculation.

distx (zz0, groups="$$host_common_name")

####  annotations providing significant group differentitation can be identified.

pp <- (rowstats (zz0, groups="$$material") $ p.value < 0.05)
pp[is.na(pp)] <- FALSE
pp

####  that information can be used to make an informative heatmap.

image (zz0 [pp,], margins=c(5,10), cexRow=0.3)

####  for comparison, here is the same heatmap, but including all annotations.

image (zz0, margins=c(5,10), cexRow=0.3)
