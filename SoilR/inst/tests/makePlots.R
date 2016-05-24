#
# vim:set ff=unix expandtab ts=2 sw=2:
#!/usr/bin/Rscript
require(RUnit)
source("prolog.R")
require(mvbutils)
pdf(file="figure.pdf",height=13.5, width=25)
foodweb(charlim=100)
dev.off()

