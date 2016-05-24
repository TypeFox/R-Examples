exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
featFile <- "gencode.ENr334.gp"
unzip(exampleArchive, featFile)
f <- read.feat(featFile)
f <- add.signals.feat(f)
# let's just look at one gene
geneNames <- tagval.feat(f, "transcript_id")
f <- f[geneNames==geneNames[1],]

# This features file already is correct, so let's mess it up to see
# how fix.start.stop can fix it:

#modify first CDS to not include start
startCodon <- f[f$feature=="start_codon",]
firstCds <- which(f$feature=="CDS" & f$start==startCodon$start)
f[firstCds,]$start <- startCodon$end+1
#modify last CDS to include stop
stopCodon <- f[f$feature=="stop_codon",]
lastCds <- which(f$feature=="CDS" & f$end+1==stopCodon$start)
f[lastCds,]$end <- stopCodon$end
# now call fix.start.stop to fix
f.fixed <- fix.start.stop.feat(f)

# first CDS has been fixed to include start codon
f[firstCds,]
f.fixed[firstCds,]
# last CDS has been fixed to not include stop codon
f[lastCds,]
f.fixed[lastCds,]

unlink(featFile)
