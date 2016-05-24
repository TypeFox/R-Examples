exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
files <- c("ENr334-100k.maf", "rev.mod", "gencode.ENr334-100k.gff")
unzip(exampleArchive, files)
# make "conserved" and "neutral" models and a phylo-HMM that describes
# transitions between them, and predict conserved elements (this
# is the same thing that phastCons does, but can be extended to general
# phylo-HMMs)
align <- read.msa("ENr334-100k.maf")
neutralMod <- read.tm("rev.mod")

# create a conserved model
conservedMod <- neutralMod
conservedMod$tree <- rescale.tree(neutralMod$tree, 0.3)

# create a simple phylo-HMM
state.names <- c("neutral", "conserved")
h <- hmm(matrix(c(0.99, 0.01, 0.01, 0.99), nrow=2,
                dimnames=list(state.names, state.names)),
                eq.freq=c(neutral=0.9, conserved=0.1))
scores <- score.hmm(align, mod=list(neutral=neutralMod,
                             conserved=conservedMod),
                    hmm=h, states="conserved")
# try an alternate approach of comparing likelihoods of genes 
feats <- read.feat("gencode.ENr334-100k.gff")
# plot in a region with some genes
plot.track(list(as.track.feat(scores$in.states, name="hmmScores"),
                as.track.feat(feats[feats$feature=="CDS",], name="genes")),
           xlim=c(41650000, 41680000))
unlink(files)
