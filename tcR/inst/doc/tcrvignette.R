## ----eval=TRUE,echo=FALSE,warning=FALSE,message=FALSE--------------------
library(tcR)
data(twa)
data(twb)

## ----eval=FALSE,echo=TRUE------------------------------------------------
#  # Load the package and load the data.
#  library(tcR)
#  data(twa)  # "twa" - list of length 4
#  data(twb)  # "twb" - list of length 4
#  
#  # Explore the data.
#  head(twa[[1]])
#  head(twb[[1]])

## ----eval=FALSE,echo=TRUE------------------------------------------------
#  # Run help to see available alphabets.
#  ?genealphabets
#  ?genesegments
#  data(genesegments)

## ----eval=FALSE,echo=TRUE------------------------------------------------
#  # Parse file in "~/mitcr/immdata1.txt" as a MiTCR file.
#  immdata1 <- parse.file("~/mitcr_data/immdata1.txt", 'mitcr')
#  # equivalent to
#  immdata1.eq <- parse.mitcr("~/mitcr_data/immdata1.txt")
#  
#  # Parse folder with MiGEC files.
#  immdata <- parse.folder("~/migec_data/", 'migec')

## ----eval=TRUE, echo=TRUE------------------------------------------------
# No D genes is available here hence "" at "D.genes" and "-1" at positions.
str(twa[[1]])

str(twb[[1]])

## ----eval=TRUE,echo=TRUE-------------------------------------------------
cloneset.stats(twb)

## ----eval=TRUE,echo=TRUE-------------------------------------------------
repseq.stats(twb)

## ----eval=TRUE,echo=TRUE-------------------------------------------------
                            # How many clonotypes fill up approximately
clonal.proportion(twb, 25)  # the 25% of the sum of values in 'Read.count'?

## ----echo=TRUE, eval=TRUE, fig=TRUE, fig.height=4, fig.width=5.5, message=FALSE, fig.align='center'----
                          # What accounts a proportion of the top-10 clonotypes' reads
top.proportion(twb, 10)   # to the overall number of reads?
vis.top.proportions(twb)  # Plot this proportions.

## ----eval=TRUE,echo=TRUE-------------------------------------------------
                                # What is a proportion of sequences which
                                # have 'Read.count' <= 100 to the
tailbound.proportion(twb, 100)  # overall number of reads?

## ----eval=TRUE, echo=TRUE, fig.height=4, fig.width=6.5, fig.align='center'----
# data(twb)
# Compute summary space of clones, that occupy
# [0, .05) and [.05, 1] proportion.
clonal.space.homeostasis(twb, c(Low = .05, High = 1))
# Use default arguments:
clonal.space.homeostasis(twb[[1]])

twb.space <- clonal.space.homeostasis(twb)
vis.clonal.space(twb.space)

## ----eval=TRUE,echo=TRUE-------------------------------------------------
imm.in <- get.inframes(twb) # Return all in-frame sequences from the 'twb'.

                            # Count the number of out-of-frame sequences
count.outframes(twb, 5000)  # from the first 5000 sequences.

## ----eval=TRUE,echo=TRUE-------------------------------------------------
imm.in <- get.frames(twb, 'in') # Similar to 'get.inframes(twb)'.

count.frames(twb[[1]], 'all')   # Just return number of rows.

flag <- 'out'
count.frames(twb, flag, 5000)   # Similar to 'count.outframes(twb, 5000)'.

## ----eval=TRUE,echo=TRUE-------------------------------------------------
cmv <- data.frame(CDR3.amino.acid.sequence = c('CASSSANYGYTF', 'CSVGRAQNEQFF', 'CASSLTGNTEAFF', 'CASSALGGAGTGELFF', 'CASSLIGVSSYNEQFF'),
                  V.genes = c('TRBV4-1', 'TRBV4-1', 'TRBV4-1', 'TRBV4-1', 'TRBV4-1'), stringsAsFactors = F)

cmv

## ----eval=TRUE,echo=TRUE-------------------------------------------------
twb <- set.rank(twb)
# Case 1.
cmv.imm.ex <- 
  find.clonotypes(.data = twb[1:2], .targets = cmv[,1], .method = 'exact',
                  .col.name = c('Read.count', 'Total.insertions'),
                  .verbose = F)
head(cmv.imm.ex)

# Case 2.
# Search for CDR3 sequences with hamming distance <= 1
# to the one of the cmv$CDR3.amino.acid.sequence with
# matching V genes. Return ranks of found sequences.
cmv.imm.hamm.v <- 
  find.clonotypes(twb[1:3], cmv, 'hamm', 'Rank', 
                  .target.col = c('CDR3.amino.acid.sequence',
                                  'V.gene'),
                  .verbose = F)
head(cmv.imm.hamm.v)

# Case 3.
# Similar to the previous example, except
# using levenshtein distance and the "Read.count" column.
cmv.imm.lev.v <- 
  find.clonotypes(twb[1:3], cmv, 'lev', 
                  .target.col = c('CDR3.amino.acid.sequence', 'V.gene'),
                  .verbose = F)
head(cmv.imm.lev.v)

## ----eval=TRUE,echo=TRUE-------------------------------------------------
imm1.vs <- geneUsage(twb[[1]], HUMAN_TRBV)
head(imm1.vs)

imm.vs.all <- geneUsage(twb, HUMAN_TRBV)
imm.vs.all[1:10, 1:4]

# Compute joint V-J counts
imm1.vj <- geneUsage(twb[[1]], list(HUMAN_TRBV, HUMAN_TRBJ))
imm1.vj[1:5, 1:5]

## ----eval=TRUE, echo=TRUE, message=FALSE, fig.align='center', fig.height=5, fig.width=7----
# Put ".dodge = F" to get distinct plot for every data frame in the given list.
vis.gene.usage(twb, HUMAN_TRBJ, .main = 'twb J-usage dodge', .dodge = T)

## ----eval=TRUE, echo=TRUE, message=FALSE, fig.align='center', fig.height=6, fig.width=9----
vis.gene.usage(twb, HUMAN_TRBJ, .main = 'twb J-usage column', .dodge = F, .ncol = 2)

## ----eval=TRUE, echo=TRUE, message=FALSE, fig.align='center', fig.height=5, fig.width=7----
vis.gene.usage(imm1.vs, NA, .main = 'twb[[1]] V-usage', .coord.flip = F)

## ----eval=T, echo=TRUE, fig.align='center'-------------------------------
                              # Transform "0:100" to distribution with Laplace correction 
entropy(0:100, .laplace = 1)  # (i.e., add "1" to every value before transformation).
                              
entropy.seg(twb, HUMAN_TRBV)  # Compute entropy of V-segment usage for each data frame.
                  
js.div.seg(twb[1:2], HUMAN_TRBV, .verbose = F)
imm.js <- js.div.seg(twb, HUMAN_TRBV, .verbose = F) 
vis.radarlike(imm.js, .ncol = 2)

## ----eval=TRUE, echo=TRUE, fig.align='center', fig.height=4.5, fig.width=6----
pca.segments(twb, .genes = HUMAN_TRBV)  # Plot PCA results of V-segment usage.

# Return object of class "prcomp"
class(pca.segments(twb, .do.plot = F, .genes = HUMAN_TRBV))

## ----eval=TRUE, echo=T, fig.align='center', warning=FALSE----------------
# Equivalent to intersect(twb[[1]]$CDR3.nucleotide.sequence,
#                         twb[[2]]$CDR3.nucleotide.sequence)
repOverlap(twb[1:2], 'exact', 'nuc', .verbose = F)

# Equivalent to intersectClonesets(twb, "n0e", .norm = T)
repOverlap(twb, 'exact', 'nuc', .norm = T, .verbose = F)
# Intersect by amino acid clonotypes + V genes
repOverlap(twb, 'exact', 'aa', .vgene = T, .verbose = F)

# Plot a heatmap of the number of shared clonotypes.
vis.heatmap(repOverlap(twb, 'exact', 'aa', .vgene = T, .verbose = F), .title = 'twb - (ave)-intersection', .labs = '')

## ----eval=TRUE, echo=TRUE------------------------------------------------
# Get logic vector of shared elements, where
# elements are tuples of CDR3 nucleotide sequence and corresponding V-segment
imm.1.2 <- intersectLogic(twb[[1]], twb[[2]],
                           .col = c('CDR3.amino.acid.sequence', 'V.gene'))  
# Get elements which are in both twb[[1]] and twb[[2]].
head(twb[[1]][imm.1.2, c('CDR3.amino.acid.sequence', 'V.gene')])

## ----eval=TRUE, echo=T, fig.align='center', fig.height=6.5, fig.width=10, warning=FALSE----
twb.top <- top.cross(.data = twb, .n = seq(500, 10000, 500), .verbose = F, .norm = T)
top.cross.plot(twb.top)

## ----eval=TRUE, echo=TRUE, results='hold'--------------------------------
# Apply the Morisitas overlap index to the each pair of repertoires.
# Use information about V genes (i.e. one CDR3 clonotype is equal to another
# if and only if their CDR3 aa sequences are equal and their V genes are equal)
repOverlap(twb, 'morisita', 'aa', 'read.count', .vgene = T, .verbose = F)

## ----eval=TRUE, echo=TRUE------------------------------------------------
# Compute shared repertoire of amino acid CDR3 sequences and V genes
# which has been found in two or more people and return the Read.count column
# of such clonotypes from each data frame in the input list.
imm.shared <- shared.repertoire(.data = twb, .type = 'avrc', .min.ppl = 2, .verbose = F)
head(imm.shared)
shared.representation(imm.shared)  # Number of shared sequences.

## ----eval=TRUE, echo=TRUE, results='hold'--------------------------------
# Evaluate the diversity of clones by the ecological diversity index.
repDiversity(twb, 'div', 'read.count')
sapply(twb, function (x) diversity(x$Read.count))

## ----eval=TRUE, echo=TRUE, results='hold'--------------------------------
# Compute the diversity as the inverse probability of choosing two similar clonotypes.
repDiversity(twb, 'inv.simp', 'read.prop')
sapply(twb, function (x) inverse.simpson(x$Read.proportion))

## ----eval=TRUE, echo=TRUE, results='hold'--------------------------------
# Evaluate the skewness of clonal distribution.
repDiversity(twb, 'gini.simp', 'read.prop')
sapply(twb, function (x) gini.simpson(x$Read.proportion))

## ----eval=TRUE, echo=TRUE, results='hold'--------------------------------
# Compute diversity of repertoire using Chao index.
repDiversity(twb, 'chao1', 'read.count')
sapply(twb, function (x) chao1(x$Read.count))

## ----eval=TRUE, echo=TRUE, fig.height=4, fig.width=5.5, fig.align='center'----
vis.count.len(twb[[1]], .name = "twb[[1]] CDR3 lengths", 
              .col = "Read.count")

## ----eval=TRUE, echo=TRUE, fig.height=4, fig.width=5.5, fig.align='center', warning=FALSE, message=FALSE----
# I comment this to avoid a strange bug in ggplot2. Will uncomment later.
# vis.number.count(twb[[1]], .name = "twb[[1]] count distribution")

## ----echo=TRUE, eval=TRUE, fig.height=4, fig.width=5.5, message=FALSE, fig.align='center'----
vis.top.proportions(twb, c(10, 500, 3000, 10000), .col = "Read.count")

## ----eval=TRUE, echo=TRUE, fig.height=4, fig.width=6.5, fig.align='center'----
twb.space <- clonal.space.homeostasis(twb)
vis.clonal.space(twb.space)

## ----eval=TRUE, echo=TRUE, fig.align='center', warning=FALSE, message=FALSE----
twb.shared <- repOverlap(twb, "exact", .norm = F, .verbose = F)
vis.heatmap(twb.shared, .title = "Twins shared nuc clonotypes", 
            .labs = c("Sample in x", "Sample in y"), .legend = "# clonotypes")

## ----eval=T, echo=TRUE, fig.align='center'-------------------------------
twb.js <- js.div.seg(twb, HUMAN_TRBV, .verbose = F) 
vis.radarlike(twb.js, .ncol = 2)

## ----eval=TRUE, echo=TRUE, message=FALSE, fig.align='center', fig.height=5, fig.width=7----
vis.gene.usage(twb[[1]], HUMAN_TRBV, .main = 'Sample I V-usage')

## ----eval=TRUE, echo=TRUE, message=FALSE, fig.align='center', fig.height=7, fig.width=5----
vis.gene.usage(twb[[2]], HUMAN_TRBV, .main = 'Sample II V-usage', .coord.flip = T)

## ----eval=TRUE, echo=TRUE, message=FALSE, fig.align='center', fig.height=5, fig.width=7----
twb.jusage <- geneUsage(twb, HUMAN_TRBJ)
vis.gene.usage(twb.jusage, .main = 'Twins J-usage', .dodge = T)

## ----eval=TRUE, echo=TRUE, message=FALSE, fig.align='center', fig.height=6, fig.width=9----
vis.gene.usage(twb, HUMAN_TRBJ, .main = 'Twins J-usage', .dodge = F, .ncol = 2)

## ----eval=TRUE, echo=TRUE, fig.align='center', fig.height=4.5, fig.width=6----
twb.pca <- pca.segments(twb, .do.plot = F) 
vis.pca(pca.segments(twb, .do.plot = F, .genes = HUMAN_TRBV), .groups = list(GroupA = c(1,2), GroupB = c(3,4)))

## ----eval=TRUE, echo=TRUE, fig.align='center', fig.width=6, fig.height=5.5, warning=FALSE, message=FALSE----
km <- get.kmers(twb[[1]]$CDR3.amino.acid.sequence, .head = 100, .k = 7, .verbose = F)
d <- kmer.profile(km)
vis.logo(d)

## ----eval=TRUE, echo=TRUE------------------------------------------------
# data(twb)
twb.shared <- shared.repertoire(twb, .head = 1000, .verbose = F)
G <- mutation.network(twb.shared)
G

## ----eval=TRUE, echo=TRUE------------------------------------------------
# data(twb)
# twb.shared <- shared.repertoire(twb, .head = 1000)
# G <- mutation.network(twb.shared)
G <- set.group.vector(G, "twins", list(A = c(1,2), B = c(3,4)))  # <= refactor this
get.group.names(G, "twins", 1)
get.group.names(G, "twins", 300)
get.group.names(G, "twins", c(1,2,3), F)
get.group.names(G, "twins", 300, F)

# Because we have only two groups, we can assign more readable attribute.
V(G)$twin.names <- get.group.names(G, "twins")
V(G)$twin.names[1]
V(G)$twin.names[300]

## ----eval=TRUE, echo=TRUE------------------------------------------------
# data(twb)
# twb.shared <- shared.repertoire(twb, .head = 1000)
# G <- mutation.network(twb.shared)
head(mutated.neighbours(G, 1)[[1]])

## ----eval=TRUE, echo=TRUE------------------------------------------------
head(get.kmers(twb[[1]]$CDR3.amino.acid.sequence, 100, .meat = F, .verbose = F))
head(get.kmers(twb[[1]], .meat = T, .verbose = F))

## ----eval=TRUE, echo=TRUE------------------------------------------------
revcomp(c('AAATTT', 'ACGTTTGGA'))
cbind(bunch.translate(twb[[1]]$CDR3.nucleotide.sequence[1:10]),
      twb[[1]]$CDR3.amino.acid.sequence[1:10])
gc.content(twb[[1]]$CDR3.nucleotide.sequence[1:10])

## ----eval=TRUE, echo=TRUE------------------------------------------------
codon.variants('LQ')
translated.nucl.sequences(c('LQ', 'CASSLQ'))
reverse.translation('LQ')
translated.nucl.sequences('LQ', 'XXXXXG')
codon.variants('LQ', 'XXXXXG')
reverse.translation('LQ', 'XXXXXG')

