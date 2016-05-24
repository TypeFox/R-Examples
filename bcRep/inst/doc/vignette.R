## ---- eval=FALSE---------------------------------------------------------
#  # Example: combine folders IMGT1a, IMGT1b and IMGT1c to a new datset NewProject
#  combineIMGT(folders = c("pathTo/IMGT1a", "pathTo/IMGT1b", "pathTo/IMGT1c"),
#              name = "NewProject")

## ---- echo=F-------------------------------------------------------------
library(bcRep)
library(pander)

## ---- echo=T-------------------------------------------------------------
# Example: 
data(aaseqtab)
aadistr<-aaDistribution(sequences = aaseqtab$CDR3_IMGT, numberSeq = TRUE)

# First 4 columns of Amino acid distribution table: 
        # for a sequence length of 13 AA (* = stop codon)
        # (aadistr$Amino_acid_distribution$`sequence length = 13`)


## ---- echo=F, results='asis'---------------------------------------------
library(pander)
pandoc.table(aadistr$Amino_acid_distribution$`sequence length = 13`[, 1:4])


## ---- echo=T, fig.width = 8, fig.height = 6, fig.align='center'----------
# Plot example for sequence lengths of 14-17 amino acids:
aadistr.part<-list(aadistr$Amino_acid_distribution[13:16], 
                   data.frame(aadistr$Number_of_sequences_per_length[13:16,]))
names(aadistr.part)<-names(aadistr)
plotAADistribution(aaDistribution.tab=aadistr.part, plotAADistr=TRUE, 
                   plotSeqN=TRUE, PDF=NULL) 

## ---- echo=F-------------------------------------------------------------
library(bcRep)
library(pander)

## ---- echo=T-------------------------------------------------------------
# Example:
data(aaseqtab)
trueDiv<-trueDiversity(sequences = aaseqtab$CDR3_IMGT, order = 1) 
        # using exponent of Shannon entropy

# True diversity of order 1 for amino acid length of 5 AA 
        # (trueDiv$True_diversity$'sequence length = 5')


## ---- echo=F, results='asis'---------------------------------------------
pandoc.table(trueDiv$True_diversity$'sequence length = 5')


## ---- eval=T, fig.width = 5, fig.height = 5, fig.align='center'----------
# True diversity for sequences of amino acid length 14-17: 
trueDiv.part<-list(trueDiv$True_diversity_order, trueDiv$True_diversity[13:16])
names(trueDiv.part)<-names(trueDiv)
plotTrueDiversity(trueDiversity.tab=trueDiv.part, mean.plot = F,color="darkblue", PDF=NULL)

## ---- eval=T, fig.width = 7, fig.height = 5, fig.align='center'----------
plotTrueDiversity(trueDiversity.tab=trueDiv, mean.plot = T,color="darkblue", PDF=NULL)

## ---- echo=F-------------------------------------------------------------
library(bcRep)

## ---- echo=T, collapse=TRUE, comment='#'---------------------------------
# Example:
data(clones.ind)
clones.giniIndex(clone.size=clones.ind$total_number_of_sequences, PDF = NULL)

## ---- fig.width = 4, fig.height = 5.3, fig.align='default'---------------
# Example:
data(summarytab)
Vgu<-geneUsage(genes = summarytab$V_GENE_and_allele, level = "subgroup", 
                functionality = summarytab$Functionality)
plotGeneUsage(geneUsage.tab = Vgu, plotFunctionality = TRUE, PDF = NULL, 
              title = "IGHV usage")

## ---- echo=F-------------------------------------------------------------
library(bcRep)

## ---- fig.width = 10, fig.height = 7, fig.align='center'-----------------
# Example: 
data(summarytab)
VDcomb.tab<-sequences.geneComb(family1 = summarytab$V_GENE_and_allele, 
                family2 = summarytab$D_GENE_and_allele, 
                level = "subgroup", abundance = "relative")
plotGeneComb(geneComb.tab = VDcomb.tab, withNA = FALSE, PDF = NULL)

## ---- echo=FALSE---------------------------------------------------------
library(bcRep)

## ---- collapse=TRUE------------------------------------------------------
# Example: filter for productive sequences
data(summarytab)
ProductiveSequences<-sequences.getProductives(summarytab)
# dimension of the summary table [rows, columns]:
dim(summarytab) 
# dimension of the summary table, filtered for productive sequences [rows, columns]:
dim(ProductiveSequences) 

## ---- echo=F-------------------------------------------------------------
library(bcRep)
library(pander)

## ---- eval=F, collapse = TRUE--------------------------------------------
#  # Example:
#  data(summarytab)
#  sequences.functionality(data = summarytab$Functionality, relativeValues=TRUE)
#  
#  # Proportion of productive and unproductive sequences:

## ---- echo=F,results='asis'----------------------------------------------
data(summarytab)
pandoc.table(sequences.functionality(data = summarytab$Functionality, relativeValues=TRUE))


## ---- echo=F-------------------------------------------------------------
# Example:
library(bcRep)
library(pander)


## ---- echo=T,results='asis'----------------------------------------------
data(mutationtab)
V.mutation<-sequences.mutation(mutationtab = mutationtab, sequence = "V", 
                               junctionFr = TRUE, rsRatio=TRUE)

# V.mutation$Number_of_mutations, first 6 lines:


## ---- echo=F, results='asis'---------------------------------------------
pandoc.table(head(V.mutation$Number_of_mutations))
# V.mutation$Number_of_mutations, first 6 lines:


## ------------------------------------------------------------------------
# V.mutation$Junction frame:


## ---- echo=F, results='asis'---------------------------------------------
pandoc.table(V.mutation$Junction_frame)

## ---- echo=F-------------------------------------------------------------
library(bcRep)
library(pander)


## ---- echo=T,results='asis', fig.align='center', fig.height=5, fig.width=9----
# Example:

data(mutationtab)
V.AAmut<-sequences.mutation.AA(mutationtab = mutationtab, sequence = "V")
plotSequencesMutationAA(mutationAAtab = V.AAmut, showChange = "hydropathy")


## ---- echo=F-------------------------------------------------------------
library(bcRep)
library(pander)


## ---- echo=T,results='asis', fig.align='center', fig.height=7, fig.width=9----
# Example:

data(mutationtab)
data(summarytab)
V.BaseMut<-sequences.mutation.base(mutationtab = mutationtab, summarytab = summarytab, sequence = "V", nrCores=1)
plotSequencesMutationBase(mutationBaseTab = V.BaseMut)


## ---- eval=FALSE---------------------------------------------------------
#  # Example:
#  data(aaseqtab)
#  data(summarytab)
#  clones.tab<-clones(aaseqtab = aaseqtab, summarytab = summarytab, ntseqtab = NULL,
#                     identity = 0.85, useJ = TRUE, dispD = FALSE, dispSeqID = FALSE,
#                     dispCDR3aa = FALSE, dispCDR3nt = FALSE,
#                     dispJunctionFr.ratio = FALSE, dispJunctionFr.list = FALSE,
#                     dispFunctionality.ratio = FALSE, dispFunctionality.list = FALSE,
#                     dispTotalSeq = FALSE, nrCores=1)
#  

## ---- eval=T, echo=T-----------------------------------------------------
# Example of a clonotype file (not from the command above)
data(clones.ind)
str(clones.ind, strict.width="cut", width=85)

## ---- echo=FALSE---------------------------------------------------------
library(bcRep)
library(pander)

## ---- collapse=TRUE------------------------------------------------------
# Example 1: Getting the 3 smallest clones (clones with 85% CDR3 identity)
clones.filtered1<-clones.filterSize(clones.tab=clones.ind, 
                                    column="total_number_of_sequences", 
                                    number=3, method="lower.tail")

    # Output of clones.filtered1[,c("total_number_of_sequences",
            # "unique_CDR3_sequences_AA", "CDR3_length_AA", "V_gene")]:

## ---- echo=F, results='asis'---------------------------------------------
pandoc.table(data.frame(clones.filtered1[,c("total_number_of_sequences","unique_CDR3_sequences_AA","CDR3_length_AA","V_gene")], row.names=NULL))


## ---- echo=T, results='markup', collapse=TRUE----------------------------
# Example 2: Getting 10% biggest and smallest clones 
# (clones with 85% CDR3 identity; column 4 = "total_number_of_sequences")
clones.filtered2<-clones.filterSize(clones.tab=clones.ind, column=4, propOfClones=0.1,
                                    method="two.tailed")
names(clones.filtered2) # a list with biggest and smallest 10% clones
dim(clones.filtered2$upper.tail) # dimension of table with biggest clones [rows, columns]
dim(clones.filtered2$lower.tail) # dimension of table with smallest clones [rows, columns]

# Example 3: Getting clones, that include 0.5% of all sequences 
# (clones with 85% CDR3 identity; column 4 = "total_number_of_sequences"")
clones.filtered3<-clones.filterSize(clones.tab=clones.ind, column=4,
                                    propOfSequences=0.005, method="two.tailed")
names(clones.filtered3) # a list with biggest and smallest 10% clones
dim(clones.ind) ## dimension of clones.ind table [rows, columns]
## >> number of rows is equal to number of rows of biggest and smallest clones
dim(clones.filtered3$upper.tail) # dimension of table with biggest clones [rows, columns]
dim(clones.filtered3$lower.tail) # dimension of table with smallest clones [rows, columns]

## ---- echo=FALSE---------------------------------------------------------
library(bcRep)

## ---- collapse=TRUE------------------------------------------------------
# Example
data(clones.ind)
productiveClones<-clones.filterFunctionality(clones.tab = clones.ind, 
                                             filter = "productive")
# dimension of clones.ind [rows, columns]:
dim(clones.ind) 
# dimension of clones.ind, filtered for productive sequences [rows, columns]:
dim(productiveClones) 

## ---- echo=FALSE---------------------------------------------------------
library(bcRep)

## ---- collapse=TRUE------------------------------------------------------
# Example
data(clones.ind)
inFrameClones<-clones.filterJunctionFrame(clones.tab = clones.ind, filter = "in-frame")
# dimension of clones.ind [rows, columns]:
dim(clones.ind) 
# dimension of clones.ind, filtered for in-frame sequences [rows, columns]:
dim(inFrameClones) 

## ---- echo=FALSE---------------------------------------------------------
library(bcRep)
library(pander)

## ---- results='markup'---------------------------------------------------
# Example:
data(clones.ind)
CDR3length<-clones.CDR3Length(CDR3Length = clones.ind$CDR3_length_AA, 
                  functionality = clones.ind$Functionality_all_sequences)

# output of some CDR3 length proportions (CDR3length$CDR3_length[1:3]): 

## ---- results='asis', echo=FALSE-----------------------------------------
pandoc.table(CDR3length$CDR3_length[1:3])

## ------------------------------------------------------------------------
# output of functionality ratios for some CDR3 lengths 
    # (CDR3length$CDR3_length_vs_functionality[,1:3)])

## ---- results='asis', echo=FALSE-----------------------------------------
pandoc.table(CDR3length$CDR3_length_vs_functionality[,1:3])
             

## ---- fig.align='center', fig.height=5, fig.width=9----------------------
plotClonesCDR3Length(CDR3Length = clones.ind$CDR3_length_AA, 
                     functionality = clones.ind$Functionality_all_sequences, 
                     title = "CDR3 length distribution", PDF = NULL)

## ---- fig.align='center', fig.width=3.5, fig.height=4.5------------------
# Example: Clone copy number distribution with and without outliers
plotClonesCopyNumber(copyNumber = clones.ind$total_number_of_sequences, 
                     withOutliers=TRUE, color = "darkblue", 
                     title = "Copy number distribution", PDF = NULL)
plotClonesCopyNumber(copyNumber = clones.ind$total_number_of_sequences, 
                     withOutliers=FALSE, color = "darkblue", 
                     title = "Copy number distribution", PDF = NULL)

## ---- echo=FALSE---------------------------------------------------------
library(bcRep)
library(pander)

## ---- collapse=TRUE------------------------------------------------------
data(aaseqtab)
data(aaseqtab2)
V.comp<-compare.geneUsage(gene.list = list(aaseqtab$V_GENE_and_allele, 
                                           aaseqtab2$V_GENE_and_allele), 
                               level = "subgroup", abundance = "relative", 
                          names = c("Individual1", "Individual2"), 
                          nrCores = 1)

## ---- echo=FALSE, results='asis'-----------------------------------------

pandoc.table(V.comp)

## ---- fig.align='center', fig.width=9, fig.height=6----------------------
plotCompareGeneUsage(comp.tab = V.comp, color = c("gray97", "darkblue"), PDF = NULL)

## ---- echo=FALSE---------------------------------------------------------
library(bcRep)

## ---- fig.align='center', fig.width=9, fig.height=6----------------------
data(aaseqtab)
data(aaseqtab2)
AAdistr.comp<-compare.aaDistribution(sequence.list = list(aaseqtab$CDR3_IMGT,
                                                          aaseqtab2$CDR3_IMGT), 
                                     names = c("Individual1", "Individual2"), 
                                     numberSeq = TRUE, nrCores = 1)

# Comparison of sequence length of 14-16 amino acids: 
## Individual 1: grep("14|15|16",names(AAdistr.comp$Amino_acid_distribution$Individual1)) 
        ## >> 13:15
## Individual 2: grep("14|15|16",names(AAdistr.comp$Amino_acid_distribution$Individual2)) 
        ## >> 12:14
AAdistr.comp.part<-list(list(AAdistr.comp$Amino_acid_distribution$Individual1[13:15], 
                             AAdistr.comp$Amino_acid_distribution$Individual2[12:14]),
                        list(AAdistr.comp$Number_of_sequences_per_length$Individual1[13:15],
                             AAdistr.comp$Number_of_sequences_per_length$Individual2[12:14]))
names(AAdistr.comp.part)<-names(AAdistr.comp)
names(AAdistr.comp.part$Amino_acid_distribution)<-
  names(AAdistr.comp$Amino_acid_distribution)
names(AAdistr.comp.part$Number_of_sequences_per_length)<-
  names(AAdistr.comp$Number_of_sequences_per_length)

plotCompareAADistribution(comp.tab = AAdistr.comp.part, plotSeqN = TRUE, 
                          colors=c("darkblue","darkred"), PDF = NULL)

## ---- echo=FALSE---------------------------------------------------------
library(bcRep)

## ---- fig.align='center', fig.width=9, fig.height=7----------------------
data(aaseqtab)
data(aaseqtab2)
trueDiv.comp<-compare.trueDiversity(sequence.list = list(aaseqtab$CDR3_IMGT, 
                                                         aaseqtab2$CDR3_IMGT), 
                                    names = c("Individual1", "Individual2"), 
                                    order = 1, nrCores = 1)

# Comparison of sequence length of 14-16 amino acids: 
grepindex1<-grep("14|15|16",names(trueDiv.comp$Individual1))
grepindex2<-grep("14|15|16",names(trueDiv.comp$Individual2))
trueDiv.comp.part<-list(trueDiv.comp$True_diversity_order, 
                        trueDiv.comp$Individual1[grepindex1],
                        trueDiv.comp$Individual2[grepindex2])
names(trueDiv.comp.part)<-names(trueDiv.comp)
plotCompareTrueDiversity(comp.tab = trueDiv.comp.part, mean.plot = F, colors=c("darkblue","darkred"), 
                         PDF = NULL)
plotCompareTrueDiversity(comp.tab = trueDiv.comp, mean.plot = T, colors = c("darkblue","darkred"), 
                         PDF = NULL)

## ---- echo=FALSE---------------------------------------------------------
library(bcRep)
library(pander)

## ---- collapse=TRUE------------------------------------------------------
# Example:
data(clones.allind) # includes 2948 clones for 7 individuals
colnames(clones.allind) # colnames of clones.allind

## ---- eval=FALSE---------------------------------------------------------
#  # Example:
#  data(clones.allind)
#  sharedclones<-clones.shared(clones.tab = clones.allind)
#                # with default parameters: identity = 0.85, useJ = TRUE, ....

## ---- eval=FALSE---------------------------------------------------------
#  # Example:
#  data(clones.allind)
#  sharedclones<-clones.shared(clones.tab = clones.allind, identity = 0.85, useJ = TRUE)
#  sharedclones.summary<-clones.shared.summary(shared.tab = sharedclones,
#                                              clones.tab = clones.allind)

## ---- eval=TRUE, echo=TRUE-----------------------------------------------
# Example:
data(clones.allind)
sharedclones<-clones.shared(clones.tab = clones.allind, identity = 0.85, useJ = TRUE)
sharedclones.summary<-clones.shared.summary(shared.tab = sharedclones, 
                                            clones.tab = clones.allind)

# sharedclones.summary:

## ---- echo=FALSE, results='asis'-----------------------------------------
pandoc.table(sharedclones.summary)

## ---- echo=T, eval=F-----------------------------------------------------
#  # Example:
#  data(vgenes)
#  vgenes[1:5,1:5]

## ---- echo=F, results='asis', warning=FALSE------------------------------
data(vgenes) 
pandoc.table(vgenes[1:5,1:5])

## ---- echo=T, eval=F-----------------------------------------------------
#  # Example for 5 samples and 5 genes:
#  data(vgenes)
#  geneUsage.distance(geneUsage.tab = vgenes[1:5,1:5], method = "bc")

## ---- echo=F, results='asis', warning=FALSE------------------------------
data(vgenes) 
pandoc.table(geneUsage.distance(geneUsage.tab = vgenes[1:5,1:5], method = "bc"))

## ---- echo=T, eval=F-----------------------------------------------------
#  # Example 1:
#  ## Divide sequences into subsets of same length and then apply Levenshtein distance:
#  data(clones.ind)
#  clones.ind<-subset(clones.ind, CDR3_length_AA==7) #take only a subset of clones.ind (CDR3 sequences of 7 AA)
#  dist1<-sequences.distance(sequences = clones.ind$unique_CDR3_sequences_AA,
#       method = "levenshtein", divLength=TRUE)
#  # Subset of output for 7AA sequences:
#  dist1[[1]][1:5,1:5]

## ---- echo=F, results='asis', warning=FALSE------------------------------
data(clones.ind)
clones.ind<-subset(clones.ind, CDR3_length_AA==7)
dist1<-sequences.distance(sequences = clones.ind$unique_CDR3_sequences_AA, 
     method = "levenshtein", divLength=TRUE)
pandoc.table(dist1[[1]][1:5,1:5])


## ---- echo=T, eval=F-----------------------------------------------------
#  # Example 2:
#  ## Use all sequences with a length of 7 or 8 AA for cosine distance matrix:
#  data(clones.allind)
#  clones.allind<-subset(clones.allind, CDR3_length_AA %in% c(9,31)) #take only a subset of clones.ind (CDR3 sequences of 9 or 31 AA)
#  dist2<-sequences.distance(sequences = clones.allind$unique_CDR3_sequences_AA,
#       groups = clones.allind$individuals, method = "cosine", divLength=FALSE)
#  # Subset of output:
#  dist2[[1]][1:4,1:4]

## ---- echo=F, results='asis', warning=FALSE------------------------------
data(clones.allind)
clones.allind<-subset(clones.allind, CDR3_length_AA %in% c(9,31))

dist2<-sequences.distance(sequences = clones.allind$unique_CDR3_sequences_AA, 
     groups = clones.allind$individuals, method = "cosine", divLength=FALSE)
pandoc.table(dist2[[1]][1:4,1:4])

## ---- echo=T, eval=T, results='asis', warning=FALSE----------------------
# Example for gene usage data:
data(vgenes) 
vgenes<-vgenes[,1:10] # include only genes 1-10
gu.dist<-geneUsage.distance(geneUsage.tab = vgenes, method = "bc")
gu.pcoa<-dist.PCoA(dist.tab = gu.dist, correction = "cailliez")

## ---- results='markup'---------------------------------------------------
str(gu.pcoa) # a PCoA object

## ---- results='asis', fig.align='center',  fig.height=4.5, fig.width=7----
# Example of a 'groups' data.frame: in the case, there are two groups:
groups.vec<-cbind(rownames(vgenes),c(rep("group1",5),rep("group2",5)))
colnames(groups.vec)<-c("Sample","Group")
pandoc.table(groups.vec)

plotDistPCoA(pcoa.tab = gu.pcoa, groups= groups.vec, axes = c(1,2), 
             names = rownames(vgenes), plotLegend = TRUE)

## ---- echo=T, eval=T, fig.align='center', fig.height=6.5, fig.width=7, results='asis', warning=FALSE----
# Example for sequence data, divided into subsets of the same length:
data(clones.ind)
clones.ind<-subset(clones.ind, CDR3_length_AA %in% c(7:10)) # include only CDR3 sequences of 7-10 AA
seq.dist<-sequences.distance(sequences = clones.ind$unique_CDR3_sequences_AA, 
     method = "levenshtein", divLength=T)
seq.pcoa<-dist.PCoA(dist.tab = seq.dist, correction = "none")

# Example of a 'groups' data.frame: 
## in the case, there are no groups (all samples have the same group):
groups.vec<-unlist(apply(data.frame(clones.ind$unique_CDR3_sequences_AA),1,
                         function(x){strsplit(x,split=", ")[[1]]}))
groups.vec<-cbind(groups.vec, 1)

plotDistPCoA(pcoa.tab = seq.pcoa, groups=groups.vec, axes = c(1,2), 
             plotCorrection = FALSE, 
             title = 'PCoA for sequences of length 7-10 AA', plotLegend=FALSE)    

