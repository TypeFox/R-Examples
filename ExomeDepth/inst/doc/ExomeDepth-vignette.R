## ----exons---------------------------------------------------------------
library(ExomeDepth)
data(exons.hg19)
print(head(exons.hg19))

## ----read.count, eval=FALSE----------------------------------------------
#  data(exons.hg19)
#  my.counts <- getBamCounts(bed.frame = exons.hg19,
#                            bam.files = my.bam,
#                            include.chr = FALSE,
#                            referenceFasta = fasta)

## ----Rsamtools.load------------------------------------------------------
library(ExomeDepth)
data(ExomeCount)
ExomeCount.dafr <- as(ExomeCount[,  colnames(ExomeCount)], 'data.frame')
ExomeCount.dafr$chromosome <- gsub(as.character(ExomeCount.dafr$space), 
                                        pattern = 'chr', 
                                        replacement = '')  ##remove the annoying chr letters
print(head(ExomeCount.dafr))

## ----<read.count---------------------------------------------------------
data(exons.hg19.X)
head(exons.hg19.X)

## ----first.test----------------------------------------------------------
test <- new('ExomeDepth',
            test = ExomeCount.dafr$Exome2,
            reference = ExomeCount.dafr$Exome3,
            formula = 'cbind(test, reference) ~ 1',
            subset.for.speed = seq(1, nrow(ExomeCount.dafr), 100))

show(test)

## ----reference.selection-------------------------------------------------
my.test <- ExomeCount$Exome4
my.ref.samples <- c('Exome1', 'Exome2', 'Exome3')
my.reference.set <- as.matrix(ExomeCount.dafr[, my.ref.samples])
my.choice <- select.reference.set (test.counts = my.test, 
                                   reference.counts = my.reference.set, 
                                   bin.length = (ExomeCount.dafr$end - ExomeCount.dafr$start)/1000, 
                                   n.bins.reduced = 10000)

print(my.choice[[1]])

## ----construct.ref-------------------------------------------------------
my.matrix <- as.matrix( ExomeCount.dafr[, my.choice$reference.choice, drop = FALSE])
my.reference.selected <- apply(X = my.matrix, 
                               MAR = 1, 
                               FUN = sum)

## ----build.complete------------------------------------------------------
all.exons <- new('ExomeDepth',
                 test = my.test,
                 reference = my.reference.selected,
                 formula = 'cbind(test, reference) ~ 1')

## ----call.CNVs-----------------------------------------------------------
all.exons <- CallCNVs(x = all.exons, 
                      transition.probability = 10^-4, 
                      chromosome = ExomeCount.dafr$space, 
                      start = ExomeCount.dafr$start, 
                      end = ExomeCount.dafr$end, 
                      name = ExomeCount.dafr$names)
head(all.exons@CNV.calls)

## ----write.results, eval = FALSE-----------------------------------------
#  output.file <- 'exome_calls.csv'
#  write.csv(file = output.file,
#            x = all.exons@CNV.calls,
#            row.names = FALSE)

## ----ranking-------------------------------------------------------------
head(all.exons@CNV.calls[ order ( all.exons@CNV.calls$BF, decreasing = TRUE),])

## ----Conrad--------------------------------------------------------------
data(Conrad.hg19)
head(Conrad.hg19.common.CNVs)

## ----check.chr-----------------------------------------------------------
levels(GenomicRanges::seqnames(Conrad.hg19.common.CNVs))

## ----anno.extra----------------------------------------------------------
all.exons <- AnnotateExtra(x = all.exons,
                           reference.annotation = Conrad.hg19.common.CNVs,
                           min.overlap = 0.5,
                           column.name = 'Conrad.hg19')

## ----anno.check----------------------------------------------------------
print(head(all.exons@CNV.calls))

## ----prep.GRanges--------------------------------------------------------
exons.hg19.GRanges <- GenomicRanges::GRanges(seqnames = exons.hg19$chromosome,
                                    IRanges::IRanges(start=exons.hg19$start,end=exons.hg19$end),
                                    names = exons.hg19$name)
all.exons <- AnnotateExtra(x = all.exons,
                           reference.annotation = exons.hg19.GRanges,
                           min.overlap = 0.0001,
                           column.name = 'exons.hg19')
all.exons@CNV.calls[3:6,]

## ----echo = TRUE, fig.width = 4, fig.height = 4--------------------------
plot (all.exons,
      sequence = '1',
      xlim = c(25598981 - 100000, 25633433 + 100000),
      count.threshold = 20,
      main = 'RHD gene',
      cex.lab = 0.8,
      with.gene = TRUE)

## ----loop, eval=FALSE----------------------------------------------------
#  
#  #### get the annotation datasets to be used later
#  data(Conrad.hg19)
#  exons.hg19.GRanges <- GRanges(seqnames = exons.hg19$chromosome,
#                                IRanges(start=exons.hg19$start,end=exons.hg19$end),
#                                names = exons.hg19$name)
#  
#  
#  ### prepare the main matrix of read count data
#  ExomeCount.mat <- as.matrix(ExomeCount.dafr[, grep(names(ExomeCount.dafr), pattern = 'Exome.*')])
#  nsamples <- ncol(ExomeCount.mat)
#  
#  ### start looping over each sample
#  for (i in 1:nsamples) {
#  
#  #### Create the aggregate reference set for this sample
#    my.choice <- select.reference.set (test.counts =  ExomeCount.mat[,i],
#                                       reference.counts = ExomeCount.mat[,-i],
#                                       bin.length = (ExomeCount.dafr$end - ExomeCount.dafr$start)/1000,
#                                       n.bins.reduced = 10000)
#  
#    my.reference.selected <- apply(X = ExomeCount.mat[, my.choice$reference.choice, drop = FALSE],
#                                   MAR = 1,
#                                   FUN = sum)
#  
#    message('Now creating the ExomeDepth object')
#    all.exons <- new('ExomeDepth',
#                     test = ExomeCount.mat[,i],
#                     reference = my.reference.selected,
#                     formula = 'cbind(test, reference) ~ 1')
#  
#  ################ Now call the CNVs
#    all.exons <- CallCNVs(x = all.exons,
#                          transition.probability = 10^-4,
#                          chromosome = ExomeCount.dafr$space,
#                          start = ExomeCount.dafr$start,
#                          end = ExomeCount.dafr$end,
#                          name = ExomeCount.dafr$names)
#  
#  ########################### Now annotate the ExomeDepth object
#    all.exons <- AnnotateExtra(x = all.exons,
#                               reference.annotation = Conrad.hg19.common.CNVs,
#                               min.overlap = 0.5,
#                               column.name = 'Conrad.hg19')
#  
#    all.exons <- AnnotateExtra(x = all.exons,
#                               reference.annotation = exons.hg19.GRanges,
#                               min.overlap = 0.0001,
#                               column.name = 'exons.hg19')
#  
#    output.file <- paste('Exome_', i, 'csv', sep = '')
#    write.csv(file = output.file, x = all.exons@CNV.calls, row.names = FALSE)
#  
#  }
#  
#  

## ----everted, eval = FALSE-----------------------------------------------
#  data(genes.hg19)
#  everted <- count.everted.reads (bed.frame = genes.hg19,
#                                 bam.files = bam.files.list,
#                                 min.mapq = 20,
#                                 include.chr = FALSE)
#  

## ----session-------------------------------------------------------------
sessionInfo()

