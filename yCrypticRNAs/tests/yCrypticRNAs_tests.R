#
# # Load yCrypticRNAs into your current workspace
#  library(yCrypticRNAs)
#
# # create a coverageDataSet containing the RNA-seq coverage values for all samples.
# # RNA-seq experiments were done in duplicates on wild-type cells (wt)
# # and cells in which Spt6, an histone chaperone was inactivated (mut).
# samples <- c("wt_rep1", "wt_rep2", "mut_rep1", "mut_rep2")
# samples_files <- system.file("extdata", paste0(samples, "_coverage.txt"), package = "yCrypticRNAs")
#
# sampleTable <- data.frame( samples, samples_files,
#                            c("wt", "wt", "mut", "mut"))
#
# data <- coverageDataSet(sampleTable = sampleTable)
#
#
# # load introns annotations
# data(introns)
#
# # Analysis on FLO8 gene
# # create geneCoverage data flor FLO8 gene
# flo8 <- gene_coverage(coverageDataSet = data, name = "YER109C", introns = introns)
# flo8
# plot(flo8)
#
# # calculate cryptic scores --------
#
# # using the ratio 3'/5' method
# ratio_score(geneCoverage = flo8)
#
# # using the 3' enrichemnt method
# enrichment_score(geneCoverage = flo8)
#
# # using the probabilistic method
# zscore_score(geneCoverage = flo8)
#
#
# # Calculating the cryptic scores for all the genes in the data
# genome_wide_scores(data, "ratio", "ratio.txt")
# genome_wide_scores(data, "enrichment", "enrichment.txt")
# genome_wide_scores(data, "probabilistic", "probabilistic.txt")
# unlink(c("ratio.txt", "enrichment.txt", "probabilistic.txt"))
#
# ## cryptic transcription start sites.
#
# # Annotation set
# #
# # data("saccer3_annotation")
# # annotations <- as.annotationsSet(saccer3_annotation)
# # fragments_files <- system.file("extdata", paste0(samples, "_fragments.txt"), package = "yCrypticRNAs")
# #
# # #all the methods
# # flo8_cTSS <- initiation_sites(name = "YER109C", fragments_files = fragments_files, annotations = annotations, introns = introns)
# # par(mfrow = c(3,2))
# # plot(flo8, cTSS = flo8_cTSS, method = "methodC_gaussian")
# # plot(flo8, cTSS = flo8_cTSS, method = "methodA")
# # plot(flo8, cTSS = flo8_cTSS, method = "methodB")
# # plot(flo8, cTSS = flo8_cTSS, method = "methodC")
# # plot(flo8, cTSS = flo8_cTSS, method = "methodD")
# #
# #
#
#
