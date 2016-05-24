## Run test file, to be sure that all functions works
# res <- gsub("[1] ", "", capture.output(source("snpEnrichment_testFile.R")), fixed = TRUE);
# table(res); # should return TRUE 113
# grep("FALSE", res) # should return integer(0)

rm(list = ls())

testPackage <- function () {
    test <- function(f, error = FALSE) {
        if (error==TRUE) {
            return(class(try(f, silent = TRUE))=="try-error")
        } else {
            return(class(try(f, silent = TRUE))!="try-error")
        }
    }
    res <- c(
        test(library("snpEnrichment"), error = FALSE), #1
        # library(parallel)
        # library(snpStats)

        test(chrA <- chromosome(), error = FALSE), #2
        test(chrC <- chromosome(), error = FALSE), #3
        test(is.chromosome(list()), error = FALSE), #4
        test(is.chromosome(1), error = FALSE), #5
        test(is.chromosome(chrA), error = FALSE), #6
        test(is.chromosome(c(chrA, chrC)), error = FALSE), #7
        test(is.chromosome(list(chrA, b = "char")), error = FALSE), #8
        test(is.chromosome(c(chrA, b = list(12, chrC))), error = FALSE), #9
        test(print(chrA, type = "eSNP"), error = FALSE), #10
        # print(10)
        test(enrichA <- enrichment(), error = FALSE), #11
        test(enrichC <- enrichment(), error = FALSE), #12
        test(is.enrichment(list()), error = FALSE), #13
        test(is.enrichment(1), error = FALSE), #14
        test(is.enrichment(enrichA), error = FALSE), #15
        test(is.enrichment(c(enrichA, enrichC)), error = FALSE), #16
        test(is.enrichment(list(enrichA, b = "char")), error = FALSE), #17
        test(is.enrichment(c(enrichA, b = list(12, enrichC))), error = FALSE), #18
        test(print(enrichA, what = "All"), error = FALSE), #19
        test(toyEnrich <- enrichment(), error = FALSE), #20
        # print(20)
        test(trash <- capture.output(show(toyEnrich)), error = FALSE), #21
        test(data(toyEnrichment), error = FALSE), #22
        test(toyEnrich["Loss"] <- toyEnrichment["Loss"], error = FALSE), #23
        test(toyEnrich["Loss"], error = FALSE), #24
        test(toyEnrich <- enrichment(Loss = toyEnrichment["Loss"], eSNP = toyEnrichment["eSNP"]), error = FALSE), #25
        test(toyEnrich <- enrichment(Loss = toyEnrichment["Loss"]), error = FALSE), #26
        test(print(toyEnrich, what = "All"), error = FALSE), #27
        test(data(toyEnrichment), error = FALSE), #28
        test(trash <- capture.output(reSample(object = toyEnrichment, nSample = 10, empiricPvalue = FALSE, MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5), mc.cores = 1, onlyGenome = FALSE)), error = FALSE), #29
        test(print(toyEnrichment, what = "All"), error = FALSE), #30
        # print(30)
        test(print(toyEnrichment, what = "Genome"), error = FALSE), #31
        test(print(toyEnrichment, what = 1), error = FALSE), #32
        test(print(toyEnrichment, what = c(2, 4)), error = FALSE), #33
        test(print(toyEnrichment, what = seq(22)), error = FALSE), #34
        test(plot(toyEnrichment, what = "Genome", type = c("eSNP", "xSNP")), error = FALSE), #35
        test(plot(toyEnrichment, what = "Genome", type = "xSNP"), error = FALSE), #36
        test(plot(toyEnrichment, what = 2, type = c("eSNP", "xSNP")), error = FALSE), #37
        test(plot(toyEnrichment, what = 22, type = "eSNP"), error = FALSE), #38
        test(plot(toyEnrichment, what = 2, type = c("eSNP", "xSNP"), ggplot = TRUE), error = FALSE), #39
        test(plot(toyEnrichment, what = 22, type = "eSNP", ggplot = TRUE), error = FALSE), #40
        # print(40)
        test(data(toyEnrichment) , error = FALSE), #41
        test(trash <- capture.output(reSample(object = toyEnrichment, nSample = 10, empiricPvalue = TRUE, MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5), mc.cores = 1, onlyGenome = TRUE)), error = FALSE), #42
        test(print(toyEnrichment, what = "All"), error = FALSE), #43
        test(print(toyEnrichment, what = "Genome"), error = FALSE), #44
        test(print(toyEnrichment, what = 1), error = FALSE), #45
        test(print(toyEnrichment, what = c(2, 4)), error = FALSE), #46
        test(print(toyEnrichment, what = seq(22)), error = FALSE), #47
        test(plot(toyEnrichment, what = "Genome", type = c("eSNP", "xSNP")), error = FALSE), #48
        test(plot(toyEnrichment, what = "Genome", type = "xSNP"), error = FALSE), #49
        test(plot(toyEnrichment, what = 2, type = c("eSNP", "xSNP")), error = TRUE), #50
        # print(50)
        test(plot(toyEnrichment, what = 22, type = "eSNP"), error = TRUE), #51
        test(a <- capture.output(dev.off()), error = FALSE), #52
        test(snpInfoDir <- system.file("extdata/snpInfo", package = "snpEnrichment"), error = FALSE), #53
        test(snpListDir <- system.file("extdata/List", package = "snpEnrichment"), error = FALSE), #54
        test(signalFile <- system.file("extdata/Signal/toySignal.txt", package = "snpEnrichment"), error = FALSE), #55
        test(data(transcript), error = FALSE), #56
        test(trash <- capture.output(initFiles(pattern = "Chrom", snpInfoDir, signalFile, mc.cores = 1)), error = FALSE), #57
        test(trash <- capture.output(toyEnrichment <- readEnrichment(pattern = "Chrom", signalFile, transcriptFile = FALSE, snpListDir, snpInfoDir, distThresh = 500, sigThresh = 0.05, LD = FALSE, ldDir = NULL, mc.cores = 1)), error = FALSE), #58
        test(print(toyEnrichment, what = "All"), error = FALSE), #59
        test(trash <- capture.output(initFiles(pattern = "Chrom", snpInfoDir, signalFile, mc.cores = 1)), error = FALSE), #60
        # print(60)
        test(trash <- capture.output(toyEnrichment <- readEnrichment(pattern = "Chrom", signalFile, transcriptFile = transcript, snpListDir, snpInfoDir, distThresh = 500, sigThresh = 0.05, LD = FALSE, ldDir = NULL, mc.cores = 1)), error = FALSE), #61
        test(print(toyEnrichment, what = "All"), error = FALSE), #62
        test(trash <- capture.output(initFiles(pattern = "Chrom", snpInfoDir, signalFile, mc.cores = 1)), error = FALSE), #63
        test(trash <- capture.output(writeLD(pattern = "Chrom", snpInfoDir, signalFile, ldDir = NULL, ldThresh = 0.8, mc.cores = 1)), error = FALSE), #64
        test(trash <- capture.output(toyEnrichment <- readEnrichment(pattern = "Chrom", signalFile, transcriptFile = transcript, snpListDir, snpInfoDir, distThresh = 500, sigThresh = 0.05, LD = TRUE, ldDir = NULL, mc.cores = 1)), error = FALSE), #65
        test(print(toyEnrichment, what = "All"), error = FALSE), #66
        test(trash <- capture.output(initFiles(pattern = "Chrom", snpInfoDir, signalFile, mc.cores = 1)), error = FALSE), #67
        test(trash <- capture.output(writeLD(pattern = "Chrom", snpInfoDir, signalFile, ldDir = NULL, ldThresh = 0.8, mc.cores = 1)), error = FALSE), #68
        test(trash <- capture.output(toyEnrichment <- readEnrichment(pattern = "Chrom", signalFile, transcriptFile = transcript, snpListDir, snpInfoDir, distThresh = 500, sigThresh = 0.05, LD = TRUE, ldDir = NULL, mc.cores = 1)), error = FALSE), #69
        test(trash <- capture.output(reSample(object = toyEnrichment, nSample = 10, empiricPvalue = FALSE, MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5), mc.cores = 1, onlyGenome = FALSE)), error = FALSE), #70
        # print(70)
        test(print(toyEnrichment, what = "All"), error = FALSE), #71
        test(trash <- capture.output(initFiles(pattern = "Chrom", snpInfoDir, signalFile, mc.cores = 1)), error = FALSE), #72
        test(trash <- capture.output(writeLD(pattern = "Chrom", snpInfoDir, signalFile, ldDir = NULL, ldThresh = 0.8, mc.cores = 1)), error = FALSE), #73
        test(trash <- capture.output(toyEnrichmentM2 <- readEnrichment(pattern = "Chrom", signalFile, transcriptFile = transcript, snpListDir, snpInfoDir, distThresh = 500, sigThresh = 0.05, LD = TRUE, ldDir = NULL, mc.cores = 1)), error = FALSE), #74
        test(trash <- capture.output(reSample(object = toyEnrichmentM2, nSample = 10, empiricPvalue = FALSE, MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5), mc.cores = 1, onlyGenome = FALSE)), error = FALSE), #75
        test(print(toyEnrichmentM2, what = "All"), error = FALSE), #76
        test(data(toyEnrichment), error = FALSE), #77
        test(excludeFile <- system.file("extdata/Exclude/toyExclude.txt", package = "snpEnrichment"), error = FALSE), #78
        test(trash <- capture.output(toyM1_exclude <- excludeSNP(toyEnrichment, excludeFile, mc.cores = 1)), error = FALSE), #79
        test(print(toyM1_exclude, what = "All"), error = FALSE), #80
        # print(80)
        test(excludeFile <- system.file("extdata/Exclude/toyExclude.txt", package = "snpEnrichment"), error = FALSE), #81
        test(trash <- capture.output(toyM1_exclude <- excludeSNP(toyEnrichment, excludeFile, mc.cores = 1)), error = FALSE), #82
        test(print(toyM1_exclude, what = "All"), error = FALSE), #83
        test(excludeFile <- c(
            "rs4376885", "rs17690928", "rs6460708", "rs13061537", "rs11769827",
            "rs12717054", "rs2907627", "rs1380109", "rs7024214", "rs7711972",
            "rs9658282", "rs11750720", "rs1793268", "rs774568", "rs6921786",
            "rs1699031", "rs6994771", "rs16926670", "rs465612", "rs3012084",
            "rs354850", "rs12803455", "rs13384873", "rs4364668", "rs8181047",
            "rs2179993", "rs12049335", "rs6079926", "rs2175144", "rs11564427",
            "rs7786389", "rs7005565", "rs17423335", "rs12474102", "rs191314",
            "rs10513168", "rs1711437", "rs1992620", "rs283115", "rs10754563",
            "rs10851727", "rs2173191", "rs7661353", "rs1342113", "rs7042073",
            "rs1567445", "rs10120375", "rs550060", "rs3761218", "rs4512977"
        ), error = FALSE), #84
        test(trash <- capture.output(toyM1_exclude <- excludeSNP(toyEnrichment, excludeFile, mc.cores = 1)), error = FALSE), #85
        test(print(toyM1_exclude, what = "All"), error = FALSE), #86
        test(data(toyEnrichment) , error = FALSE), #87
        test(trash <- capture.output(reSample(object = toyEnrichment, nSample = 10, empiricPvalue = FALSE, MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5), mc.cores = 1, onlyGenome = FALSE)), error = FALSE), #88
        test(excludeFile <- system.file("extdata/Exclude/toyExclude.txt", package = "snpEnrichment"), error = FALSE), #89
        test(trash <- capture.output(toyM1_exclude <- excludeSNP(toyEnrichment, excludeFile, mc.cores = 1)), error = FALSE), #90
        # print(90)
        test(trash <- capture.output(reSample(object = toyM1_exclude, nSample = 10, empiricPvalue = FALSE, MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5), mc.cores = 1, onlyGenome = FALSE)), error = FALSE), #91
        test(trash <- capture.output(res_toyM1 <- compareEnrichment(object.x = toyEnrichment, object.y = toyM1_exclude, pattern = "Chrom", nSample = 10, empiricPvalue = FALSE, mc.cores = 1, onlyGenome = FALSE)), error = FALSE), #92
        test(data(toyEnrichment) , error = FALSE), #93
        test(trash <- capture.output(reSample(object = toyEnrichment, nSample = 10, empiricPvalue = FALSE, MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5), mc.cores = 1, onlyGenome = FALSE)), error = FALSE), #94
        test(excludeFile <- system.file("extdata/Exclude/toyExclude.txt", package = "snpEnrichment"), error = FALSE), #95
        test(trash <- capture.output(toyM1_exclude <- excludeSNP(toyEnrichment, excludeFile, mc.cores = 1)), error = FALSE), #96
        test(trash <- capture.output(reSample(object = toyM1_exclude, nSample = 10, empiricPvalue = FALSE, MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5), mc.cores = 1, onlyGenome = FALSE)), error = FALSE), #97
        test(trash <- capture.output(reSample(object = toyM1_exclude, nSample = 10, empiricPvalue = TRUE, MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5), mc.cores = 1, onlyGenome = FALSE)), error = FALSE), #98
        test(trash <- capture.output(res_toyM1 <- compareEnrichment(object.x = toyEnrichment, object.y = toyM1_exclude, pattern = "Chrom", nSample = 10, empiricPvalue = FALSE, mc.cores = 1, onlyGenome = FALSE)), error = FALSE), #99
        test(snpInfoDir <- system.file("extdata/snpInfo", package = "snpEnrichment"), error = FALSE), #100
        # print(100)
        test(snpListDir <- system.file("extdata/List", package = "snpEnrichment"), error = FALSE), #101
        test(signalFile <- system.file("extdata/Signal/toySignal.txt", package = "snpEnrichment"), error = FALSE), #102
        test(data(transcript), error = FALSE), #103
        test(trash <- capture.output(initFiles(pattern = "Chrom", snpInfoDir, signalFile, mc.cores = 2)), error = FALSE), #104
        test(trash <- capture.output(writeLD(pattern = "Chrom", snpInfoDir, signalFile, ldDir = NULL, ldThresh = 0.8, mc.cores = 2)), error = FALSE), #105
        test(trash <- capture.output(toyEnrichment <- readEnrichment(pattern = "Chrom", signalFile, transcriptFile = transcript, snpListDir, snpInfoDir, distThresh = 500, sigThresh = 0.05, LD = TRUE, ldDir = NULL, mc.cores = 2)), error = FALSE), #106
        test(trash <- capture.output(reSample(object = toyEnrichment, nSample = 10, empiricPvalue = FALSE, MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5), mc.cores = 2, onlyGenome = FALSE)), error = FALSE), #107
        test(excludeFile <- system.file("extdata/Exclude/toyExclude.txt", package = "snpEnrichment"), error = FALSE), #108
        test(trash <- capture.output(toyM1_exclude <- excludeSNP(toyEnrichment, excludeFile, mc.cores = 2)), error = FALSE), #109
        test(trash <- capture.output(res_toyM1 <- compareEnrichment(object.x = toyEnrichment, object.y = toyM1_exclude, pattern = "Chrom", nSample = 10, empiricPvalue = FALSE, mc.cores = 2, onlyGenome = FALSE)), error = FALSE), #110
        # print(110)
        test(data(toyEnrichment), error = FALSE), #111
        test(trash <- getEnrichSNP(toyEnrichment, type = "eSNP"), error = FALSE), #112
        test(trash <- getEnrichSNP(toyEnrichment, type = "xSNP"), error = FALSE) #113
    )
    if (any(res%in%FALSE)) {
        errMessage <- paste("The following expressions failed test:\n", paste(which(res%in%FALSE), collapse = ", "))
        stop(errMessage)
    } else {}
    return(invisible())
}