### R code from vignette source 'RPPanalyzer.rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: loadPackage (eval = FALSE)
###################################################
## library(RPPanalyzer)


###################################################
### code chunk number 2: loadSD (eval = FALSE)
###################################################
## ## define path to example files
## dataDir <- system.file("extdata",package="RPPanalyzer")
## ## change working directory
## setwd(dataDir)
## ## store example sample description in a variable
## sampledescription <- read.delim("sampledescription.txt")
## ## show sample description header
## head(sampledescription)


###################################################
### code chunk number 3: loadSlides (eval = FALSE)
###################################################
## ## change directory to example files
## dataDir <- system.file("extdata",package="RPPanalyzer")
## setwd(dataDir)
## ## store example slide description in a variable
## slidedescription <- read.delim("slidedescription.txt")
## ## show slide description header
## head(slidedescription)


###################################################
### code chunk number 4: dataPreprochelp (eval = FALSE)
###################################################
## ?dataPreproc


###################################################
### code chunk number 5: dataPreprocEx (eval = FALSE)
###################################################
## ## change directory to example files
## dataDir <- system.file("extdata",package="RPPanalyzer")
## setwd(dataDir)
## ## pre-process data
## preprocessedData <- dataPreproc(dataDir=dataDir,blocks=12,
##                                 spot="aushon",exportNo=2,correct="both")


###################################################
### code chunk number 6: changedirectory (eval = FALSE)
###################################################
## setwd(dataDir)


###################################################
### code chunk number 7: readdata (eval = FALSE)
###################################################
## dataDir <- system.file("extdata",package="RPPanalyzer")
## setwd(dataDir)
## rawdata <- read.Data(blocksperarray=12,spotter="aushon",
##                      printFlags=FALSE)


###################################################
### code chunk number 8: writeText (eval = FALSE)
###################################################
## write.Data(rawdata,FileNameExtension="test_data")


###################################################
### code chunk number 9: correctBackgroundDil (eval = FALSE)
###################################################
## ## import raw data
## fgRaw.tmp <- read.delim("test_dataexpression.txt",stringsAsFactors=FALSE,
##                         row.names=NULL,header=TRUE)
## fgRaw <- read.delim("test_dataexpression.txt",
##                     skip=max(which(fgRaw.tmp[,1]==""))+1,
##                     stringsAsFactors=FALSE,row.names=NULL,header=TRUE)
## ## remove NAs
## fgNAVec <- which(is.na(fgRaw[,"ID"]))
## if(length(fgNAVec) > 0){
##   fgRaw <- fgRaw[-fgNAVec,]
## }
## colnames(fgRaw) <- sub("X","",gsub("\\.","-",colnames(fgRaw)))
## ## correct data for BG noise
## correctedData <- correctDilinterc(
##   dilseries=fgRaw[which(fgRaw$sample_type=="control" & 
##                           !is.na(fgRaw$dilSeriesID)),],
##   arraydesc=rawdata$arraydescription,
##   timeseries=fgRaw[which(fgRaw$sample_type=="measurement"),],exportNo=2)
## ## correct negative values
## if(min(correctedData[,colnames(rawdata$arraydescription)]) < 0){
##   correctedData[,colnames(rawdata$arraydescription)] <- 
##     correctedData[,colnames(rawdata$arraydescription)] + 
##     abs(min(correctedData[,colnames(rawdata$arraydescription)]))+1
## }


###################################################
### code chunk number 10: loadHKdata (eval = FALSE)
###################################################
## ## load data set
## dataDir <- system.file("data",package="RPPanalyzer")
## setwd(dataDir)
## load("HKdata.RData")
## data(HKdata)


###################################################
### code chunk number 11: plotQC (eval = FALSE)
###################################################
## plotQC(HKdata,file="control_samples.pdf",arrays2rm=c("protein"))


###################################################
### code chunk number 12: plotmeasurements (eval = FALSE)
###################################################
## plotMeasurementsQC(HKdata,file="control_measurements.pdf",
##                    arrays2rm=c("protein"))


###################################################
### code chunk number 13: plotqq (eval = FALSE)
###################################################
## plotqq(HKdata,fileName="qqplot_measurements.pdf")


###################################################
### code chunk number 14: correctBackground (eval = FALSE)
###################################################
## dataBGcorrected <- correctBG(HKdata,method="normexp")


###################################################
### code chunk number 15: calculateMedians (eval = FALSE)
###################################################
## ## To run the serial dilution curve algorithm it is neccessary
## ## to aggregate replicate spots first.
## dataDir <- system.file("data",package="RPPanalyzer")
## setwd(dataDir)
## load("ser.dil.samples.RData")
## data(ser.dil.samples)
## ser.dil_median <- sample.median(ser.dil.samples)
## ## calculate concentration (for the attributes see help pages)
## c_Values <- calcSdc(ser.dil_median,D0=2,sel=c("measurement"), 
##                     dilution="dilution")


###################################################
### code chunk number 16: calculateSdcHelp (eval = FALSE)
###################################################
## ?calcSdc


###################################################
### code chunk number 17: normalizeProtDye (eval = FALSE)
###################################################
## ## load data set
## dataDir <- system.file("data",package="RPPanalyzer")
## setwd(dataDir)
## load("HKdata.RData")
## data(HKdata)
## ## normalize
## norm_values_pd <- normalizeRPPA(HKdata,method="proteinDye",
##                                 vals="logged")


###################################################
### code chunk number 18: normalizeHK (eval = FALSE)
###################################################
## norm_values_hk <- normalizeRPPA(HKdata,method="housekeeping",
##                                 normalizer="housekeeping",vals="logged")


###################################################
### code chunk number 19: normalizeSpotwise (eval = FALSE)
###################################################
## norm_values_hk_sbs <- normalizeRPPA(HKdata,method="spotbyspot",
##                                     normalizer="housekeeping",vals="logged")


###################################################
### code chunk number 20: normalizeMedian (eval = FALSE)
###################################################
## norm_values_row <- normalizeRPPA(HKdata,method="row")


###################################################
### code chunk number 21: normalizeBCA (eval = FALSE)
###################################################
## norm_values_eV <- normalizeRPPA(HKdata,method="extValue",
##                                 useCol="concentration",vals="logged")


###################################################
### code chunk number 22: aggReplis (eval = FALSE)
###################################################
## ## all normalization methods were performed on a sample set that was 
## ## spotted in replicates (not serially diluted). 
## ## In this case you can aggregate the replicate spots after the 
## ## normalization:
## norm_data <- sample.median(norm_values_pd)


###################################################
### code chunk number 23: selSamples (eval = FALSE)
###################################################
## selectedSamples <- select.sample.group(norm_data,
##                                        params=list("replicate"=c("1")))


###################################################
### code chunk number 24: selArrays (eval = FALSE)
###################################################
## selectedData <- remove.arrays(selectedSamples,param="target",
##                               arrays2rm=c("protein","blank","housekeeping"))


###################################################
### code chunk number 25: plotTC (eval = FALSE)
###################################################
## ## load a time course data set
## dataDir <- system.file("data",package="RPPanalyzer")
## setwd(dataDir)
## load("dataII.RData")
## data(dataII)
## ## plot time course data
## plotTimeCourse(dataII,tc.identifier=c("sample","stimulation",
##                                       "inhibition","stim_concentration"),
##                plot.split="experiment",file="Timeplot.pdf",
##                arrays2rm=c("protein","Blank"),plotformat="spline")


###################################################
### code chunk number 26: plotTC2 (eval = FALSE)
###################################################
## ## pre-process the data
## dataDir <- system.file("extdata",package="RPPanalyzer")
## setwd(dataDir)
## res <- dataPreproc(dataDir=dataDir,blocks=12,spot="aushon",exportNo=2)
## ## remove arrays
## normdat_rm <- remove.arrays(res$normdat,param="target",
##                             arrays2rm=c("protein","blank"))  
## ## select samples and export data
## sel_sampels_A549 <- select.sample.group(normdat_rm,
##                                         params=list("cell_line"="A549"),
##                                         combine=F)
## write.Data(sel_sampels_A549,FileNameExtension="HGF_sample_data_A549")
## ## read selected data
## dataexpression_1 <- read.table("HGF_sample_data_A549expression.txt")
## ## use getErrorModel function
## dataexpression_2 <- getErrorModel(dataexpression_1,verbose=FALSE)
## ## use averageData function
## dataexpression_3 <- averageData(dataexpression_2,
##                                 scaling=c("slide","replicate"),
##                                 distinguish=c("cell_line","treatment"))
## ## plot time course data
## plotTimeCourseII(dataexpression_3,
##                  filename="timecourse_HGF_sample_data_A549.pdf",
##                  legpos="top",xname="time [min]",yname="signal [a.u.]",
##                  linecolor=c("red","green","blue","black","orange","grey"))


###################################################
### code chunk number 27: boxplot (eval = FALSE)
###################################################
## ## load data set
## dataDir <- system.file("data",package="RPPanalyzer")
## setwd(dataDir)
## load("dataIII.RData")
## data(dataIII)
## ## aggregate replicates
## dataIII_median <- sample.median(dataIII)
## ## draw simple boxplot and generate PDF
## simpleBoxplot(x=dataIII_median,param="rank",
##               orderGrp=c("vx","zx","yzr","rxi"),
##               file="simpleBoxplot.pdf")
## ## draw boxplot, test for (differential) expression in comparison to 
## ## control group "vx", and generate PDF
## rppa2boxplot(x=dataIII_median,param="rank",control="vx",
##              orderGrp=c("vx","zx","yzr","rxi"),
##              file="wilcoxonBoxplot.pdf")
## ## draw boxplot, test for general differences in group expressions, 
## ## and generate PDF
## rppa2boxplot(x=dataIII_median,param="rank",control=NULL,
##              orderGrp=c("vx","zx","yzr","rxi"),
##              file="kruskalBoxplot.pdf")


###################################################
### code chunk number 28: corTest (eval = FALSE)
###################################################
## ## load data set
## dataDir <- system.file("data",package="RPPanalyzer")
## setwd(dataDir)
## load("dataIII.RData")
## data(dataIII)
## ## normalize data
## n.data <- normalizeRPPA(dataIII,method="row")
## ## aggregate replicates
## cl.data <- sample.median(n.data)
## ## test correlation
## test.correlation(cl.data,param="concentration",method.cor="kendall",
##                  method.padj="BH",file="correlation_plot.pdf")


###################################################
### code chunk number 29: heatmap (eval = FALSE)
###################################################
## ## load data set
## dataDir <- system.file("data",package="RPPanalyzer")
## setwd(dataDir)
## load("dataIII.RData")
## data(dataIII)
## ## aggregate replicates
## dataIII_median <- sample.median(dataIII)
## ## plot heatmap
## rppaList2Heatmap(dataIII_median)


