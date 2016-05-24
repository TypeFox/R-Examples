### R code from vignette source 'MPAgenomics.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: MPAgenomics.Rnw:77-78 (eval = FALSE)
###################################################
## install.packages("MPAgenomics", repos="http://R-Forge.R-project.org")


###################################################
### code chunk number 2: MPAgenomics.Rnw:81-82 (eval = FALSE)
###################################################
## vignette("MPAgenomics")


###################################################
### code chunk number 3: MPAgenomics.Rnw:88-96 (eval = FALSE)
###################################################
## source("http://www.braju.com/R/hbLite.R")
## hbLite("sfit")
## source("http://bioconductor.org/biocLite.R")
## biocLite("affxparser")
## biocLite("DNAcopy")
## biocLite("aroma.light")
## install.packages("aroma.affymetrix")
## install.packages("aroma.cn")


###################################################
### code chunk number 4: MPAgenomics.Rnw:113-114 (eval = FALSE)
###################################################
## celPATH="/home/user/Documents/workdir/CELdata/cel"


###################################################
### code chunk number 5: MPAgenomics.Rnw:130-140 (eval = FALSE)
###################################################
## #set your working directory (replace with the appropriate path)
## workdir="/home/user/Documents/workdir"
## setwd(workdir)
## #download file
## download.file("http://www.broadinstitute.org/mpg/birdsuite/downloads/birdsuite_inputs_1.5.3.tgz",
## destfile="./CELdata.tgz")
## #untar the file
## untar("./CELdata.tgz",files="cel",exdir=".")
## #indicate the path containing .cel files
## celPATH="./cel"


###################################################
### code chunk number 6: MPAgenomics.Rnw:153-154 (eval = FALSE)
###################################################
## chipPATH="/home/user/Documents/workdir/CD_GenomeWideSNP_6_rev3/Full/GenomeWideSNP_6/LibFiles/"


###################################################
### code chunk number 7: MPAgenomics.Rnw:160-165 (eval = FALSE)
###################################################
## #unzip required files
## unzip("./genomewidesnp6_libraryfile.zip",
##  files=c("CD_GenomeWideSNP_6_rev3/Full/GenomeWideSNP_6/LibFiles/GenomeWideSNP_6.Full.cdf"),exdir=".")
## #indicate the path containing .cdf files
## chipPATH="./home/user/Documents/workdir/CD_GenomeWideSNP_6_rev3/Full/GenomeWideSNP_6/LibFiles/"


###################################################
### code chunk number 8: MPAgenomics.Rnw:195-213 (eval = FALSE)
###################################################
## #set the directory where the .cdf files are as your working directory 
## setwd(chipPATH) 	
## ##download the 3 files .ufl, .ugp, .acs
## download.file("http://www.aroma-project.org/data/annotationData/chipTypes/GenomeWideSNP_6/GenomeWideSNP_6,Full,na31,hg19,HB20110328.ufl.gz",
## destfile="GenomeWideSNP_6,Full,na31,hg19,HB20110328.ufl.gz")
## download.file("http://www.aroma-project.org/data/annotationData/chipTypes/GenomeWideSNP_6/GenomeWideSNP_6,Full,na31,hg19,HB20110328.ugp.gz",
## destfile="GenomeWideSNP_6,Full,na31,hg19,HB20110328.ugp.gz")
## download.file("http://www.aroma-project.org/data/annotationData/chipTypes/GenomeWideSNP_6/GenomeWideSNP_6,HB20080710.acs.gz", 
## destfile="GenomeWideSNP_6,HB20080710.acs.gz")
## #unzip the gz files
## #install R.utils package containing the gunzip function
## install.packages("R.utils")
## library("R.utils")
## gunzip("GenomeWideSNP_6,Full,na31,hg19,HB20110328.ufl.gz",
## destname="GenomeWideSNP_6,Full,na31,hg19,HB20110328.ufl")
## gunzip("GenomeWideSNP_6,Full,na31,hg19,HB20110328.ugp.gz",
## destname="GenomeWideSNP_6,Full,na31,hg19,HB20110328.ugp")
## gunzip("GenomeWideSNP_6,HB20080710.acs.gz",destname="GenomeWideSNP_6,HB20080710.acs")


###################################################
### code chunk number 9: MPAgenomics.Rnw:227-237 (eval = FALSE)
###################################################
## #list files to check that your path to annotation files (.cdf, .ugp, .ufl, .acs) is correctly set
## dir(chipPATH)
## # list files to check that your path to cel files (.cel) is correctly set
## dir(celPATH) 
## #set your working directory (where you have rights to write)
## setwd(workdir) 
## #normalize data (might take several hours)
## signalPreProcess(dataSetName="datatest1", chipType="GenomeWideSNP_6",
## dataSetPath=celPATH,chipFilesPath=chipPATH, path=".",
## createArchitecture=TRUE, savePlot=TRUE, tags="Full")


###################################################
### code chunk number 10: MPAgenomics.Rnw:242-245 (eval = FALSE)
###################################################
## segcall=cnSegCallingProcess("datatest1",chromosome=c(1,5),method="PELT")
## #summary of segmentation and calling process
## segcall


###################################################
### code chunk number 11: MPAgenomics.Rnw:253-255 (eval = FALSE)
###################################################
## callfiltered=filterSeg(segcall,minLength=10,minProbes=2,keptLabel=c("gain","loss"))
## head(callfiltered)


###################################################
### code chunk number 12: MPAgenomics.Rnw:260-264 (eval = FALSE)
###################################################
## dataResponse=data.frame(files=getListOfFiles("datatest1"),
## response=c(2.105092,1.442868,1.952103,1.857819,2.047897,1.654766,2.385327,2.113406))
## res=markerSelection("datatest1",dataResponse,chromosome=21:22,signal="CN",
## onlySNP=TRUE,loss="linear")


###################################################
### code chunk number 13: MPAgenomics.Rnw:276-280 (eval = FALSE)
###################################################
## #get the file names of our data-set
## files=getListOfFiles("datatest1")
## #create the data.frame linking normal and tumor files
## normalTumorArray=data.frame(normal=rep(files[1],7),tumor=files[2:8])


###################################################
### code chunk number 14: MPAgenomics.Rnw:294-297 (eval = FALSE)
###################################################
## addData(dataSetName="datatest2",dataPath=celPATH,chipType="GenomeWideSNP_6")
## signalPreProcess(dataSetName="datatest2", chipType="GenomeWideSNP_6",
## normalTumorArray=normalTumorArray, createArchitecture=FALSE, savePlot=TRUE, tags="Full")


###################################################
### code chunk number 15: MPAgenomics.Rnw:312-316 (eval = FALSE)
###################################################
## #run the segmentation
## segfracB=segFracBSignal("datatest1",chromosome=c(1,5))
## #print summary of segmentation
## segfracB


###################################################
### code chunk number 16: MPAgenomics.Rnw:382-383 (eval = FALSE)
###################################################
## createArchitecture("datatest1","GenomeWideSNP_6",celPATH,chipPATH,".",TRUE,"Full")


###################################################
### code chunk number 17: MPAgenomics.Rnw:694-701 (eval = FALSE)
###################################################
## #normal-tumor study
## CNdata2=getCopyNumberSignal("datatest2",5,normalTumorArray=normalTumorArray,TRUE)
## fracBdata2=getFracBSignal("datatest2",5,normalTumorArray=normalTumorArray)
## 
## #study without reference
## CNdata1=getCopyNumberSignal("datatest1",5,onlySNP=TRUE)	
## fracBdata1=getFracBSignal("datatest1",5)		


###################################################
### code chunk number 18: MPAgenomics.Rnw:705-709 (eval = FALSE)
###################################################
## CNdata2$chr5
## fracBdata2$chr5$tumor
## fracBdata2$chr5$normal
## fracBdata1$chr5$tumor


###################################################
### code chunk number 19: MPAgenomics.Rnw:791-794 (eval = FALSE)
###################################################
## file="GIGAS_g_GAINmixHapMapAffy2_GenomeWideEx_6_A02_31234"
## seg1=segmentationAroma("datatest1",chromosome=21:22,onlySNP=TRUE,plot=TRUE,
## listOfFiles=file,method="PELT")


###################################################
### code chunk number 20: MPAgenomics.Rnw:880-884 (eval = FALSE)
###################################################
## file="GIGAS_g_GAINmixHapMapAffy2_GenomeWideEx_6_A07_31314"
## CNdata1=getCopyNumberSignal("datatest1",20,onlySNP=TRUE,listOfFiles=file)
## copyNumber=CNdata1$chr20$GIGAS_g_GAINmixHapMapAffy2_GenomeWideEx_6_A07_31314
## position=CNdata$chr20$position


###################################################
### code chunk number 21: MPAgenomics.Rnw:888-889 (eval = FALSE)
###################################################
## seg=segmentation(copyNumber,position=position,plot=TRUE,verbose=TRUE,method="PELT")


###################################################
### code chunk number 22: MPAgenomics.Rnw:895-896 (eval = FALSE)
###################################################
## seg$segment


###################################################
### code chunk number 23: MPAgenomics.Rnw:1016-1017 (eval = FALSE)
###################################################
## seg2=cnSegCallingProcess("datatest1",chromosome=21:22,method="PELT")


###################################################
### code chunk number 24: MPAgenomics.Rnw:1092-1099 (eval = FALSE)
###################################################
## #create the segmentData object
## callobj= callingObject(copynumber=seg$signal, segmented=seg$segmented,
##  chromosome=rep(20,length(seg$signal)), position=seg$startPos, 
##  sampleNames="sample1")
## #run the calling
## call=callingProcess(callobj,nclass=3,cellularity=1,verbose=TRUE)
## call$segment


###################################################
### code chunk number 25: MPAgenomics.Rnw:1164-1166 (eval = FALSE)
###################################################
## segmentfilter=filterSeg(call$segment,keptLabel="gain")
## segmentfilter


###################################################
### code chunk number 26: MPAgenomics.Rnw:1259-1263 (eval = FALSE)
###################################################
## dataResponse=data.frame(files=getListOfFiles("datatest1"),
## response=c(2.105092,1.442868,1.952103,1.857819,2.047897,1.654766,2.385327,2.113406))
## res=markerSelection(dataSetName="datatest1",dataResponse,chromosome=21:22,signal="CN",
## onlySNP=TRUE,loss="linear")


###################################################
### code chunk number 27: MPAgenomics.Rnw:1311-1314 (eval = FALSE)
###################################################
## dataMatrix=matrix(rnorm(5000,0,0.5),nrow=50)
## dataResponse=drop(dataMatrix%*%sample(c(rep(0,90),rep(1,10))))
## res=variableSelection(dataMatrix,dataResponse,nbFolds=5,loss="linear",plot=TRUE)


