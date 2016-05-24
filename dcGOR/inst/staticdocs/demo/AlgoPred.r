# This is a demo for creating and using dcGO
# 
# There are two key concepts behind dcGO. The first concept is to label protein domains with ontology, for example, with Gene Ontology. That is why it is called dcGO, domain-centric Gene Ontology. The second concept is to use ontology-labeled protein domains for, for example, protein function prediction. Put it in a simple way, the first concept is about how to create dcGO resource, and the second concept is about how to use dcGO resource.
# Accordingly, this demo contains two parts. The first part is to creat/label SCOP domains/supra-domains with Human Phenotype (HP) Ontology (actually a subontology called Phenotypic Abnormality: HPPA). The second part is to use the built resource for phenotype predictions in human genes. 
# The demo assumes that the user has collected two bits (files) of information:
## 1) anno.file: an annotation file containing annotations between proteins/genes and ontology terms. For example, a file containing annotations between human genes and HP terms can be found in <a href="http://dcgor.r-forge.r-project.org/data/Algo/HP_anno.txt">HP_anno.txt</a>, with two columns (1st for 'SeqID', 2nd for 'termID')
## 2) architecture.file: an architecture file containing domain architectures (including individual domains) for proteins/genes. For example, a file containing human genes and domain architectures can be found in <a href="http://dcgor.r-forge.r-project.org/data/Algo/SCOP_architecture.txt">SCOP_architecture.txt</a>, with two columns (1st for 'SeqID', 2nd for 'Architecture', that is, a string of comma-separated SCOP domains)
#
###############################################################################
library(dcGOR)


#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
# Creat domains/supra-domains with Human Phenotype Phenotypic Abnormality (HPPA) terms
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------

## 1) prepare input and output files
### input files
anno.file <- "http://dcgor.r-forge.r-project.org/data/Algo/HP_anno.txt"
architecture.file <- "http://dcgor.r-forge.r-project.org/data/Algo/SCOP_architecture.txt"
### output files
output.file <- "sf2HPPA.txt"
RData.HIS.customised <- "sf2HPPA.HIS.RData"
## 2) apply dcGO algorithm to infer domain-centric ontology
res <- dcAlgo(anno.file, architecture.file, output.file, ontology="HPPA", feature.mode="supra", fdr.cutoff=0.05, parallel=FALSE)
res[1:5,]
## 3) propagate ontology annotations
res_RData <- dcAlgoPropagate(input.file=output.file, ontology="HPPA", output.file=RData.HIS.customised)
## In your working directory, you should see a RData-formatted file: "sf2HPPA.HIS.RData". It is an object of S3 class 'HIS' (see <a href="http://supfam.org/dcGOR/dcAlgoPropagate.html">dcAlgoPropagate</a> for details). It contains created domain-centric annotations, which will be used subsequently for prediction.
names(res_RData)

#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
# Predict human genes with Human Phenotype Phenotypic Abnormality (HPPA) terms
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------

## 1) prepare input and output files
### input file
architecture.file <- "http://dcgor.r-forge.r-project.org/data/Algo/SCOP_architecture.txt"
### output file
prediction.file <- "SCOP_architecture.HPPA_predicted.txt"
## 2) do predictions
res_Pred <- dcAlgoPredictMain(input.file=architecture.file, output.file=prediction.file, feature.mode="supra", parallel=FALSE, RData.HIS.customised=RData.HIS.customised)
res_Pred[1:5,]
## 3) look at the prediction performance via Precision-Recall (PR) analysis
GSP.file <- "http://dcgor.r-forge.r-project.org/data/Algo/HP_anno.txt"
res_PR <- dcAlgoPredictPR(GSP.file=GSP.file, prediction.file=prediction.file, ontology="HPPA")
res_PR
### plot PR curve
plot(res_PR[,2], res_PR[,1], xlim=c(0,1), ylim=c(0,1), type="b", xlab="Recall", ylab="Precision")

#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
# Compared to RWR-based predictions and naive predictions
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------

## 1) RWR-based predictions, please see: dcRWRpredict
## 2) Naive predictions, please see: dcNaivePredict
