## ----installAndLoadPackages,eval=TRUE------------------------------------
#install.packages("MM2S", repos="http://cran.r-project.org")
suppressPackageStartupMessages(library(MM2S))
#install.packages("MM2Sdata", repos="http://cran.r-project.org")
suppressPackageStartupMessages(library(MM2Sdata))

## ----findMouseModelSubtypes----------------------------------------------
data(GSE36594Expr)
ExprMat<-exprs(GSE36594Expr)
GTML<-ExprMat[,grep("GTML_MB",(colnames(exprs(GSE36594Expr))))]

#Change mouse sample names for clarity
for(sample in 1:ncol(GTML))
{
  newnames<-strsplit(x=(colnames(GTML)[sample]),split="_")[[1]][1]
  colnames(GTML)[sample]<-newnames
}

# Conduct Subtype Predictions for those particular replicates, save results in a XLS file
GTMLPreds<-MM2S.mouse(InputMatrix=GTML,xls_output=TRUE,parallelize=1)

## ----GeneratePredictionHeatmap,echo=TRUE---------------------------------
# Now generate a heatmap of the predictions and save the results in a PDF file.  
# This indicates MM2S confidence perdictions for each sample replicate of the GMTL model. 
# We view the samples here. 
PredictionsHeatmap(InputMatrix=GTMLPreds$Predictions[1:20,],pdf_output=TRUE,pdfheight=12,pdfwidth=10)

# NB: Output may appear on multiple pages

## ----GeneratePredictionBarplot,eval=FALSE,echo=TRUE----------------------
#  # To run the function all the GTML sample replicates, please run:
#  # PredictionsBarplot(InputMatrix=GTMLPreds$Predictions[1:20,],pdf_output=TRUE,pdfheight=5,pdfwidth=12)
#  # NB: Output may appear on multiple pages

## ----PredictionDistributionPie,echo=TRUE---------------------------------
PredictionsDistributionPie(InputMatrix=GTMLPreds,pdf_output=TRUE,pdfheight=5,pdfwidth=5)

## ----PredictionDistributionBoxplot,echo=TRUE-----------------------------
PredictionsDistributionBoxplot(InputMatrix=GTMLPreds,pdf_output=FALSE)

## ----PCARenderingOfPredictions,echo=TRUE---------------------------------
PCARender(GSVAmatrixTesting=GTMLPreds$RankMatrixTesting, 
          GSVAmatrixTraining=GTMLPreds$RankMatrixTraining)

## ----findHumanModelSubtypes----------------------------------------------
data(GSE37418Expr)
HumanExpr<-exprs(GSE37418Expr)
# Conduct Subtype Predictions for all samples, save results in a XLS file
# [This will take a few minutes to compute]
HumanPreds<-MM2S.human(InputMatrix=HumanExpr,xls_output=TRUE,parallelize=1)

## ----ComparePredictions,echo=TRUE----------------------------------------
# We first assess the distribution of the known subtypes for the 76 samples.
table(pData(GSE37418Expr)$characteristics_ch1)
# We now assess the distribtuion of MM2S predicted subtypes for the 76 samples. 
table(HumanPreds$MM2S_Subtype[,2])
# Side-by-side comparison of MM2S predictions and pre-determined subtypes across all samples
# first check that all samples are matching in the pData and MM2S
all(HumanPreds$MM2S_Subtype[,1] == rownames(pData(GSE37418Expr)))
# then generate comparisons
ComparisonTable<-cbind(Sample=rownames(pData(GSE37418Expr)),
                       Original=as.character(pData(GSE37418Expr)$characteristics_ch1),MM2S=HumanPreds$MM2S_Subtype[,2])
# We view the first 15 samples here
ComparisonTable[1:10,]

## ----GeneratePredictionHeatmapAndPCARendering,echo=TRUE------------------
# Now generate a heatmap of the predictions and save the results in a PDF file.  
# This indicates MM2S confidence perdictions for each sample. 
# We can view the first 10 samples. 
PredictionsHeatmap(InputMatrix=HumanPreds$Predictions[1:10,],pdf_output=TRUE,pdfheight=10,pdfwidth=5)

# NB: Output may appear on multiple pages

# We can graphically visualize different sample replicates and their nearest human MB neighbors 
# from the MM2S training set using Principal Component Analysis (PCA). 
PCARender(GSVAmatrixTesting=HumanPreds$RankMatrixTesting, 
          GSVAmatrixTraining=HumanPreds$RankMatrixTraining)

## ----InstallingFromGithubExample,echo=TRUE-------------------------------
# library(Biobase)
# library(devtools)
# install_github(repo="DGendoo/MM2S")
# install_github(repo="DGendoo/MM2Sdata")

## ----sessionInfo,echo=FALSE,results="asis"-------------------------------
utils::toLatex(sessionInfo())

