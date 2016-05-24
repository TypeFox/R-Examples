## ----installAndLoadPackages,eval=FALSE,echo=TRUE-------------------------
#  
#  #Bioconductor package installation
#  source("http://bioconductor.org/biocLite.R")
#  biocLite(c("GEOquery","Biobase"))
#  install.packages("MM2S", repos="http://cran.r-project.org")
#  
#  #CDF installation
#  download.file(
#    url = "http://mbni.org/customcdf/20.0.0/entrezg.download/hgu133plus2hsentrezgcdf_20.0.0.tar.gz",
#                method = "auto",destfile = "hgu133plus2hsentrezgcdf_20.0.0.tar.gz")
#  install.packages("hgu133plus2hsentrezgcdf_20.0.0.tar.gz",type = "source",repos=NULL)
#  
#  #Modified Affy package
#  download.file(
#    url = "http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/20.0.0/affy_1.48.0.tar.gz",
#                method = "auto",destfile = "affy_1.48.0.tar.gz")
#  install.packages("affy_1.48.0.tar.gz",type = "source",repos=NULL)

## ----Load Libraries,eval=FALSE-------------------------------------------
#  suppressPackageStartupMessages(library(MM2S))
#  suppressPackageStartupMessages(library(affy))
#  suppressPackageStartupMessages(library(Biobase))
#  suppressPackageStartupMessages(library(GEOquery))
#  suppressPackageStartupMessages(library(hgu133plus2hsentrezgcdf))

## ----getDataFromGEO,eval=FALSE-------------------------------------------
#  gse<-getGEOSuppFiles(GEO = "GSE37418")
#  untar(tarfile = "./GSE37418/GSE37418_RAW.tar",exdir = "CelFiles")

## ----cleanAndNormalize,eval=FALSE----------------------------------------
#  
#  # Generate the Affy Expression Object
#  affyRaw <- ReadAffy(celfile.path = "CelFiles",verbose = F,
#                      cdfname="hgu133plus2hsentrezgcdf",compress = T)
#  
#  # View object
#  affyRaw
#  
#  #Perform Data Background Correction and Normalization
#  eset <- expresso(affyRaw,bgcorrect.method="rma",normalize.method="quantiles",
#                   pmcorrect.method="pmonly",summary.method="medianpolish",verbose = FALSE)
#  #Obtain the Microarray Expression Dataset
#  datamatrix<-exprs(eset)
#  
#  # Polish the rownames (remove the _at from the Entrez IDs)
#  rownames(datamatrix)<-gsub(rownames(datamatrix),pattern="_at",replacement="")
#  
#  # Create a new variable representing the cleaned microarray data that will be used in MM2S
#  ExprMatrix<-datamatrix

## ----findHumanModelSubtypes,eval=FALSE-----------------------------------
#  # Conduct Subtype Predictions the samples, save results in a XLS file
#  HumanPreds<-MM2S.human(InputMatrix=ExprMatrix[,1:10],xls_output=FALSE,parallelize=4)

## ----GeneratePredictionHeatmap,echo=TRUE,eval=FALSE----------------------
#  # Now generate a heatmap of the predictions and save the results in a PDF file.
#  # This indicates MM2S confidence perdictions for each sample .
#  # We view the samples here.
#  PredictionsHeatmap(InputMatrix=HumanPreds$Predictions,pdf_output=TRUE,pdfheight=12,pdfwidth=10)

## ----InstallingFromGithubExample,echo=TRUE-------------------------------
# library(Biobase)
# library(devtools)
# install_github(repo="DGendoo/MM2S")
# install_github(repo="DGendoo/MM2Sdata")

## ----sessionInfo,echo=FALSE,results="asis",eval=FALSE--------------------
#  utils::toLatex(sessionInfo())

