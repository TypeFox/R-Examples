## ----echo=FALSE, message=FALSE-------------------------------------------
devtools::load_all(".")
knitr::opts_chunk$set(echo = TRUE, fig.retina=2, fig.width=7, fig.height=5)

## ---- message=FALSE------------------------------------------------------
library(radiomics)

## ---- eval = FALSE-------------------------------------------------------
#  test <- matrix(sample(1:10, 25, replace=TRUE), ncol=5)
#  calc_features(test)

## ---- echo=FALSE, results='asis', fig.cap = "First-order features."------

knitr::kable(data.frame("Feature"=c("Energy", "Entropy", "Kurtosis", 
                                    "Mean Deviation","Skewness","Uniformity",
                                    "Mean","Median","Maximum","Minimum",
                                    "Variance","Root Mean Square","Standard Deviation"),
                        "Argument"=c("calc_energy",
                                "calc_entropy",
                                "calc_kurtosis",
                                "calc_meanDeviation",
                                "calc_skewness",
                                "calc_uniformity",
                                "calc_mean",
                                "calc_median",
                                "calc_max",
                                "calc_min",
                                "calc_variance",
                                "calc_RMS",
                                "calc_sd"
            )))

## ---- eval=FALSE---------------------------------------------------------
#  calc_features(test, features = c("calc_energy", "calc_mean"))

## ---- message=FALSE------------------------------------------------------
#Load the dataset from the Hallbey tutorial:
data(hallbey)
(hbGLCM <- glcm(hallbey, angle=0, d=1))

## ---- eval=FALSE---------------------------------------------------------
#  calc_features(hbGLCM)

## ----fig.width=5---------------------------------------------------------
image(hbGLCM)

## ---- echo=FALSE, results='asis', fig.cap = "GLCM features."-------------

knitr::kable(data.frame("Feature"=c("Mean", "Variance", "Auto Correlation" ,"Cluster Prominence", "Cluster Shade" ,
                                    "Cluster Tendency", "Contrast" ,"Correlation", "Difference Entropy",
                                    "Dissimilarity", "Energy","Entropy", "Homogeneity1", "Homogeneity2", 
                                    "Inverse Difference Moment (Normalized)", "Inverse Difference Moment",
                                    "Inverse Variance", "Maximum Probability", "Sum Average", "Sum Entropy",
                                    "Sum Variance"),
                        "Argument"=c(
              "glcm_mean", "glcm_variance", "glcm_autoCorrelation",
              "glcm_cProminence", "glcm_cShade", "glcm_cTendency",
              "glcm_contrast", "glcm_correlation", "glcm_differenceEntropy",
              "glcm_dissimilarity", "glcm_energy", "glcm_entropy", 
              "glcm_homogeneity1", "glcm_homogeneity2", "glcm_IDMN",
              "glcm_IDN", "glcm_inverseVariance", "glcm_maxProb", 
              "glcm_sumAverage", "glcm_sumEntropy", "glcm_sumVariance"
            )), caption="GLCM features")

## ----echo=FALSE----------------------------------------------------------
s <- matrix(c(0,1,2,3,0,2,3,3,2,1,1,1,3,0,3,0), nrow=4, byrow=T) #from Galloway 1974
s

## ------------------------------------------------------------------------
glrlm(s, angle=0, verbose=F)

## ---- echo=FALSE, results='asis', fig.cap = "GLCM features."-------------

knitr::kable(data.frame("Feature"=c("Grey Level Non-uniformity", "High Grey Level Run Emphasis", "Long Run Emphasis", 
                                "Long Run High Grey Level Emphasis", "Long Run Low Grey Level Emphasis",
                                "Low Grey Level Run Emphasis", "Run Length Non-uniformity", "Run Percentage",
                                "Short Run Emphasis", "Short Run High Grey Level Emphasis", "Short Run Low Grey Level Emphasis"
            ),
                        "Argument"=c("glrlm_GLN", "glrlm_HGLRE", "glrlm_LRE", 
                                "glrlm_LRHGLE", "glrlm_LRLGLE",
                                "glrlm_LGLRE", "glrlm_RLN", "glrlm_RP",
                                "glrlm_SRE", "glrlm_SRHGLE", "glrlm_SRLGLE"
            )), caption="GLRLM features")

## ---- fig.width=5--------------------------------------------------------
discTumor <- discretizeImage(radiomics::tumor, n_grey=16)
image(discTumor, axes=F, col=viridis::viridis(16))

## ---- message=FALSE------------------------------------------------------
image(glszm(discTumor))

## ---- message=FALSE------------------------------------------------------
image(mglszm(tumor))

## ---- echo=FALSE, results='asis', fig.cap = "GLCM features."-------------

knitr::kable(data.frame("Feature"=c("Small Area Emphasis", "Large Area Emphasis", "Intensity Variability", 
                                "Size Zone Variance", "Zone Percentage", "Low Intensity",
                                "High Intensity Emphasis", "Low Intensity Small Area Emphasis", "High Intensity Small Area Emphasis", 
                                "Low Intensity Large Area Emphasis", "High Intensity Large Area Emphasis"
            ),
                        "Argument"=c("glszm_SAE", "glszm_LAE", "glszm_IV", 
                                "glszm_SZV", "glszm_ZP", "glszm_LIE",
                                "glszm_HIE", "glszm_LISAE", "glszm_HISAE", 
                                "glszm_LILAE", "glszm_HILAE"
            )), caption="GLSMZ and MGLSZM features")

