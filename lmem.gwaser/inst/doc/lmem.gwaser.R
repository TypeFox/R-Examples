## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
library(lmem.gwaser)

## ---- eval = FALSE-------------------------------------------------------
#  
#  library(lmem.gwaser)
#  
#  data ("QA_geno")
#  data ("QA_map")
#  data ("QA_pheno")
#  
#  P.data <- data ("QA_pheno")
#  G.data <- data ("QA_geno")
#  map.data <- data ("QA_map")
#  
#  gwas.cross (P.data, G.data, map.data,cross='gwas', heterozygotes=FALSE)
#  
#  summary (cross.data)
#  

## ---- eval = FALSE-------------------------------------------------------
#  
#  # Marker Quality
#  mq.g.diagnostics (crossobj=cross.data,I.threshold=0.1, p.val=0.01,na.cutoff=0.1)
#  

## ---- eval = FALSE-------------------------------------------------------
#  
#  pca <- pca.analysis(crossobj=cross.data, p.val=0.05)
#  

## ---- eval = FALSE-------------------------------------------------------
#  qk.GWAS <- gwas.analysis (crossobj=cross.data, method='QK', provide.K=FALSE,
#                            covariates=pca$scores, trait='yield',
#                            threshold='Li&Ji', p=0.05,
#                            out.file='GWAS Q + K model')
#  

## ---- eval = FALSE-------------------------------------------------------
#  pcaR.GWAS <- gwas.analysis(crossobj=cross.data, method='eigenstrat',
#                             provide.K=FALSE, covariates=pca$scores,
#                             trait='yield',
#                             threshold='Li&Ji', p=0.05,
#                             out.file='GWAS PCA as Random model')

## ---- eval = FALSE-------------------------------------------------------
#  k.GWAS <- gwas.analysis(crossobj=cross.data, method='kinship',
#                            provide.K=FALSE, covariates=FALSE, trait='yield',
#                            threshold='Li&Ji', p=0.05,
#                            out.file =' GWAS K as Random model ')

## ---- eval = FALSE-------------------------------------------------------
#  data (QA_pheno)
#  P.data.1 <- QA_pheno
#  covariate <- P.data.1 [,2]
#  
#  g.GWAS <- gwas.analysis (crossobj=cross.data,
#                           method='fixed', provide.K=FALSE,
#                           covariates=covariate,
#                           trait='yield', threshold='Li&Ji', p=0.05,
#                           out.file='GWAS fixed Groups model')

## ---- eval = FALSE-------------------------------------------------------
#  naive.GWAS <- gwas.analysis(crossobj=cross.data4, method='naive',
#                              provide.K=FALSE, covariates=FALSE,
#                              trait='yield', threshold='Li&Ji',
#                               p=0.05, out.file='GWAS naive model')

