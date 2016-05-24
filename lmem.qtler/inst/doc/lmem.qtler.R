## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
library(lmem.qtler)

## ---- eval = FALSE-------------------------------------------------------
#  library(lmem.qtler)
#  data (DHpop_pheno)
#  data (DHpop_geno)
#  data (DHpop_map)
#  G.data <- DHpop_geno
#  map.data <- DHpop_map
#  P.data <- DHpop_pheno
#  
#  cross.data <- qtl.cross (P.data, G.data, map.data, cross='dh',
#                           heterozygotes=FALSE)
#  summary (cross.data)
#  
#  ## Pheno Quality
#  pq.diagnostics (crossobj=cross.data)
#  
#  ## Marker Quality
#  mq.diagnostics (crossobj=cross.data,I.threshold=0.1,
#                p.val=0.01,na.cutoff=0.1)
#  

## ---- eval = FALSE-------------------------------------------------------
#  ### QTL_SMA
#  QTL.result <- qtl.analysis (crossobj=cross.data,step=0, method='SIM',
#                              trait="height", threshold="Li&Ji", distance=30,
#                              cofactors=NULL, window.size=30)
#  
#  ### QTL_SIM
#  QTL.result <- qtl.analysis ( crossobj=cross.data, step=5,
#                               method='SIM',trait="height", threshold="Li&Ji",
#                               distance=30,cofactors=NULL, window.size=30)

## ---- eval = FALSE-------------------------------------------------------
#  ### QTL CIM
#  cofactors <- as.vector (QTL.result$selected$marker)
#  QTL.result <- qtl.analysis ( crossobj=cross.data, step=5,
#                               method='CIM', trait="height", threshold="Li&Ji",
#                               distance=30, cofactors=cofactors, window.size=30)
#  

## ---- eval = FALSE-------------------------------------------------------
#  data (SxM_geno)
#  data (SxM_map)
#  data (SxMxE_pheno)
#  
#  P.data <- SxMxE_pheno
#  G.data <- SxM_geno
#  map.data <- SxM_map
#  
#  cross.data <- qtl.cross (P.data, G.data, map.data, cross='dh',
#                            heterozygotes=FALSE)
#  
#  summary (cross.data)
#  jittermap(cross.data)
#  
#  Pheno Quality
#  pq.diagnostics (crossobj=cross.data, boxplot=FALSE)
#  
#  Marker Quality
#  mq.diagnostics (crossobj=cross.data,I.threshold=0.1,
#                p.val=0.01,na.cutoff=0.1)
#  
#  ### QTL_SIM
#  QTL.result <- qtl.memq (crossobj = cross.data, P.data = P.data,
#                         env.label = c('ID91','ID92','MAN92','MTd91',
#                          'MTd92','MTi91','MTi92','SKs92','WA91','WA92'),
#                          trait = 'yield', step = 10, method = 'SIM',
#                          threshold = 'Li&Ji', distance = 50, cofactors = NULL,
#                          window.size = 50)
#  
#  ### QTL_CIM
#  QTL.result <- qtl.memq (crossobj = cross.data, P.data = P.data,
#                          env.label = c('ID91','ID92','MAN92','MTd91','MTd92',
#                          'MTi91','MTi92','SKs92','WA91','WA92'),
#                          trait = 'yield', step = 10, method = 'CIM',
#                          threshold = 'Li&Ji', distance = 50,
#                          cofactors = QTL.result$selected$marker, window.size = 50)
#  

