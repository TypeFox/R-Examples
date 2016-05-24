## ----echo = F------------------------------------------------------------
library(knitr)
options(digits = 3)

## ----echo = F------------------------------------------------------------
options(digits = 3)

## ----results='asis'------------------------------------------------------
library(codyn)
library(knitr)
data(knz_001d)
kable(head(knz_001d))

## ----results='asis'------------------------------------------------------
KNZ_stability <- community_stability(knz_001d, 
                                   time.var = "year", 
                                   abundance.var = "abundance", 
                                   replicate.var = "subplot")
kable(head(KNZ_stability))

## ----results='asis'------------------------------------------------------
KNZ_A1_stability <- community_stability(df = subset(knz_001d, subplot=="A_1"),  
                                      time.var = "year", 
                                      abundance.var = "abundance")
KNZ_A1_stability

## ----results='asis', eval=FALSE------------------------------------------
#  KNZ_variance_ratio <- variance_ratio(df = knz_001d,
#                                     species.var = "species",
#                                     time.var = "year",
#                                     abundance.var = "abundance",
#                                     bootnumber = 10,
#                                     replicate.var = "subplot")
#  
#  kable(KNZ_variance_ratio)

## ----results='asis', eval=FALSE------------------------------------------
#  KNZ_variance_ratio_avgrep <- variance_ratio(knz_001d,
#                                            time.var = "year",
#                                            species.var = "species",
#                                            abundance.var = "abundance",
#                                            bootnumber = 10,
#                                            replicate.var = "subplot",
#                                            average.replicates = FALSE)
#  
#  kable(head(KNZ_variance_ratio_avgrep))

## ----results='asis', eval=FALSE------------------------------------------
#  KNZ_A1_variance_ratio <- variance_ratio(df = subset(knz_001d, subplot=="A_1"),
#                                        time.var = "year",
#                                        species.var = "species",
#                                        abundance.var = "abundance",
#                                        bootnumber = 10)
#  kable(KNZ_A1_variance_ratio)

## ----results='asis'------------------------------------------------------
KNZ_synchrony_Loreau <- synchrony(df = knz_001d, 
                         time.var = "year", 
                         species.var = "species", 
                         abundance.var = "abundance", 
                         replicate.var = "subplot")
kable(head(KNZ_synchrony_Loreau))

## ----results='asis'------------------------------------------------------
KNZ_A1_synchrony_Loreau <- synchrony(df = subset(knz_001d, subplot=="A_1"),
                            time.var = "year",
                            species.var = "species", 
                            abundance.var = "abundance")
KNZ_A1_synchrony_Loreau

## ----results='asis'------------------------------------------------------
KNZ_synchrony_Gross <- synchrony(df = knz_001d, 
                           time.var = "year", 
                           species.var = "species",  
                           abundance.var = "abundance", 
                           metric = "Gross", 
                           replicate.var = "subplot")
kable(head(KNZ_synchrony_Gross))

## ----results='asis'------------------------------------------------------
KNZ_A1_synchrony_Gross <- synchrony(df = subset(knz_001d, subplot=="A_1"),
                              time.var = "year", 
                              species.var = "species",  
                              abundance.var = "abundance", 
                              metric = "Gross")
KNZ_A1_synchrony_Gross

