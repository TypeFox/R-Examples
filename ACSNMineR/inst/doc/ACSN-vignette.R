## ----gmt_example---------------------------------------------------------
# Retrieve path of the example gmt
file<-system.file("extdata", "cellcycle_short.gmt", package = "ACSNMineR")
# Then import it
gmt<-ACSNMineR::format_from_gmt(file)

## ----gmt_display, echo = FALSE-------------------------------------------
knitr::kable(gmt[,1:10])

## ----gmt_map_code, eval = FALSE------------------------------------------
#  # Name of available maps:
#  names(ACSNMineR::ACSN_maps)
#  

## ----gmt_map_show, echo = FALSE------------------------------------------
knitr::kable(names(ACSNMineR::ACSN_maps))

## ----gmt_access_code, eval = FALSE---------------------------------------
#  #And accessing them:
#  ACSNMineR::ACSN_maps$CellCycle

## ----gmt_access_display, echo = FALSE------------------------------------
knitr::kable(head(ACSNMineR::ACSN_maps$CellCycle[,1:10]), row.names = FALSE)

## ----test_genes,eval = FALSE---------------------------------------------
#  ACSNMineR::genes_test

## ----test_genes_show, echo = FALSE---------------------------------------
knitr::kable(ACSNMineR::genes_test)

## ----Analysis_code-------------------------------------------------------
Example<-ACSNMineR::enrichment(ACSNMineR::genes_test,
    min_module_size = 10, 
    threshold = 0.05,
    maps = list(cellcycle = ACSNMineR::ACSN_maps$CellCycle))

## ----Analysis_show, echo = FALSE-----------------------------------------
knitr::kable(Example,row.names = FALSE)

## ----multi_analysis_code-------------------------------------------------
Example<-ACSNMineR::multisample_enrichment(Genes_by_sample = list(set1 = ACSNMineR::genes_test[-1],
                                                              set2 = ACSNMineR::genes_test[-2]),
    maps = ACSNMineR::ACSN_maps$CellCycle,
    min_module_size = 10,
    cohort_threshold = FALSE)

## ----multi_ana_code_1, eval = FALSE--------------------------------------
#  print(Example[[1]])

## ----multi_ana_show_1, echo = FALSE--------------------------------------
knitr::kable(Example[[1]],row.names = FALSE)

## ----multi_ana_code_2, eval = FALSE--------------------------------------
#  print(Example[[2]])

## ----multi_ana_show_2, echo = FALSE--------------------------------------
knitr::kable(Example[[2]],row.names = FALSE)


## ----heatmap_plot--------------------------------------------------------
ACSNMineR::represent_enrichment(enrichment = list(
    SampleA = ACSNMineR::enrichment_test[1:10,], 
    SampleB = ACSNMineR::enrichment_test[3:10,]),
    plot = "heatmap", 
    scale = "reverselog",
    low = "steelblue" , high ="white",
    na.value = "grey")

## ----barplot, warning = FALSE--------------------------------------------
ACSNMineR::represent_enrichment(enrichment = list(
    SampleA = ACSNMineR::enrichment_test[1:10,], 
    SampleB = ACSNMineR::enrichment_test[3:10,]),
    plot = "bar", 
    scale = "log")


