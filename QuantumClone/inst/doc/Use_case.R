## ----QuantumCat_example, eval = FALSE------------------------------------
#      # Example was generated calling:
#     Input_Example<-QuantumCat(number_of_clones = 4,
#                               number_of_mutations = 100,
#                               ploidy = "AB",depth = 150,
#                               number_of_samples = 2,
#                               contamination = c(0,0))

## ----echo=FALSE----------------------------------------------------------
  knitr::kable(head(QuantumClone::Input_Example[[1]]))

## ----One_step_example, eval = FALSE--------------------------------------
#    One_step_clustering(SNV_list, FREEC_list = NULL, contamination,
#    nclone_range = 2:5, clone_priors = NULL, prior_weight = NULL,
#    maxit = 1 , preclustering = T, simulated = F, epsilon = 5 * (10^(-3)),
#    save_plot = T, ncores = 1, plot_3D = F, plot_3D_before_clustering = F,
#    restrict.to.AB = F, output_directory = NULL)

## ----out,echo= FALSE-----------------------------------------------------
  knitr::kable(head(QuantumClone::QC_output$filtered.data[[1]]))

## ----plot, echo= FALSE,warning=FALSE-------------------------------------
  QuantumClone::plot_QC_out(QuantumClone::QC_output)

## ----margin, echo= FALSE,warning=FALSE-----------------------------------
  QuantumClone::plot_with_margins_densities(QuantumClone::QC_output)

## ----evol, echo = FALSE, warning = FALSE---------------------------------
QuantumClone::evolution_plot(QuantumClone::QC_output,Sample_names = c("Timepoint_1","Timepoint_2"))


## ----Tree, echo = TRUE, warning = FALSE,eval=TRUE------------------------
Cellularities<-cbind(QuantumClone::QC_output$EM.output$centers[[1]],QuantumClone::QC_output$EM.output$centers[[2]])
Tree<-QuantumClone::Tree_generation(Cellularities)


## ----ShowTree, echo = FALSE----------------------------------------------
knitr::kable(Tree[[1]][[1]])

## ----TreePlot, echo = TRUE, warning = FALSE,eval=TRUE--------------------
QuantumClone::multiplot_trees(Tree,d = 4)


