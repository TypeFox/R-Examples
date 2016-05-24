plot.NHEMOtree <-
function(x, data, which = c("Pareto", "S_Metric", "VIMS", "Tree"), vim=0,...){
  which <- match.arg(which)
  switch(which,
         Pareto = { # Pareto front   
           misclassification_rate<- sapply(x$Trees, "[[", "Misclassification")
           costs                 <- sapply(x$Trees, "[[", "Costs")  
           FP_ges                <- cbind(misclassification_rate, costs) 
           ii                    <- order(FP_ges[,1], FP_ges[,2], decreasing = FALSE)
           FP_ges                <- FP_ges[ii,]
           FP_ges                <- cbind(FP_ges, t(t(nds_rank(t(FP_ges)))))
           Misclassification_rate<- FP_ges[which(FP_ges[,3]==1),1]  # Just non-dominated individuals
           Costs                 <- FP_ges[which(FP_ges[,3]==1),2]
           par(mfrow=c(1,1))
           plot(Misclassification_rate, Costs, main="Pareto front approximation", ...)
         },
         
         S_Metric = { # Dominated hypevolumne of all generations   
           TITEL.Konv           <- "Dominated hypervolumes of all generations"
           Generations          <- 0:(length(x$S_Metrik_temp)-1)
           Dominated_hypervolume<- x$S_Metrik_temp
           par(mfrow=c(1,1))
           plot(Generations, Dominated_hypervolume, main=TITEL.Konv, ...)
         },
         
         VIMS = { # Variable importance measures
           
           if (x$method!="Wrapper"){
             par(mfrow=c(1,1))
             if (vim==1){
               Barplot_Sim(WMasse=x$VIMS, vim=1, ...) 
             } else if (vim==2){
               Barplot_Sim(WMasse=x$VIMS, vim=2, ...) 
             } else if (vim==3){
               Barplot_Sim(WMasse=x$VIMS, vim=3, ...) 
             } else if (vim==4){
               Barplot_Sim(WMasse=x$VIMS, vim=4, ...) 
             } else if (vim==5){
               Barplot_Sim(WMasse=x$VIMS, vim=5, ...) 
             } else if (vim==6){
               Barplot_Sim(WMasse=x$VIMS, vim=6, ...) 
             } else {
               par(mfrow=c(nrow(x$VIMS),1), mar=c(4,4,2,1)) 
               for (j in 1:nrow(x$VIMS)) Barplot_Sim(WMasse=x$VIMS, vim=j, ...) 
             }}else{ 
               print("No VIMs available for Wrapper!")
             }},
         
         Tree = { # Tree with lowest misclassification rate
           
           if (x$method!="Wrapper"){
             Misclassification<- sapply(x$Trees, "[[", "Misclassification")
             Costs            <- sapply(x$Trees, "[[", "Costs")       
             FP_ges           <- cbind(1:length(x$Trees), Misclassification, Costs) 
             ii               <- order(FP_ges[,2], FP_ges[,3], decreasing = FALSE)
             FP_ges           <- FP_ges[ii,]   
             Best_Tree        <- Treeplot_prep(Tree=x$Trees[[FP_ges[1,1]]], data=data)
             Titel            <- paste("Tree with lowest misclassification rate = ", round(FP_ges[1,2], 3), "%, Costs = ", round(FP_ges[1,3], 3), sep="")
             plot(Best_Tree, main=Titel, ...)
           }else{
             Titel<- paste("Tree with lowest misclassification rate = ", round(x$Best_Tree[[2]], 3), "%, Costs = ", round(x$Best_Tree[[3]], 3), sep="")
             plot(x$Best_Tree[[1]], main=Titel, ...)
             text(x$Best_Tree[[1]])
           }
         })
}
