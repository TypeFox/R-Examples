Barplot_Sim <-
function(WMasse, vim){
  N_Varis<- ncol(WMasse)

  # Formatieren der VIMs
    Namen <- paste("X", 1:N_Varis, sep = "")
    WMassV<- c()
    for (i in 1:N_Varis) WMassV[i]<- WMasse[vim,i]
    Max   <- max(WMassV)

  # Graphik
    if (vim==1) TITEL<- "VI - Simple absolute frequency"
    if (vim==2) TITEL<- "VI - Simple relative frequency"
    if (vim==3) TITEL<- "VI - Relative frequency"
    if (vim==4) TITEL<- "VI - Linear weigthed rel. frequency"
    if (vim==5) TITEL<- "VI - Exponential weigthed rel. frequency"
    if (vim==6) TITEL<- "VI - Permutation accuracy importance"
    
    barplot(WMassV, beside = TRUE, main = TITEL, names.arg=Namen, ylim=c(0,Max),
            cex.axis=0.9, cex.names = 0.56, ylab = "Variable improtance", xlab="Features")
}
