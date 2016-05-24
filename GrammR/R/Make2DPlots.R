Make2DPlots <-
function(GraphQuant){
           message(paste("Constructing 2D", substr(GraphQuant$Name,15,18), "model"));
           PlotTabs <- NULL;
           if (!is.null(GraphQuant)){
                    ## MAKE ALL THE NECESSARY PLOTS AND SAVE AS PDFs
                    ## A TOTAL OF FOUR PLOTS WILL BE CONSTRUCTED:
                    ## SILHOUETTE PLOT, TRUE MODEL(BLACK), ESTIMATED CLUSTER MODEL AND TRUE CLUSTER MODEL(OPTIONAL)
                   
                    ## DETERMINE IF IT IS A PCoA MODEL OR A MDS MODEL AND CREATE A DIRECTORY ACCORDINGLY
                    ModType <- substr(GraphQuant$Name, 15, 18);
                    if (substr(ModType,4,4) == ""){
                                    ModType <- substr(ModType,1,3);
                    }
                    DirName <- paste(ModType, "2DModel",sep="");
                    dir.create(DirName, showWarnings = FALSE);
                  
                    # THE SILHOUETTE PLOT
                    pdf(paste(DirName,"/Silhouette.pdf", sep=""));
                    plot(GraphQuant$PamClRange, GraphQuant$SilPlot,'o', pch=18, main="Silhouette plot", xlab="Number of cluster", ylab="Silhouette width") 
                    dev.off()
                    
                      
                    # THE MODEL WITH ESTIMATED CLUSTERS
                    pdf(paste(DirName, "/EstimatedClusters.pdf", sep = ""))
                    Colors <- rainbow(GraphQuant$OptimClust);
                    plot(GraphQuant$Coords[,1], GraphQuant$Coords[,2], col = Colors[GraphQuant$ClusMem], 'p', pch = 20, cex =1, xlab="X",ylab="Y",main="Estimated models showing clusters obtained from PAM");
                    dev.off()
                    
                    ## SAVE THE ESTIMATED CLUSTERS AS A TEXT FILE.
                    SaveClus <- matrix(as.vector(GraphQuant$ClusMem), nrow = length(GraphQuant$ClusMem), ncol = 1);
                    rownames(SaveClus) <- paste("Sample",  as.character(seq(1, nrow(SaveClus))), sep="");
                    write.table(GraphQuant$ClusMem, file = paste(DirName,"/ClusterMembership.txt", sep = ""), sep = ",", row.names = TRUE);
                   
                    # THE MODEL WITH TRUE CLUSTERS
                    if (!is.null(GraphQuant$TrueMem)){
                    pdf(paste(DirName, "/TrueClusters.pdf", sep = ""))
                    NClust <- length(table(GraphQuant$TrueMem));
                    Colors <- rainbow(NClust);
                    plot(GraphQuant$Coords[,1], GraphQuant$Coords[,2], col = Colors[GraphQuant$TrueMem], 'p', pch = 20, cex =1, xlab="X",ylab="Y",main="Estimated models showing clusters obtained from PAM");
                    dev.off()
                     }
        }
 }
