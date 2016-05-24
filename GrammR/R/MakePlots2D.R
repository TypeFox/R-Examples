MakePlots2D <- function(GraphQuant){
           PlotTabs <- NULL;
           if (!is.null(GraphQuant)){
                    ## MAKE ALL THE NECESSARY PLOTS AND PRESENT THEM IN A NEW WINDOW
                    ## THE WINDOW WILL HAVE FOUR(WITH AN OPTIONAL FIFTH) WINDOWS IN THE FOLLOWING ORDER:
                    ## SUMMARY, SILHOUETTE PLOT, TRUE MODEL(BLACK), ESTIMATED CLUSTER MODEL AND TRUE CLUSTER MODEL(OPTIONAL)

                    PlotTabs <- gnotebook(tab.pos = 1, closebuttons = FALSE);
                    Tab <- list();
                    
                    ModType <- substr(GraphQuant$Name, 15, 18);
                     if (substr(ModType,4,4) == " "){
                            ModType <- substr(ModType,1,3);
                    }
                    DirName <- paste(ModType, "2DModel",sep="");
                    dir.create(DirName, showWarnings = FALSE);
                    
                    # THE SUMMARY TAB
                    Tab$Summary <- gtext(label = "Summary", Summarize(GraphQuant), container = PlotTabs);

                    # THE SILHOUETTE PLOT TAB
                    SilFile <- paste(getwd(),"/", DirName, "/Sihouette.jpg", sep = "");
                    jpeg(SilFile, width = 600, height = 600, quality = 100)
                    plot(GraphQuant$PamClRange, GraphQuant$SilPlot,'o', pch=18, main="Silhouette plot", xlab="Number of cluster", ylab="Silhouette width") 
                    dev.off()
                    Tab$Silhouette <- gimage(label = "Silhouette Plot", filename = SilFile, container = PlotTabs)

                   ## THE MODEL WITH ESTIMATED CLUSTERS TAB
                    ClustFile <- paste(getwd(),"/", DirName,"/ClusterModel.jpg", sep = "");
                    jpeg(ClustFile, width = 600, height = 600, quality = 100)
                    Colors <- rainbow(GraphQuant$OptimClust);
                    plot(GraphQuant$Coords[,1], GraphQuant$Coords[,2], col = Colors[GraphQuant$ClusMem], 'p', pch = 20, cex =1, xlab="X",ylab="Y",main="Estimated models showing clusters obtained from PAM");
                    dev.off()
                    Tab$EstClust <- gimage(label = "Estimated Clusters", filename = ClustFile, container = PlotTabs);
                    
                    ## SAVE THE ESTIMATED CLUSTERS AS A TEXT FILE.
                    SaveClus <- matrix(as.vector(GraphQuant$ClusMem), nrow = length(GraphQuant$ClusMem), ncol = 1);
                    rownames(SaveClus) <- paste("Sample",  as.character(seq(1, nrow(SaveClus))), sep="");
                    write.table(GraphQuant$ClusMem, file = paste(DirName,"/ClusterMembership.txt", sep = ""), sep = ",", row.names = TRUE);
                    
                    # THE MODEL WITH TRUE CLUSTERS TAB
                    if (!is.null(GraphQuant$TrueMem)){
                            TrueFile <- paste(getwd(),"/", DirName,"/TrueModel.jpg", sep = "");
                            jpeg(TrueFile, width = 600, height = 600, quality = 100);
                            NClust <- length(table(GraphQuant$TrueMem));
                            Colors <- rainbow(NClust);
                            plot(GraphQuant$Coords[,1], GraphQuant$Coords[,2], col = Colors[GraphQuant$TrueMem], 'p', pch = 20, cex =1, xlab="X",ylab="Y",main="Estimated models showing true clusters");
                            dev.off()
                            Tab$EstClust <- gimage(label = "True Clusters", filename = TrueFile, container = PlotTabs);
                    }

                    svalue(PlotTabs) <- 1;         
          }
          return(PlotTabs);
}