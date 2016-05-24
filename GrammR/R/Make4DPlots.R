Make4DPlots <-
function(GraphQuant){
         message(paste("Constructing 4D", substr(GraphQuant$Name,15,18), "model"));
        PlotTabs <- NULL;
        if (!is.null(GraphQuant)){
                    ## MAKE ALL THE NECESSARY PLOTS AND SAVE IN SUBDIRECTORIES.
                    ## THREE SUBDIRECTORIES WILL BE CREATED ALONG WITH A SILHOUETTE PDF PLOT.
                    ## TRUE MODEL(BLACK), ESTIMATED CLUSTER MODEL AND TRUE CLUSTER MODEL(OPTIONAL)
                                       
                    ## DETERMINE IF IT IS A PCoA MODEL OR A MDS MODEL AND CREATE A DIRECTORY ACCORDINGLY.
                         ModType <- substr(GraphQuant$Name, 15, 18);
                    if (substr(ModType,4,4) == ""){
                                    ModType <- substr(ModType,1,3);
                   }
                    DirName <- paste(ModType, "4DModel",sep="");
                    dir.create(DirName, showWarnings = FALSE);
                                        
                   
                    # THE SILHOUETTE PLOT TAB
                    pdf(paste(DirName,"/Silhouette.pdf", sep=""));
                    plot(GraphQuant$PamClRange, GraphQuant$SilPlot,'o', pch=18, main="Silhouette plot", xlab="Number of cluster", ylab="Silhouette width") 
                    dev.off()
                    
                   ## THE MODEL WITH ESTIMATED CLUSTERS TAB
                     
                   ## CREATE A SUBDIRECTORY FOR THE WEBPAGES
                   DirName_Clust <- paste(DirName, "/EstimatedClusters", sep="");
                   dir.create(DirName_Clust, showWarnings = FALSE);    
                   Colors <- rainbow(GraphQuant$OptimClust);                           
                    
                    ## Tab$ClustModel TAB CONTAINS 4 TABS, EACH WITH PLOTS CONSTRUCTED BY REMOVING ONE DIMENSION AT A TIME.
                    open3d(useNULL = TRUE) #Estimated cluster MDS plot showing XYZ axis
                    plot3d(GraphQuant$Coords[,-4], col = Colors[GraphQuant$ClusMem], type='s',size=1,xlab="X",ylab="Y",zlab="Z",main="Cluster model");
                    writeWebGL(dir = DirName_Clust, filename = file.path(DirName_Clust, "Model_XYZ.html"), width=800, height=800, snapshot = FALSE)
                    rgl.close();
                  
                    open3d(useNULL = TRUE) #Estimated cluster MDS plot showing YZW axis
                    plot3d(GraphQuant$Coords[,-1], col = Colors[GraphQuant$ClusMem], type='s',size=1,xlab="Y",ylab="Z",zlab="W",main="Cluster model");
                    writeWebGL(dir = DirName_Clust, filename = file.path(DirName_Clust, "Model_YZW.html"), width=800, height=800, snapshot = FALSE)
                    rgl.close();
                  
                    open3d(useNULL = TRUE) #Estimated cluster MDS plot showing XZW axis
                    plot3d(GraphQuant$Coords[,-2], col = Colors[GraphQuant$ClusMem], type='s',size=1,xlab="X",ylab="Z",zlab="W",main="Cluster model");
                    writeWebGL(dir = DirName_Clust, filename = file.path(DirName_Clust, "Model_XZW.html"), width=800, height=800, snapshot = FALSE)
                    rgl.close();
                  
                    open3d(useNULL = TRUE) #Estimated cluster MDS plot showing XYW axis
                    plot3d(GraphQuant$Coords[,-3], col = Colors[GraphQuant$ClusMem], type='s',size=1,xlab="X",ylab="Y",zlab="W",main="Cluster model");
                    writeWebGL(dir = DirName_Clust, filename = file.path(DirName_Clust, "Model_XYW.html"), width=800, height=800, snapshot = FALSE)
                    rgl.close();
                    
                    
                    ## SAVE THE ESTIMATED CLUSTERS AS A TEXT FILE.
                    SaveClus <- matrix(as.vector(GraphQuant$ClusMem), nrow = length(GraphQuant$ClusMem), ncol = 1);
                    rownames(SaveClus) <- paste("Sample",  as.character(seq(1, nrow(SaveClus))), sep="");
                    write.table(GraphQuant$ClusMem, file = paste(DirName,"/ClusterMembership.txt", sep = ""), sep = ",", row.names = TRUE);
                    
                    # THE MODEL WITH THE TRUE CLUSTERS TAB
                    if (!is.null(GraphQuant$TrueMem)){
                  
                     ## CREATE A SUBDIRECTORY FOR THE WEBPAGES
                    
                    DirName_True <- paste(DirName, "/TrueModel", sep="");
                    dir.create(DirName_True, showWarnings = FALSE);           
                    NClust <- length(table(GraphQuant$TrueMem));
                    Colors <- rainbow(NClust);
                    
                    open3d(useNULL = TRUE) #Estimated cluster MDS plot showing XYZ axis
                    plot3d(GraphQuant$Coords[,-4], col = Colors[GraphQuant$TrueMem], type='s',size=1,xlab="X",ylab="Y",zlab="Z",main="True Clusters");
                    writeWebGL(dir = DirName_True, filename = file.path(DirName_True, "Model_XYZ.html"), width=800, height=800, snapshot = FALSE)
                    rgl.close();
          
                    open3d(useNULL = TRUE) #Estimated cluster MDS plot showing YZW axis
                    plot3d(GraphQuant$Coords[,-1], col = Colors[GraphQuant$TrueMem], type='s',size=1,xlab="Y",ylab="Z",zlab="W",main="True Clusters");
                    writeWebGL(dir = DirName_True, filename = file.path(DirName_True, "Model_YZW.html"), width=800, height=800, snapshot = FALSE)
                    rgl.close();
          
                    open3d(useNULL = TRUE) #Estimated cluster MDS plot showing XZW axis
                    plot3d(GraphQuant$Coords[,-2], col = Colors[GraphQuant$TrueMem], type='s',size=1,xlab="X",ylab="Z",zlab="W",main="True Clusters");
                    writeWebGL(dir = DirName_True, filename = file.path(DirName_True, "Model_XZW.html"), width=800, height=800, snapshot = FALSE)
                    rgl.close();
                    
                    open3d(useNULL = TRUE) #Estimated cluster MDS plot showing XYW axis
                    plot3d(GraphQuant$Coords[,-3], col = Colors[GraphQuant$TrueMem], type='s',size=1,xlab="X",ylab="Y",zlab="W",main="True Clusters");
                    writeWebGL(dir = DirName_True, filename = file.path(DirName_True, "Model_XYW.html"), width=800, height=800, snapshot = FALSE)
                    rgl.close();
                 }
        }   
}
