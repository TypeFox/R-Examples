MakePlots4D <- function(GraphQuant){
           PlotTabs <- NULL;
           if (!is.null(GraphQuant)){
                    ## MAKE ALL THE NECESSARY PLOTS AND PRESENT THEM IN A NEW WINDOW
                    ## THE WINDOW WILL HAVE FOUR(WITH AN OPTIONAL FIFTH) WINDOWS IN THE FOLLOWING ORDER:
                    ## SUMMARY, SILHOUETTE PLOT, TRUE MODEL(BLACK), ESTIMATED CLUSTER MODEL AND TRUE CLUSTER MODEL(OPTIONAL)
                    
                    ## FOR 4D MODELS, THE ESTIMATED MODEL, ESTIMATED CLUSTER MODEL AND TRUE CLUSTER MODEL TABS WILL CONTAIN 
                    ## HTML LINKS WHICH WHEN CLICKED WILL BE OPENED IN A WEB BROWSER.

                    PlotTabs <- gnotebook(tab.pos = 3, closebuttons = FALSE);
                    Tab <- list();
                    
                    ## DETERMINE IF IT IS A PCoA MODEL OR A MDS MODEL AND CREATE A DIRECTORY ACCORDINGLY.
                    ModType <- substr(GraphQuant$Name, 15, 18);
                    if (substr(ModType,4,4) == " "){
                            ModType <- substr(ModType,1,3);
                    }
                    DirName <- paste(ModType, "4DModel",sep="");
                    dir.create(DirName, showWarnings = FALSE);
                                        
                    # THE SUMMARY TAB
                    Tab$Summary <- gtext(label = "Summary", Summarize(GraphQuant), container = PlotTabs);

                    # THE SILHOUETTE PLOT TAB
                    SilFile <- paste(getwd(),"/",DirName,"/Silhouette.jpg", sep = "");
                    jpeg(SilFile, width = 600, height = 600, quality = 100)
                    plot(GraphQuant$PamClRange, GraphQuant$SilPlot,'o', pch=18, main="Silhouette plot", xlab="Number of cluster", ylab="Silhouette width") 
                    dev.off()
                    Tab$Silhouette <- gimage(label = "Silhouette Plot", filename = SilFile, container = PlotTabs)
                    
                   ## THE MODEL WITH ESTIMATED CLUSTERS TAB
                    DirName_Clust <- paste(DirName, "/EstimatedClusters", sep="");
                    dir.create(DirName_Clust, showWarnings = FALSE);           
                    Colors <- rainbow(GraphQuant$OptimClust);
                    
                    Tab$ClustModel <- ggroup(tab.pos = 3, label = "Estimated Clusters", container = PlotTabs, horizontal = FALSE);
                    
                    ## Tab$ClustModel TAB CONTAINS 4 TABS, EACH WITH PLOTS CONSTRUCTED BY REMOVING ONE DIMENSION AT A TIME.
                    Tab$ClustMod1 <- gbutton("Show the model - XYZ Coordinates", container = Tab$ClustModel);
                    addSpace(Tab$ClustModel, 20, horizontal = FALSE);
                    Tab$ClustMod2 <- gbutton("Show the model - YZW Coordinates", container = Tab$ClustModel);
                    addSpace(Tab$ClustModel, 20, horizontal = FALSE);
                    Tab$ClustMod3 <- gbutton("Show the model - XZW Coordinates", container = Tab$ClustModel);
                    addSpace(Tab$ClustModel, 20, horizontal = FALSE);
                    Tab$ClustMod4 <- gbutton("Show the model - XYW Coordinates", container = Tab$ClustModel);
                    
                    open3d(useNULL = TRUE) #Estimated cluster MDS plot showing XYZ axis
                    plot3d(GraphQuant$Coords[,-4], col = Colors[GraphQuant$ClusMem], type='s',size=1,xlab="X",ylab="Y",zlab="Z",main="Cluster model");
                    writeWebGL(dir = DirName_Clust, filename = file.path(DirName_Clust, "Model_XYZ.html"), width=800, height=800, snapshot = FALSE)
                    rgl.close();
                    addHandlerClicked(Tab$ClustMod1, handler = function(h,...){ browseURL(paste(getwd(),"/",DirName_Clust,"/Model_XYZ.html",sep="")) });
        
                    open3d(useNULL = TRUE) #Estimated cluster MDS plot showing YZW axis
                    plot3d(GraphQuant$Coords[,-1], col = Colors[GraphQuant$ClusMem], type='s',size=1,xlab="Y",ylab="Z",zlab="W",main="Cluster model");
                    writeWebGL(dir = DirName_Clust, filename = file.path(DirName_Clust, "Model_YZW.html"), width=800, height=800, snapshot = FALSE)
                    rgl.close();
                    addHandlerClicked(Tab$ClustMod2, handler = function(h,...){ browseURL(paste(getwd(),"/",DirName_Clust,"/Model_YZW.html",sep="")) }); 
                    
                    open3d(useNULL = TRUE) #Estimated cluster MDS plot showing XZW axis
                    plot3d(GraphQuant$Coords[,-2], col = Colors[GraphQuant$ClusMem], type='s',size=1,xlab="X",ylab="Z",zlab="W",main="Cluster model");
                    writeWebGL(dir = DirName_Clust, filename = file.path(DirName_Clust, "Model_XZW.html"), width=800, height=800, snapshot = FALSE)
                    rgl.close();
                     addHandlerClicked(Tab$ClustMod3, handler = function(h,...){ browseURL(paste(getwd(),"/",DirName_Clust,"/Model_XZW.html",sep="")) });
                    
                    open3d(useNULL = TRUE) #Estimated cluster MDS plot showing XYW axis
                    plot3d(GraphQuant$Coords[,-3], col = Colors[GraphQuant$ClusMem], type='s',size=1,xlab="X",ylab="Y",zlab="W",main="Cluster model");
                    writeWebGL(dir = DirName_Clust, filename = file.path(DirName_Clust, "Model_XYW.html"), width=800, height=800, snapshot = FALSE)
                    rgl.close();
                    addHandlerClicked(Tab$ClustMod4, handler = function(h,...){ browseURL(paste(getwd(),"/",DirName_Clust,"/Model_XYW.html",sep="")) });
                    
                    ## SAVE THE ESTIMATED CLUSTERS AS A TEXT FILE.
                    SaveClus <- matrix(as.vector(GraphQuant$ClusMem), nrow = length(GraphQuant$ClusMem), ncol = 1);
                    rownames(SaveClus) <- paste("Sample",  as.character(seq(1, nrow(SaveClus))), sep="");
                    write.table(GraphQuant$ClusMem, file = paste(DirName,"/ClusterMembership.txt", sep = ""), sep = ",", row.names = TRUE);
                    
                    # THE MODEL WITH THE TRUE CLUSTERS TAB
                    if (!is.null(GraphQuant$TrueMem)){
                                DirName_True <- paste(DirName, "/TrueModel", sep="");
                                dir.create(DirName_True, showWarnings = FALSE);           
                                NClust <- length(table(GraphQuant$TrueMem));
                                Colors <- rainbow(NClust);
                            
                                Tab$TrueModel <- ggroup(tab.pos = 3, label = "True Clusters", container = PlotTabs, horizontal = FALSE);
                            
                                ## Tab$TrueModel TAB CONTAINS 4 TABS, EACH WITH PLOTS CONSTRUCTED BY REMOVING ONE DIMENSION AT A TIME.
                                Tab$TrueMod1 <- gbutton("Show the model - XYZ Coordinates", container = Tab$TrueModel);
                                addSpace(Tab$TrueModel, 20, horizontal = FALSE);
                                Tab$TrueMod2 <- gbutton("Show the model - YZW Coordinates", container = Tab$TrueModel);
                                addSpace(Tab$TrueModel, 20, horizontal = FALSE);
                                Tab$TrueMod3 <- gbutton("Show the model - XZW Coordinates", container = Tab$TrueModel);
                                addSpace(Tab$TrueModel, 20, horizontal = FALSE);
                                Tab$TrueMod4 <- gbutton("Show the model - XYW Coordinates", container = Tab$TrueModel);
                            
                                open3d(useNULL = TRUE) #Estimated cluster MDS plot showing XYZ axis
                                plot3d(GraphQuant$Coords[,-4], col = Colors[GraphQuant$TrueMem], type='s',size=1,xlab="X",ylab="Y",zlab="Z",main="True Clusters");
                                writeWebGL(dir = DirName_True, filename = file.path(DirName_True, "Model_XYZ.html"), width=800, height=800, snapshot = FALSE)
                                rgl.close();
                                addHandlerClicked(Tab$TrueMod1, handler = function(h,...){ browseURL(paste(getwd(),"/",DirName_True,"/Model_XYZ.html",sep="")) });
                    
                                open3d(useNULL = TRUE) #Estimated cluster MDS plot showing YZW axis
                                plot3d(GraphQuant$Coords[,-1], col = Colors[GraphQuant$TrueMem], type='s',size=1,xlab="Y",ylab="Z",zlab="W",main="True Clusters");
                                writeWebGL(dir = DirName_True, filename = file.path(DirName_True, "Model_YZW.html"), width=800, height=800, snapshot = FALSE)
                                rgl.close();
                                addHandlerClicked(Tab$TrueMod2, handler = function(h,...){ browseURL(paste(getwd(),"/",DirName_True,"/Model_YZW.html",sep="")) }); 
                                
                                open3d(useNULL = TRUE) #Estimated cluster MDS plot showing XZW axis
                                plot3d(GraphQuant$Coords[,-2], col = Colors[GraphQuant$TrueMem], type='s',size=1,xlab="X",ylab="Z",zlab="W",main="True Clusters");
                                writeWebGL(dir = DirName_True, filename = file.path(DirName_True, "Model_XZW.html"), width=800, height=800, snapshot = FALSE)
                                rgl.close();
                                addHandlerClicked(Tab$TrueMod3, handler = function(h,...){ browseURL(paste(getwd(),"/",DirName_True,"/Model_XZW.html",sep="")) });
                                
                                open3d(useNULL = TRUE) #Estimated cluster MDS plot showing XYW axis
                                plot3d(GraphQuant$Coords[,-3], col = Colors[GraphQuant$TrueMem], type='s',size=1,xlab="X",ylab="Y",zlab="W",main="True Clusters");
                                writeWebGL(dir = DirName_True, filename = file.path(DirName_True, "Model_XYW.html"), width=800, height=800, snapshot = FALSE)
                                rgl.close();
                                addHandlerClicked(Tab$TrueMod4, handler = function(h,...){ browseURL(paste(getwd(),"/",DirName_True,"/Model_XYW.html",sep="")) });
                           }
                    svalue(PlotTabs) <- 1;         
                  }   
          return(PlotTabs);
}