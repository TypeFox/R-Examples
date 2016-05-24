MakePlots3D <- function(GraphQuant){
           PlotTabs <- NULL;
           if (!is.null(GraphQuant)){
                    ## MAKE ALL THE NECESSARY PLOTS AND PRESENT THEM IN A NEW WINDOW
                    ## THE WINDOW WILL HAVE FOUR(WITH AN OPTIONAL FIFTH) WINDOWS IN THE FOLLOWING ORDER:
                    ## SUMMARY, SILHOUETTE PLOT, TRUE MODEL(BLACK), ESTIMATED CLUSTER MODEL AND TRUE CLUSTER MODEL(OPTIONAL)
                    
                    ## FOR 3D MODELS, THE ESTIMATED MODEL, ESTIMATED CLUSTER MODEL AND TRUE CLUSTER MODEL TABS WILL CONTAIN 
                    ## HTML LINKS WHICH WHEN CLICKED WILL BE OPENED IN A WEB BROWSER.

                    PlotTabs <- gnotebook(tab.pos = 3, closebuttons = FALSE);
                    Tab <- list();
                    
                    ## DETERMINE IF IT IS A PCoA MODEL OR A MDS MODEL AND CREATE A DIRECTORY ACCORDINGLY.
                    ModType <- substr(GraphQuant$Name, 15, 18);
                     if (substr(ModType,4,4) == " "){
                            ModType <- substr(ModType,1,3);
                    }
                    DirName <- paste(ModType, "3DModel",sep="");
                    dir.create(DirName, showWarnings = FALSE);
                                        
                    # THE SUMMARY TAB
                    Tab$Summary <- gtext(label = "Summary", Summarize(GraphQuant), container = PlotTabs);

                    # THE SILHOUETTE PLOT TAB
                    SilFile <- paste(getwd(),"/",DirName,"/Silhouette.jpg", sep = "");
                    jpeg(SilFile, width = 600, height = 600, quality = 100)
                    plot(GraphQuant$PamClRange, GraphQuant$SilPlot,'o', pch=18, main="Silhouette plot", xlab="Number of cluster", ylab="Silhouette width") 
                    dev.off()
                    Tab$Silhouette <- gimage(label = "Silhouette Plot", filename = SilFile, container = PlotTabs)

                   # THE MODEL WITH ESTIMATED CLUSTERS TAB
                     DirName_Clust <- paste(DirName, "/EstimatedClusters", sep="");
                    dir.create(DirName_Clust, showWarnings = FALSE);
                    open3d(useNULL = TRUE) #Estimated MDS plot
                    Colors <- rainbow(GraphQuant$OptimClust);
                    plot3d(GraphQuant$Coords, type='s', col = Colors[GraphQuant$ClusMem], size=1,xlab="X",ylab="Y",zlab="Z",main="Estimated Clusters");
                    writeWebGL(dir = DirName_Clust, filename = file.path(DirName_Clust, "Model.html"), width=800, height=800, snapshot = FALSE)
                    rgl.close();
                    Tab$Clust <- ggroup(horizontal  = FALSE, container = PlotTabs, label = "Model showing estimated clusters");
                    addSpace(Tab$Clust, 20, horizontal = FALSE);
                    Tab$ClustWeb <- gbutton("Show the 3D Model", container = Tab$Clust);
                    addHandlerClicked(Tab$ClustWeb, handler = function(h,...){ browseURL(paste(getwd(),"/",DirName_Clust,"/Model.html",sep="")) });
                    
                    ## SAVE THE ESTIMATED CLUSTERS AS A TEXT FILE.
                    SaveClus <- matrix(as.vector(GraphQuant$ClusMem), nrow = length(GraphQuant$ClusMem), ncol = 1);
                    rownames(SaveClus) <- paste("Sample",  as.character(seq(1, nrow(SaveClus))), sep="");
                    write.table(GraphQuant$ClusMem, file = paste(DirName,"/ClusterMembership.txt", sep = ""), sep = ",", row.names = TRUE);
                    
                    # THE MODEL WITH THE TRUE CLUSTERS TAB
                    if (!is.null(GraphQuant$TrueMem)){
                                DirName_True <- paste(DirName, "/TrueModel", sep="");
                                dir.create(DirName_True, showWarnings = FALSE);
                                open3d(useNULL = TRUE) #Estimated MDS plot
                                NClust <- length(table(GraphQuant$TrueMem));
                                Colors <- rainbow(NClust);
                                plot3d(GraphQuant$Coords, type='s', col = Colors[GraphQuant$TrueMem], size=1,xlab="X",ylab="Y",zlab="Z",main="Estimated metric MDS model");
                                writeWebGL(dir = DirName_True, filename = file.path(DirName_True, "Model.html"), width=800, height=800, snapshot = FALSE)
                                rgl.close();
                                Tab$True <- ggroup(horizontal  = FALSE, container = PlotTabs, label = "True clusters");
                                addSpace(Tab$True, 20, horizontal = FALSE);
                                Tab$TrueWeb <- gbutton("Show the 3D Model", container = Tab$True);
                                addHandlerClicked(Tab$TrueWeb, handler = function(h,...){ browseURL(paste(getwd(),"/",DirName_True,"/Model.html",sep="")) });
                    }

                    svalue(PlotTabs) <- 1;         
                    }          
          return(PlotTabs);
}