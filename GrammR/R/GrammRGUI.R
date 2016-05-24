GrammRGUI <- function(Direc = getwd()){
## CHANGE DIRECTORY TO THE DIRECTORY SPECIFIED. IF NOT SPECIFIED, USE THE CURRENT DIRECTORY.
setwd(Direc)
## CREATE THE WINDOW TO COLLECT ALL INFORMATION NEEDED TO PROCESS THE MODEL.
DataWindow <- gwindow(title="MDS representations",width=800,height=600);
MainGroup <- ggroup(horizontal = FALSE, container = DataWindow);

DATA <- glabel(text = "Data Selection", container = MainGroup);
font(DATA) = c(weight = "bold", size = 14);
addSpace(MainGroup, 10, horizontal = FALSE);

## LOAD DATA FILE.
glabel(text = "Select the data file containing the metagenomic count data", container = MainGroup);
DataBrowse = gfilebrowse( text = "Select a file", type="open", container=MainGroup);
## OPTIONS WHICH SELECT THE FORMAT OF THE DATA FILE.
addSpace(MainGroup, 10, horizontal = FALSE)
DataFormat <- ggroup(horizontal = TRUE, container = MainGroup);
Delimiter.Label1 <- glabel("Select the delimiter in the data file", container = DataFormat);
addSpring(DataFormat);
Delimiter.button1 <- gdroplist(items = c("Space","Tab","Comma"), selected = 1, container = DataFormat);

addSpace(MainGroup, 20, horizontal = FALSE);

## CHOOSE WHAT TO CLUSTER - ROWS OR COLUMNS.
ClustWhat <- ggroup(horizontal = TRUE, container = MainGroup);
ClusWhat1 <- glabel("Choose what is to be clustered", container = ClustWhat);
ClusWhat2 <- gimage(filename = "info", dirname = "stock", container = ClustWhat);
addSpace(ClustWhat, 5);
addSpring(ClustWhat, horizontal = TRUE);
ClusWhat3 <- gradio(items =c("Rows", "Columns"), container = ClustWhat, horizontal = TRUE);

addHandlerClicked(ClusWhat2, handler = function(h,...){
gmessage("Choose the quantity that you would like to cluster. \n
If the data file contains counts with samples as rows and taxa as cluster, select rows to cluster the samples and columns to cluster the taxa.
If the data file contains distances, the choice is redundant.")
});

## LOAD THE CLUSTER FILE
addSpace(MainGroup, 20, horizontal = FALSE)
CLUST <- glabel(text="Select the data file containing pre-specified cluster membership vector(Optional) ", container = MainGroup);
ClustBrowse = gfilebrowse(text = "Select a file", type = "open", container = MainGroup);

## OPTIONS WHICH SELECT THE FORMAT OF THE DATA FILE.
addSpace(MainGroup, 10, horizontal = FALSE)
ClustFormat <- ggroup(horizontal = TRUE, container = MainGroup);
Delimiter.Label2 <- glabel("Select the delimiter in the cluster file", container = ClustFormat);
addSpring(ClustFormat);
Delimiter.button2 <- gdroplist(items = c("Space","Tab","Comma"), selected = 1, container = ClustFormat);



## ENTER THE MODEL PARAMETERS
ModelFormat <- ggroup(horizontal = FALSE, container = MainGroup);
addSpace(ModelFormat,10,horizontal = FALSE);
ModelParam <- glabel(text="Modeling parameters");
font(ModelParam) = c(weight="bold", size = 14);
add(ModelFormat,ModelParam);
addSpace(ModelFormat,10,horizontal = FALSE);

## DETERMINE THE CONTENTS OF THE DATA FILE ENTERED.
Contents <- ggroup(horizontal = TRUE, container = ModelFormat);
glabel(text = "Select the contents of the data matrix to be analyzed", container = Contents);
addSpring(Contents);
Cont <- gdroplist(items = c("\t\t\t\t","Counts", "Distance"), selected = 0,  container = Contents);

## SELECT THE DISTANCE MEASURE
Distance <- ggroup(horizontal = TRUE, visible = FALSE);
glabel(text = "Select the distance measure", container = Distance);
addSpring(Distance);
Dist <- gdroplist(items = c("\t\t\t\t","Kendall's tau-distance", "UniFrac"), selected = 0, container = Distance);    

## SPECIFY THE KENDALL'S TAU PENALTY
Penalty <- ggroup(horizontal = TRUE, visible = FALSE);
glabel(text = "Specify the penalty for ties", container = Penalty);
addSpring(Penalty);
PenSlider <- gslider(from = 0, to = 1, by = 0.01, value = 0.5, horizontal = TRUE, container = Penalty, expand = TRUE);
addSpace(Penalty, 3, horizontal = TRUE);

## SUPPLY THE UNIFRAC TREE FILE        
Unif <- ggroup(horizontal = FALSE, visible = FALSE);
glabel("Select the file containing the phylogenetic tree for calculating UniFrac", container = Unif);
PhyTreeBrowse <- gfilebrowse(text = "Phylogenetic Tree", type="open", container = Unif);

## OPTIONS FOR GUniFrac DISTANCE CALCULATION.
addSpace(Unif, 10);
UnifType <- ggroup(horizontal = TRUE, container = Unif);
glabel("Select the type of UniFrac distance to be calculated", container = UnifType);
addSpring(UnifType, horizontal = TRUE);
GunifType <- gdroplist(items = c("Unweighted", "Variance Adjusted", "Generalized"), selected=1, container = UnifType);
addSpace(Unif, 10);
UnifWt <- ggroup(horizontal = TRUE, container = Unif);
glabel("Specify the weight to calculate Generalized UniFrac", container = UnifWt);
GunifWeight <- gslider(from = 0, to = 1, by = 0.01, value = 0.5, horizontal = TRUE, container = UnifWt, expand = TRUE);

addHandlerChanged(Cont, handler = function(h,...){ 
                                delete(ModelFormat, Distance);
                                delete(ModelFormat, Penalty);
                                delete(ModelFormat, Unif);
                                if(svalue(Cont) == "Counts"){
                                        add(ModelFormat, Distance);
                                }
});
addHandlerChanged(Dist, handler = function(h,...){
                                delete(ModelFormat, Penalty);
                                delete(ModelFormat, Unif);
                                if(svalue(Dist) == "Kendall's tau-distance"){
                                        delete(ModelFormat, Unif);
                                        add(ModelFormat, Penalty);
                                        PhyTreeYN <- "None";
                                 }else if(svalue(Dist) == "UniFrac"){
                                        delete(ModelFormat, Penalty);
                                        add(ModelFormat, Unif);
                                        PhyTreeYN <- "Yes";
                                 }else{
                                        delete(ModelFormat, Penalty);
                                        delete(ModelFormat, Unif);
                                }
});

## SELECT THE DIMENSIONS TO BE USED
addSpace(ModelFormat, 10, horizontal = FALSE);
ModelSpecs <- ggroup(horizontal = FALSE, container = MainGroup);
Dimension <- ggroup(horizontal = TRUE, container = ModelSpecs);
glabel(text = "Dimensions for which models are to be constructed", container = Dimension);
addSpring(Dimension);
Dim <- gcheckboxgroup(items = c("2","3","4"), horizontal = TRUE, container = Dimension);
## NORM TO BE USED FOR THE MDS MODELS
Norm <- ggroup(horizontal = TRUE, container = ModelSpecs);
glabel(text="Norm to be used for the MDS models", container = Norm);
addSpring(Norm);
LpNorm <-  gcheckboxgroup(items =c("1","2"), horizontal = TRUE, container = Norm);
## SPECIFY THE MINIMUM NUMBER OF CLUSTERS TO USE. DEFAULT IS 2
ClustSpec <- ggroup(horizontal = TRUE, container = ModelSpecs)
glabel(text = "(Optional) Specify the minimum number of clusters(Default value is 2)", container = ClustSpec);
addSpring(ClustSpec);
MinClust <- gtext(text = "2", container = ClustSpec, width = 50, height = 27);
#font(MinClust) <- c(weights = "bold");

## LOAD ALL DATA FILES AND CONSTRUCTING PLOTS USING ALL THE INFORMATION PROVIDED

FileRead <- function(h,...){
         ## READ THE DATA FILE USING THE SPECIFIED OPTIONS
         
         ## CREATE A NEW WINDOW SHOWING THE PROGRESS FOR EACH OF THE STEPS
         ProgWind <- gwindow(title = "Progress");
         ProgWindow <- ggroup(horizontal = FALSE, container = ProgWind);
         ReadProgGroup <- ggroup(container = ProgWindow, horizontal = TRUE);
         ReadProg <- glabel("Reading the data files");
         add(ReadProgGroup, ReadProg);
         
        Data.FileName <- svalue(DataBrowse);
        Data.FileName <- substr(Data.FileName, 2, nchar(Data.FileName)-1);
        
        ## DEFINE THE FOLDER NAME WHERE ALL THE SUBDIRECTORIES 
        FindPeriod <- gregexpr("\\.", Data.FileName)[[1]];
        Data.FolderName <- substr(Data.FileName, 1, FindPeriod[length(FindPeriod)]-1);
        Data.FolderName <- paste(Data.FolderName, svalue(Dist),sep="_");
        dir.create(Data.FolderName, showWarnings = FALSE);
        setwd(Data.FolderName);
        
        if (svalue(Delimiter.button1) == "Space"){
                 X <- as.matrix(read.table(Data.FileName));
         }
        if (svalue(Delimiter.button1) == "Tab"){
                X <- as.matrix(read.table(Data.FileName, sep = "\t"));              
         }
         if (svalue(Delimiter.button1) == "Comma"){
                X <- as.matrix(read.table(Data.FileName, sep = ","));
         }
         ## READ THE CLUSTER FILE USING THE SPECIFIED OPTIONS
         Clust.FileName <- svalue(ClustBrowse);
         Y <- NULL;
         if ( file.exists(Clust.FileName)){
                    Clust.FileName <- substr(Clust.FileName, 2, nchar(Clust.FileName) - 1);
                    if (svalue(Delimiter.button2) == "Space"){
                            Y <- unlist(read.table(Clust.FileName));
                    }
                    if (svalue(Delimiter.button2) == "Tab"){
                            Y <- unlist(read.table(Clust.FileName, sep = "\t"));
                    }
                    if (svalue(Delimiter.button2) == "Comma"){
                            Y <- unlist(read.table(Clust.FileName, sep = ","));
                    }
         }
                    
          ## READ THE PHYLOGENETIC TREE FROM THE FILE SELECTED
          PhyTree <- NULL;
          if (file.exists(svalue(PhyTreeBrowse))){
                    PhyTree.FileName <- svalue(PhyTreeBrowse);
                    PhyTree.FileName <- substr(PhyTree.FileName, 2, nchar(PhyTree.FileName)-1);
                    PhyTree <- read.tree(PhyTree.FileName);
          }
          
          ## MDSdata IS THE LIST WE PROVIDE AS THE ARGUMENT FOR MDS CONSTRUCTION AND PLOTTING
          MDSdata <- list();
          if (svalue(ClusWhat3) == "Rows"){
                    MDSdata$Contents <- X;
          }
          if (svalue(ClusWhat3) == "Columns"){
                    MDSdata$Contents <- t(X);
          }
          MDSdata$Clust <- Y;
          MDSdata$DataType <- svalue(Cont);
          MDSdata$DistType <- svalue(Dist);
          MDSdata$PhyTree$Tree <- PhyTree;
          MDSdata$PhyTree$Type <- svalue(GunifType);
          MDSdata$PhyTree$Weight <- as.numeric(svalue(GunifWeight));
          
          ## COLLECT THE INFORMATION ENTERED FOR MODELLING.
          MDSdata$Dimensions <- svalue(Dim);
          MDSdata$Norms <- svalue(LpNorm);
          MDSdata$Penalty <- svalue(PenSlider);
          MDSdata$MinClust <- as.numeric(svalue(MinClust));
          
          ## PROGRESS UPDATE
          delete(ReadProgGroup, ReadProg);
          ReadProg <- glabel("Reading the data file...... DONE \n");
          add(ReadProgGroup, ReadProg);
          
          ## CALCULATE THE QUANTITY NEEDED FOR PLOTTING THE GRAPHICAL MODEL.
          ModelProgGroup <- ggroup(container = ProgWindow, horizontal = TRUE);
          ModelProg <- glabel("Constructing the graphical models")
          add(ModelProgGroup, ModelProg);
          GraphQuants = GraphMetagen(MDSdata);
          delete(ModelProgGroup, ModelProg);
          ModelProg <- glabel("Constructing the graphical models...... DONE \n")
          add(ModelProgGroup, ModelProg);
          
          ## GRAPHQUANTS RETURNS ALL THE INFORMATION WHICH IS CALCULATED USING GRAPHMETAGEN FUNCTION.
          ## THESE QUANTITIES WILL NOW BE USED TO CONSTRUCT GRAPHICAL REPRESENTATIONS.
          
          FinalWindow <- gwindow(title="GrammR - Results",width=600,height=600, visible = FALSE);
          FWGroups <-  ggroup(container = FinalWindow, horizontal = TRUE);
          
          
          ## CREATE BUTTONS TO BE SHOWN ON THE LEFT HAND SIDE
          FWTabs <- ggroup(container = FWGroups, horizontal = FALSE);
          FWButt1 <- gbutton(text = "2D PCoA Model")
          FWButt2 <- gbutton(text = "2D MDS Model")
          
          FWButt3 <- gbutton(text = "3D PCoA Model")
          FWButt4 <- gbutton(text = "3D MDS Model")
          
          FWButt5 <- gbutton(text = "4D PCoA Model")
          FWButt6 <- gbutton(text = "4D MDS Model")
          
          ## INCLUDE ONLY THOSE BUTTONS CORRESPONDING TO THE MODELS SELECTED
          if (2 %in% MDSdata$Dimensions && 2 %in% MDSdata$Norms){
                        add(FWTabs, FWButt1);}
          if (2 %in% MDSdata$Dimensions && 1 %in% MDSdata$Norms){
                        add(FWTabs, FWButt2); }                
          if (3 %in% MDSdata$Dimensions && 2 %in% MDSdata$Norms){
                        add(FWTabs, FWButt3);}
          if (3 %in% MDSdata$Dimensions && 1 %in% MDSdata$Norms){
                        add(FWTabs, FWButt4);}
          if (4 %in% MDSdata$Dimensions && 2 %in% MDSdata$Norms){
                        add(FWTabs, FWButt5);}
          if (4 %in% MDSdata$Dimensions && 1 %in% MDSdata$Norms){
                        add(FWTabs, FWButt6);}
          
          ## CONSTRUCT AND DISPLAY PROGRESS FOR 2D MODELS.
          Prog2DGroup <- ggroup(container = ProgWindow, horizontal = TRUE);
          Prog2D <- glabel("Adding 2D Models to the output window");
          add(Prog2DGroup, Prog2D);
          Mod1 <- MakePlots2D(GraphQuants$pcoa2d);
          Mod2 <- MakePlots2D(GraphQuants$mds2d);
          delete(Prog2DGroup, Prog2D);
          Prog2D <- glabel("Adding 2D Models to the output window...... DONE");
          add(Prog2DGroup, Prog2D);
          
          ## CONSTRUCT AND DISPLAY PROGRESS FOR 3D MODELS.
          Prog3DGroup <- ggroup(container = ProgWindow, horizontal = TRUE);
          Prog3D <- glabel("Adding 3D Models to the output window");
          add(Prog3DGroup, Prog3D);
          Mod3 <- MakePlots3D(GraphQuants$pcoa3d)
          Mod4 <- MakePlots3D(GraphQuants$mds3d)
          delete(Prog3DGroup, Prog3D);
          Prog3D <- glabel("Adding 3D Models to the output window...... DONE");
          add(Prog3DGroup, Prog3D);
          
          ## CONSTRUCT AND DISPLAY PROGRESS FOR 3D MODELS.
          Prog4DGroup <- ggroup(container = ProgWindow, horizontal = TRUE);
          Prog4D <- glabel("Adding 4D Models to the output window");
          add(Prog4DGroup, Prog4D);
          Mod5 <- MakePlots4D(GraphQuants$pcoa4d)
          Mod6 <- MakePlots4D(GraphQuants$mds4d)
          delete(Prog4DGroup, Prog4D);
          Prog4D <- glabel("Adding 4D Models to the output window...... DONE");
          add(Prog4DGroup, Prog4D);
          
          ## ADD A BUTTON TO SHOW THE RESULTS
          addSpace(ProgWindow, 20);
          ProgFinalGroup <- ggroup(container = ProgWindow, horizontal = TRUE);
          FinalText <- paste("Analysis Complete. \n
          Plots are saved in the following folder:", getwd());
          glabel(text = FinalText, container = ProgFinalGroup);
          
          addSpring(ProgWindow);
          FinalButton <- gbutton(text = "Show Results", container = ProgWindow);
          
          addHandlerClicked(FWButt1, handler = function(h,...){ 
                   if (!is.null(Mod1)){delete(FWGroups, Mod1); }
                   if (!is.null(Mod2)){delete(FWGroups, Mod2); }
                   if (!is.null(Mod3)){delete(FWGroups, Mod3); }
                   if (!is.null(Mod4)){delete(FWGroups, Mod4); }
                   if (!is.null(Mod5)){delete(FWGroups, Mod5); }
                   if (!is.null(Mod6)){delete(FWGroups, Mod6); }
                   add(FWGroups, Mod1);
           })
          addHandlerClicked(FWButt2, handler = function(h,...){ 
                   if (!is.null(Mod1)){delete(FWGroups, Mod1); }
                   if (!is.null(Mod2)){delete(FWGroups, Mod2); }
                   if (!is.null(Mod3)){delete(FWGroups, Mod3); }
                   if (!is.null(Mod4)){delete(FWGroups, Mod4); }
                   if (!is.null(Mod5)){delete(FWGroups, Mod5); }
                   if (!is.null(Mod6)){delete(FWGroups, Mod6); }
                   add(FWGroups, Mod2);
           })
          addHandlerClicked(FWButt3, handler = function(h,...){ 
                   if (!is.null(Mod1)){delete(FWGroups, Mod1); }
                   if (!is.null(Mod2)){delete(FWGroups, Mod2); }
                   if (!is.null(Mod3)){delete(FWGroups, Mod3); }
                   if (!is.null(Mod4)){delete(FWGroups, Mod4); }
                   if (!is.null(Mod5)){delete(FWGroups, Mod5); }
                   if (!is.null(Mod6)){delete(FWGroups, Mod6); }
                   add(FWGroups, Mod3);
           })
          addHandlerClicked(FWButt4, handler = function(h,...){ 
                   if (!is.null(Mod1)){delete(FWGroups, Mod1); }
                   if (!is.null(Mod2)){delete(FWGroups, Mod2); }
                   if (!is.null(Mod3)){delete(FWGroups, Mod3); }
                   if (!is.null(Mod4)){delete(FWGroups, Mod4); }
                   if (!is.null(Mod5)){delete(FWGroups, Mod5); }
                   if (!is.null(Mod6)){delete(FWGroups, Mod6); }
                   add(FWGroups, Mod4);
           })
          addHandlerClicked(FWButt5, handler = function(h,...){ 
                   if (!is.null(Mod1)){delete(FWGroups, Mod1); }
                   if (!is.null(Mod2)){delete(FWGroups, Mod2); }
                   if (!is.null(Mod3)){delete(FWGroups, Mod3); }
                   if (!is.null(Mod4)){delete(FWGroups, Mod4); }
                   if (!is.null(Mod5)){delete(FWGroups, Mod5); }
                   if (!is.null(Mod6)){delete(FWGroups, Mod6); }
                   add(FWGroups, Mod5);
           })
         addHandlerClicked(FWButt6, handler = function(h,...){ 
                   if (!is.null(Mod1)){delete(FWGroups, Mod1); }
                   if (!is.null(Mod2)){delete(FWGroups, Mod2); }
                   if (!is.null(Mod3)){delete(FWGroups, Mod3); }
                   if (!is.null(Mod4)){delete(FWGroups, Mod4); }
                   if (!is.null(Mod5)){delete(FWGroups, Mod5); }
                   if (!is.null(Mod6)){delete(FWGroups, Mod6); }
                   add(FWGroups, Mod6);
           })
                
         addHandlerClicked(FinalButton, handler = function(h,...){visible(FinalWindow) <- TRUE;
         delete(ProgWindow, FinalButton);});
         visible(ProgWind) <- TRUE;
}
addSpring(horizontal = FALSE, MainGroup);

## Move to the next window - Presentation of the plots.
Analyze <- gaction(label="Analyze", handler= FileRead)
Next <- gbutton(text = "Analyze", action = Analyze, container = MainGroup);
}
