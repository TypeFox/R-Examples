GrammRServ <-
function(Data = NULL, Cluster = NULL, DataType = "Counts", DistType = "Kendall's tau-distance", PhyTree = NULL, GunifType = NULL, GunifWeight = 0, Dim = c(2,3,4), LpNorm = c(1), Penalty = 0.5, MinClust = 2){
         
         ## CHECK IF THE DATA IS ENTERED PROPERLY.
        
        MDSdata <- list();
        ## THIS IS THE INFORMATION PERTAINING TO THE DATA ENTERED.
        MDSdata$Contents <- Data;
        MDSdata$Clust <- Cluster;
        MDSdata$DataType <- DataType;
        MDSdata$DistType <- DistType;
        MDSdata$PhyTree$Tree <- PhyTree;
	MDSdata$PhyTree$Type <- GunifType;
	MDSdata$PhyTree$Weight <- GunifWeight;

        ## COLLECT THE INFORMATION ENTERED FOR MODELLING.
        MDSdata$Dimensions <- Dim;
        MDSdata$Norms <- LpNorm;
        MDSdata$Penalty <- Penalty;
        MDSdata$MinClust <- MinClust;
        GraphQuants = GraphMetagen(MDSdata);

        ## MAKE DIRECTORIES AND CONSTRUCT PLOTS.
        ## DEFINE THE FOLDER NAME WHERE ALL THE SUBDIRECTORIES 
        
        Data.FolderName <- paste("Data", DistType,sep="_");
        dir.create(Data.FolderName, showWarnings = FALSE);
        setwd(Data.FolderName);
        
        ## ALL THE PLOTS WILL BE SAVED AS PDFs AND THE WEBPAGES IN SEPARATE DIRECTORIES.
        Mod1 <- Make2DPlots(GraphQuants$pcoa2d)
        Mod2 <- Make2DPlots(GraphQuants$mds2d)
        Mod3 <- Make3DPlots(GraphQuants$pcoa3d)
        Mod4 <- Make3DPlots(GraphQuants$mds3d)
        Mod5 <- Make4DPlots(GraphQuants$pcoa4d)
        Mod6 <- Make4DPlots(GraphQuants$mds4d)
        setwd("../")
        message(paste("Analysis Complete. \nGraphical models are saved in subdirectories within ",  paste(getwd(), "/",Data.FolderName,sep=""),  " \nCurrent working directory:", getwd()));
        
 }
