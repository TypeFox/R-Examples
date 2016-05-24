GraphMetagen <-
function(MDSdata){
           ## DETERMINE THE SAMPLE SIZE AND NUMBER OF FEATURES.
           N <- nrow(MDSdata$Contents);
           Features <- ncol(MDSdata$Contents);
           
           ## OUTPUT IS THE LIST CALLED GRAPHQUANTS. IT CONTAINS ALL THE MODELS(POINTS) AND MISCLASSIFICATION ERRORS.
           GraphQuants <- list();
           GraphQuants$mds2d <- NULL;
           GraphQuants$mds3d <- NULL;
           GraphQuants$mds4d <- NULL;
           GraphQuants$pcoa2d <- NULL;
           GraphQuants$pcoa3d <- NULL;
           GraphQuants$pcoa4d <- NULL;
           
           ##DEFINE GLOBAL PARAMETERS
           PamCutoff <- 2*floor(N^(1/2));
           PamClRange <- seq(MDSdata$MinClust, PamCutoff);
           PCRLength <- length(PamClRange);
           Eps <- 0.1;
           MCE <- matrix(0,6,3);
           rownames(MCE) <- c("2D MDS model", "3D MDS model", "4D MDS model", "2D PCoA model", "3D PCoA model", "4D PCoA model");
           colnames(MCE) <- c("Error Proportion", "MCE", "Optimal Clusters");
           
           ## CONVERT THE GIVEN DATA TO DISTANCE MATRIX(IF REQUIRED)
           if (MDSdata$DataType == "Distance"){
                        ## DISTANCE ALREADY CALCULATED.
                        X <- MDSdata$Contents; 
            }
            if (MDSdata$DataType == "Counts"){
                        ## CALCULATE THE DISTANCE USING THE DISTANCE TYPE SPECIFIED 
                        if (MDSdata$DistType == "Kendall's tau-distance"){
                                X <- Count2Distance(Data = MDSdata$Contents, Distance = MDSdata$DistType, Penalty = MDSdata$Penalty);
                        }
                        ## CALCULATE THE UNIFRAC USING THE OPTIONS SPECIFIED
                        if (MDSdata$DistType == "UniFrac"){
				## ADD TIP LABELS OF THE TREE AS COLUMN NAMES.
                                Counts <- MDSdata$Contents; 
                                colnames(Counts) <- MDSdata$PhyTree$Tree$tip.label;
                                ## CALCULATE THE UNIFRAC DISTANCE USING COUNT2DISTANCE
                                X <- Count2Distance(Data = Counts, Distance = MDSdata$DistType, PhyTree = MDSdata$PhyTree$Tree, UnifOpts = c(MDSdata$PhyTree$Weight, MDSdata$PhyTree$Type));
			}
            }
            KDMatrix <- X;
    #### CONSTRUCTION OF GRAPHS - EACH OPTION CREATES A NEW WINDOW.
       if (any(MDSdata$Dimensions == 2) && any(MDSdata$Norms == 1)){
                  ################################################ 2D metric MDS MODEL
                  Y1 <- isoMDS(KDMatrix, k=2, p=1, tol=1e-6, maxit=1000,trace=FALSE);
                  PamWidth1 = rep(0, PCRLength);
                  for (i in 1:PCRLength){
                        PamWidth1[i] = pam(Y1$points, k= PamClRange[i], metric="manhattan")$silinfo$avg.width;
                  }
                  OptimClust <- PamClRange[OptimClusts(PamWidth1, Eps)];
                  P = pam(Y1$points, k=OptimClust, metric="manhattan");
                  D1 <- MatrixkNorm(Y1$points,1);
                  MCE[1,] <- c(sum((KDMatrix - D1)^2)/sum(D1^2), MCError(MDSdata$Clust,P$clustering),OptimClust);
                  GraphQuants$mds2d$Name <- "2 Dimensional MDS model";
                  GraphQuants$mds2d$Coords <- Y1$points;
                  GraphQuants$mds2d$ClusMem <- P$clustering;
                  GraphQuants$mds2d$TrueMem <- MDSdata$Clust;
                  GraphQuants$mds2d$PamClRange <- PamClRange;
                  GraphQuants$mds2d$OptimClust <- OptimClust; 
                  GraphQuants$mds2d$SilPlot <- PamWidth1;
       }
       ################################################ 3D metric MDS MODEL
       if (any(MDSdata$Dimensions == 3) && any(MDSdata$Norms == 1)){
                  Y2 <- isoMDS(KDMatrix,k=3,p=1,tol=1e-6,maxit=1000,trace=FALSE);
                  PamWidth2 = rep(0, PCRLength);
                  for (i in 1:PCRLength){
                        PamWidth2[i] = pam(Y2$points, k = PamClRange[i], metric = "manhattan")$silinfo$avg.width;
                  }
                  OptimClust <- PamClRange[OptimClusts(PamWidth2, Eps)];
                  P = pam(Y2$points, k = OptimClust, metric = "manhattan");
                  D2 <- MatrixkNorm(Y2$points,1);
                  MCE[2,] <- c(sum((KDMatrix - D2)^2)/sum(D2^2), MCError(MDSdata$Clust,P$clustering),OptimClust);
                  GraphQuants$mds3d$Name <- "3 Dimensional MDS model";
                  GraphQuants$mds3d$Coords <- Y2$points;
                  GraphQuants$mds3d$ClusMem <- P$clustering;
                  GraphQuants$mds3d$TrueMem <- MDSdata$Clust;
                  GraphQuants$mds3d$PamClRange <- PamClRange;
                  GraphQuants$mds3d$OptimClust <- OptimClust; 
                  GraphQuants$mds3d$SilPlot <- PamWidth2;
       }
       ###############################################4D metric MDS MODEL
       if (any(MDSdata$Dimensions == 4) && any(MDSdata$Norms == 1)){
                  Y3 = isoMDS(KDMatrix,k=4,p=1,tol=1e-6,maxit=1000,trace=FALSE);
                  PamWidth3 = rep(0, PCRLength);
                  for (i in 1:PCRLength){
                        PamWidth3[i] = pam(Y3$points, k = PamClRange[i], metric="manhattan")$silinfo$avg.width;
                  }
                  ##Plots
                  OptimClust <- PamClRange[OptimClusts(PamWidth3, Eps)];
                  P = pam(Y3$points,k=OptimClust,metric="manhattan")
                  D3 <- MatrixkNorm(Y3$points,1);
                  MCE[3,] <- c(sum((KDMatrix - D3)^2)/sum(D3^2), MCError(MDSdata$Clust,P$clustering),OptimClust);
                  GraphQuants$mds4d$Name <- "4 Dimensional MDS model";
                  GraphQuants$mds4d$Coords <- Y3$points;
                  GraphQuants$mds4d$ClusMem <- P$clustering;
                  GraphQuants$mds4d$TrueMem <- MDSdata$Clust;
                  GraphQuants$mds4d$PamClRange <- PamClRange;
                  GraphQuants$mds4d$OptimClust <- OptimClust; 
                  GraphQuants$mds4d$SilPlot <- PamWidth3;
       }
       ################################################2D PCoA MODEL
       if (any(MDSdata$Dimensions == 2) && any(MDSdata$Norms == 2)){
                  Y4 <- cmdscale(KDMatrix,k=2);
                  PamWidth4 = rep(0, PCRLength);
                  for (i in 1:PCRLength){
                        PamWidth4[i-1] = pam(Y4, k = PamClRange[i], metric="manhattan")$silinfo$avg.width;
                  }
                  OptimClust <- PamClRange[OptimClusts(PamWidth4, Eps)];
                  P = pam(Y4,k=OptimClust,metric="manhattan");
                  D4 <- MatrixkNorm(Y4,1);
                  MCE[4,] <- c(sum((KDMatrix - D4)^2)/sum(D4^2), MCError(MDSdata$Clust,P$clustering),OptimClust);
                  GraphQuants$pcoa2d$Name <- "2 Dimensional PCoA model";
                  GraphQuants$pcoa2d$Coords <- Y4;
                  GraphQuants$pcoa2d$ClusMem <- P$clustering;
                  GraphQuants$pcoa2d$TrueMem <- MDSdata$Clust;
                  GraphQuants$pcoa2d$PamClRange <- PamClRange;
                  GraphQuants$pcoa2d$OptimClust <- OptimClust; 
                  GraphQuants$pcoa2d$SilPlot <- PamWidth4;
        }          
       ###############################################3D PCoA MODEL
       if (any(MDSdata$Dimensions == 3) && any(MDSdata$Norms == 2)){
                  Y5 <- cmdscale(KDMatrix,k=3);
                  PamWidth5 = rep(0, PCRLength);
                  for (i in 1:PCRLength){
                        PamWidth5[i-1] = pam(Y5, k = PamClRange[i], metric="manhattan")$silinfo$avg.width;
                  }
                  OptimClust <- PamClRange[OptimClusts(PamWidth5, Eps)];
                  P = pam(Y5,k=OptimClust,metric="manhattan");
                  D5 <- MatrixkNorm(Y5,1);
                  MCE[5,] <-c(sum((KDMatrix - D5)^2)/sum(D5^2), MCError(MDSdata$Clust,P$clustering),OptimClust);
                  GraphQuants$pcoa3d$Name <- "3 Dimensional PCoA model";
                  GraphQuants$pcoa3d$Coords <- Y5;
                  GraphQuants$pcoa3d$ClusMem <- P$clustering;
                  GraphQuants$pcoa3d$TrueMem <- MDSdata$Clust;
                  GraphQuants$pcoa3d$PamClRange <- PamClRange;
                  GraphQuants$pcoa3d$OptimClust <- OptimClust; 
                  GraphQuants$pcoa3d$SilPlot <- PamWidth5;
       }
       ###############################################4D PCoA MODEL
       if (any(MDSdata$Dimensions == 4) && any(MDSdata$Norms == 2)){
                  Y6 = cmdscale(KDMatrix,k=4);
                  PamWidth6 = rep(0, PCRLength);
                  for (i in 1:PCRLength){
                        PamWidth6[i-1] = pam(Y6, k = PamClRange[i], metric="manhattan")$silinfo$avg.width;
                  }
                  OptimClust <- PamClRange[OptimClusts(PamWidth6, Eps)];
                  P = pam(Y6,k=OptimClust,metric="manhattan")
                  D6 <- MatrixkNorm(Y6,1);
                  MCE[6,] <- c(sum((KDMatrix - D6)^2)/sum(D6^2), MCError(MDSdata$Clust,P$clustering),OptimClust);
                  GraphQuants$pcoa4d$Name <- "4 Dimensional PCoA model";
                  GraphQuants$pcoa4d$Coords <- Y6;
                  GraphQuants$pcoa4d$ClusMem <- P$clustering;
                  GraphQuants$pcoa4d$TrueMem <- MDSdata$Clust;
                  GraphQuants$pcoa4d$PamClRange <- PamClRange;
                  GraphQuants$pcoa4d$OptimClust <- OptimClust; 
                  GraphQuants$pcoa4d$SilPlot <- PamWidth6;
       }
       GraphQuants$MCE <- MCE;
       return(GraphQuants);
}
