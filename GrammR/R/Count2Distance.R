Count2Distance <-
function(Data, Distance, Penalty = NULL, PhyTree = NULL, UnifOpts = NULL, Adjust = TRUE){
        ##Data <- The count matrix
        ##Distance <- Kendall's tau-distance or UniFrac, as specified in the file.
        ##PhyTree <- Phylogenetic tree provided by the user for calculating UniFrac.
        ##If not provided, Kendall's tau-distance is calculated by default.
        N <- nrow(Data); ## Sample size.
        Features <- ncol(Data);  ##Number of taxa/features.
        D <- matrix(0, N, N);
        MinCorr <- .Machine$double.eps;   

        ##REMOVE SAMPLES WHICH HAVE ZERO READS.
        AnyZeros <- which(rowSums(Data) == 0);
        if (length(AnyZeros) != 0){
                   Data <- Data[-AnyZeros,];
          }        
        
         ##CONVERTING COUNTS TO RELATIVE FREQUENCIES.
         DataRF <- Data/rowSums(Data);
         if (Distance == "Kendall's tau-distance"){
                  D = KenDist(Data,  F);
          }
          if (Distance == "UniFrac" && !is.null(PhyTree) && length(UnifOpts) == 2){
                    X <- GUniFrac(otu.tab = Data, tree = PhyTree, alpha = as.numeric(UnifOpts[1]))$unifracs;
                    UFOpts <-  c("Unweighted", "Variance Adjusted", "Generalized")
                    UnifOptions <- c("d_UW", "d_VAW", paste("d_",as.character(UnifOpts[1]), sep = ""));                                
                    D <- as.matrix(X[,, UnifOptions[which(UnifOpts[2] == UFOpts)]]);
          }
          
          ## ADJUSTING THE OFF-DIAGONAL ZEROS BY ADDING AN INIFINITESIMAL VALUE SO THAT ISOMDS WORKS
          if (Adjust){
                    D <- (D + MinCorr*(D == 0)) - MinCorr*diag(1,N);
          }
          return(D);
}
