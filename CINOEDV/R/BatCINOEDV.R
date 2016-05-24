BatCINOEDV <- function(FileName,MaxOrder=3,RatioThreshold=c(1,1,1),
                       NumberThreshold=c(10,10,10),measure=1,Strategy=2,
                       Population=1000,Iteration=100,SNPNameFileName){
  
  # Batch mode for using CINOEDV
  #
  # input
  #    FileName: names of SNP files.
  #              for example, FileName <- c("test.mat","test1.mat")
  #    MaxOrder: The specified maximum order, must be setted as 1,2,3,4 or 5.
  #              for example, MaxOrder <- 3
  #    RatioThreshold: Ratio thresholds control the numbers of retained top SNPs 
  #                   or SNP-combinations. Each of them should be setted in [0,1],
  #                   By default, it is setted as 1. the length of RatioThreshold
  #                   should be equal to MaxOrder.
  #                   for example, RatioThreshold <- c(1,1,1)
  #    NumberThreshold: Number thresholds control the numbers of retained top SNPs 
  #                     or SNP-combinations.Each of them should be setted as an 
  #                     integer, By default, it is setted as 10. The length of 
  #                     NumberThreshold should be equal to MaxOrder.
  #                     for example, NumberThreshold <- c(10,10,10)
  #    measure: the label of current used evaluation measure
  #             1 -> The classic co-information measure
  #             2 -> The Normalized co-information measure
  #             3 -> TingHu's Co-Information Measure
  #             others -> The classic co-information measure
  #             for example, measure <- 1
  #    Strategy: The search strategy
  #             1 -> The exhaustive search strategy
  #             2 -> The PSO search strategy
  #    Population: The number of particles
  #                 for example, Population <- 1000
  #    Iteration: The iteration number
  #               for example, Iteration <- 100
  #    SNPNameFileName: names of SNP-name files.
  #                     the length of SNPNameFileName should be equal to FileName.
  #                     If not exist such SNP-name file, please input NA
  #                     for example, SNPNameFileName <- c("test_Name.mat","test1_Name.mat")
  #
  # Junliang Shang
  # 4.22/2014
  
  # library(R.matlab)
  
  FileNum <- length(FileName)
  
  if(FileNum!=length(SNPNameFileName)){
    stop(" the length of FileName is not equal to the length of SNPNameFileName !\n")
  }
  
  if((MaxOrder!=length(RatioThreshold)) || (MaxOrder!=length(NumberThreshold))){
    stop(" the lengths of RatioThreshold and NumberThreshold are not equal to MaxOrder !\n")
  }
  
  # Check MaxOrder
  TestMaxOrder(as.character(MaxOrder))
  MaxOrder <- as.numeric(MaxOrder)
  
  # Check RatioThreshold
  TestRatioThreshold(MaxOrder,as.character(RatioThreshold))
  RatioThreshold <- as.numeric(RatioThreshold)
  
  # Check NumberThreshold
  TestNumberThreshold(MaxOrder,as.character(NumberThreshold))
  NumberThreshold <- as.numeric(NumberThreshold)
  
  # Check measure
  if (as.character(measure) %in% c("1","2","3")){
    measure <- as.numeric(measure)
  }else
  {
    measure <- 1
  }
  
  # check strategy
  if (as.character(Strategy) %in% c("1","2")){
    Strategy <- as.numeric(Strategy)
  }else
  {
    Strategy <- 1
  }
  
  if (Strategy==2){
    # check PSO based parameters, that is, Population, Iteration
    TestPSOParameters(as.character(Population),as.character(Iteration))
    
    Population <- as.numeric(Population)
    Iteration <- as.numeric(Iteration)
  }
  
  for (kk in 1:FileNum){
    
    ########################
    # Check Parameters
    ########################
    
    # Check FileName
    Data <- InputData(as.character(FileName[kk]))
    pts <- Data$pts
    class <- Data$class

    # Check SNPNameFileName
    SNPNames <- TestSNPNameFile(ncol(pts),as.character(SNPNameFileName[kk]))
    SNPNames <- SNPNames$SNPNames
    
    #############################
    # Define file name
    #############################
    # Define file name which is used for saving results.
    
    RatioN <- "("
    for (i in 1:MaxOrder){
      RatioN <- paste(RatioN,RatioThreshold[i])
    }
    RatioN <- paste(RatioN,")")
    
    
    NumberN <- "("
    for (i in 1:MaxOrder){
      NumberN <- paste(NumberN,NumberThreshold[i])
    }
    NumberN <- paste(NumberN,")")
    
    if (Strategy==1) {
      SaveFileName <- paste(substr(FileName[kk],1,nchar(FileName[kk])-4),MaxOrder,
                            RatioN,NumberN,Strategy,sep="_")
    }
    if (Strategy==2) {
      SaveFileName <- paste(substr(FileName[kk],1,nchar(FileName[kk])-4),MaxOrder,
                            RatioN,NumberN,Strategy,Population,Iteration,sep="_")
    }
    
    #############################
    # Search Strategies
    #############################
    # Search SNPs with high main effects and SNP-combinations with high interaction
    # effects. In current version, only one strategy is provided, that is, 
    # exhausative search.
    
    cat("#### Search Results ####\n\n")
    cat(" Please waiting for the search ...\n") 
    
    if(Strategy==1){
      # Search SNPs with high main effects and SNP-combinations with high interaction
      # effects. In current version, only one strategy is provided, that is, 
      # exhausative search.
      Effect <- ExhaustiveSearch(pts,class,MaxOrder,measure,0)
    }
    
    if(Strategy==2){
      # PSO search strategy
      Effect <- PSOSearch(pts,class,MaxOrder,Population,Iteration,c1=2,c2=2,TopSNP=10,measure,0)
    }
    
    # from 1-order to 5-order
    SingleEffect <- Effect$SingleEffect
    TwoEffect <- Effect$TwoEffect
    ThreeEffect <- Effect$ThreeEffect
    FourEffect <- Effect$FourEffect
    FiveEffect <- Effect$FiveEffect
    
    #############################
    # Normalization
    #############################
    # Normalization of SingleEffect, TwoEffect, ThreeEffect, FourEffect and FiveEffect
    #
    # Notice:
    #
    #       In current version, three provided evaluation measures are sensitive to
    # the orders of SNP combinations. For alleviating this situation, an normalization
    # strategy is provided.
    #
    #       If new evaluation measures provided Subsequent versions are not sensitive
    # to the orders of SNP combinations, this funcution can be commented.
    #
    
    cat("#### Normalization ####\n")
    Effect <- NormalizationEffect(MaxOrder,SingleEffect,TwoEffect,ThreeEffect,FourEffect
                                  ,FiveEffect,SaveFileName)
    
    SingleEffect <- Effect$SingleEffect
    TwoEffect <- Effect$TwoEffect
    ThreeEffect <- Effect$ThreeEffect
    FourEffect <- Effect$FourEffect
    FiveEffect <- Effect$FiveEffect
    
    #############################
    # SNP Name Notation
    #############################
    # Notation of real SNP Name
    
    cat("#### SNP Name Notation ####\n")
    Effect <- NotationName(MaxOrder,SingleEffect,TwoEffect,ThreeEffect,FourEffect
                           ,FiveEffect,SNPNames)
    
    SingleEffect <- Effect$SingleEffect
    TwoEffect <- Effect$TwoEffect
    ThreeEffect <- Effect$ThreeEffect
    FourEffect <- Effect$FourEffect
    FiveEffect <- Effect$FiveEffect
    
    #############################
    #  Collect Vertices and Edges
    #############################
    # Collect vertices and edges for network construction  
    
    cat("#### Collect Vertices and Edges ####\n")
    GraphData <- NetworkData(SingleEffect,TwoEffect,ThreeEffect,FourEffect,
                             FiveEffect,RatioThreshold,NumberThreshold)
    Edges <- GraphData$edges
    Vertices <- GraphData$vertices
    
    #############################
    # Construct Complete Graph
    #############################
    # Construct Complete Graph by using all vertices and edges
    
    cat("#### Construct Complete Graph ####\n")
    ConstructCompleteGraph(Vertices,Edges,6,SaveFileName)
    
    #############################
    #  Plot Top Effects
    #############################
    # Top effects of SNPs or SNP-combinations
    
    if (nrow(Vertices)>1){
      cat("#### Plot Top Effects ####\n")
      TpEffect <- PlotTopEffects(Vertices,20,SaveFileName)
      TopEffect <- TpEffect$TopEffect
      CombinationEffect <- TpEffect$CombinationEffect
    }else
    {
      TopEffect <- Vertices
      CombinationEffect <- Vertices
    }
    
    #############################
    #  Degree Analysis
    #############################
    # Degree Analysis of real vertices
    
    if (nrow(Vertices)>1){
      cat("#### Degree Analysis ####\n")
      Degrees <- DegreeAnalysis(Vertices,Edges,SaveFileName)
      Degrees <- Degrees$Degrees
      print(Degrees)
    }else
    {
      Degrees <- 0
    }
    
    #############################
    #  Split subgraphs
    #############################
    # Split subgraphs using walktrap.community algorithm
    
    if (nrow(Vertices)>1){
      cat("#### Split subgraphs ####\n")
      SubgraphSNPs <- SubgraphSplit(Vertices,Edges)
      SubgraphSNPs <- SubgraphSNPs$SubgroupSNPs
    }else
    {
      SubgraphSNPs <- Vertices
    }
    
    #############################
    # heatmap Factor
    #############################
    # heatmap the top factor
    
    if (nrow(Vertices)>1){
      if (is.na(SNPNameFileName[kk])){
        
        SNPCombination <- TopEffect[length(TopEffect)]
        SNPCombination <- names(SNPCombination)
        Title <- SNPCombination
        SNPCombination <- strsplit(SNPCombination,split=":")
        factor <- as.numeric(c(SNPCombination, recursive=T))
        
        HeatMapFactors <- HeatMapFactor(pts,class,factor,SaveFileName,Title)
        
      }else
      {
        SNPCombination <- TopEffect[length(TopEffect)]
        SNPCombination <- names(SNPCombination)
        Title <- SNPCombination
        SNPCombination <- strsplit(SNPCombination,split=":")
        SNPCombination <- c(SNPCombination, recursive=T)
        
        factorNum <- length(SNPCombination)
        factor <- rep(0,factorNum)
        
        SNPNamesNum <- length(SNPNames)
        
        for (i in 1:factorNum){
          for (j in 1:SNPNamesNum){
            if (SNPCombination[i]==SNPNames[j]){
              factor[i]=j
              break
            }
          }
        }
        
        HeatMapFactors<- HeatMapFactor(pts,class,factor,SaveFileName,Title)
      }
      
      HeatMapFactors <- HeatMapFactors$HeatMapFactors
    }else
    {
      HeatMapFactors <- 0
    }
    
    #############################
    # Return Results
    #############################
    # save the parameters
    parameters <- list(FileName=FileName[kk],MaxOrder=MaxOrder,RatioThreshold=
                         RatioThreshold,NumberThreshold=NumberThreshold,
                       SNPNameFileName=SNPNameFileName[kk],measure=measure)
    
    # save the data
    data <- list(pts=pts,class=class,SNPNames=SNPNames)
    
    # save the results
    results <- list(SingleEffect=SingleEffect,TwoEffect=TwoEffect,
                    ThreeEffect=ThreeEffect,FourEffect=FourEffect,FiveEffect=
                      FiveEffect)
    
    # save graphs
    graphs <- list(Edges=Edges,Vertices=Vertices,TopEffect=TopEffect,
                   Degrees=Degrees,SubgraphSNPs=SubgraphSNPs,CombinationEffect=
                     CombinationEffect,HeatMapFactors=HeatMapFactors)
    
    # return
    writeMat(paste(SaveFileName,".mat",sep=""),parameters=parameters,data=data,
             results=results,graphs=graphs)
    
    save(parameters=parameters,data=data,results=results,graphs=graphs,
         file=paste(SaveFileName,".RData",sep=""))
    
    graphics.off()
  }
  
  cat("############################\n")
  cat("#####      Finish      #####\n")
  cat("############################\n")
}
