CINOEDV_Main <-
function(){
  # Main function of R package CINOEDV (Co-Information based N-Order
  # Epistasis Detector and Visualizer).
  #
  # author: Junliang Shang
  #
  # Email: jlshang@mail.xidian.edu.cn
  #        shangjunliang110@163.com
  #
  # 3.24/2014
  
  #############################
  # Welcome
  #############################
  cat("#### Welcome ####\n\n")
  cat(" Hello, Welcome to CINOEDV !\n")
  cat(" CINOEDV: Co-Information based N-Order Epistasis Detector and Visualizer.\n\n")
  rm(list=ls())
  
  #############################
  # Install Packages
  #############################  
  cat("#### Install Packages ####\n\n")
  InstallPackage()
  cat(" Install Packages OK !\n\n") 
  
  #############################
  # Set and check Parameters
  #############################
  # Two input ways are provided in this version.
  # First, they are setted by the console.
  # Second, they are setted by a file "Parameters.R", which is recommended.
  
  cat("#### Set and Check Parameters ####\n\n")
  # cat(" Please select the input way that you wish to use (1/2).\n")
  # cat(" 1: The Console.\n")
  # cat(" 2: The Parameter file (Recommendation Option).\n")
  # way <- readline()
  #
  # # The parameters are setted by the console.
  # if (way=="1"){
  #   cat(" Please set parameters using the console.\n\n")
    
    ########################
    # File name with (.mat) format that saves SNP data.
    # The file has two variables, i.e., pts and class.
    #
    # pts: Row -> Sample, Column -> SNP
    #      1 -> AA
    #      2 -> Aa
    #      3 -> aa
    #
    # class: Row -> 1, Column -> class label
    #      1 -> case
    #      2 -> control
    #
    # for example, test.mat
    
    cat(" Please input the file name with its format (.mat) that saves SNP data.\n")
    cat(" The file has two variables, i.e., pts and class.\n\n")
    cat(" pts: Row -> Sample, Column -> SNP\n")
    cat("    1 -> AA\n")
    cat("    2 -> Aa\n")
    cat("    3 -> aa\n\n")
    cat(" class: Row -> 1, Column -> class label\n")
    cat("    1 -> case\n")
    cat("    2 -> control\n\n")
    cat(" for example, test.mat\n")
    FileName <- readline()
    
    # Input and Check SNP data
    Data <- InputData(FileName)
    pts <- Data$pts
    class <- Data$class
  
  ########################
  # If there are real SNP names which will be used for constructing graphs and further
  # analysis, the name of file that saves real SNP names should be provided.
  
  cat(" If there are real SNP names which will be used for constructing graphs and further
       analysis, the name of file that saves real SNP names should be provided. \n\n")
  cat(" Please input the name of such file with (.mat) format.\n")
  cat(" The file has only one variable, i.e., Name.\n")
  cat(" Name: Row -> 1, Column -> SNP Name\n\n")
  cat(" If not exist such file, please input NA\n\n")
  cat(" For example, test_Name.mat\n\n")
  SNPNameFileName <- readline()
  if (SNPNameFileName=="NA"){
    SNPNameFileName <- NA  
  } 
  SNPNames <- TestSNPNameFile(ncol(pts),SNPNameFileName)
  SNPNames <- SNPNames$SNPNames
  
    ########################
    # The specified maximum order, must be setted as 2,3,4 or 5.
    
    cat(" Please input the maximum order (2/3/4/5), and 3 is the Recommendation Option.\n")
    MaxOrder <- readline()
    TestMaxOrder(MaxOrder)
    MaxOrder <- as.numeric(MaxOrder)
    
    ########################
    # Ratio Thresholds.
    #
    # Ratio thresholds control the numbers of retained top SNPs or SNP-combinations.
    # For each order, from 1 to MaxOrder (User Specified), there is a ratio threshold.
    # 
    # for example
    # RatioThreshold[1] means that top (RatioThreshold[1]*nrow(SingleEffect)) SNPs 
    # with the biggest main effects should be retained.
    # RatioThreshold[2] means that top (RatioThreshold[2]*nrow(TwoEffect)) SNP-combinations
    # with the biggest 2-order interaction effects should be retained.
    # RatioThreshold[3] means that top (RatioThreshold[3]*nrow(ThreeEffect)) SNP-combinations
    # with the biggest 3-order interaction effects should be retained.
    # RatioThreshold[4] means that top (RatioThreshold[4]*nrow(FourEffect)) SNP-combinations
    # with the biggest 4-order interaction effects should be retained.
    # RatioThreshold[5] means that top (RatioThreshold[5]*nrow(FiveEffect)) SNP-combinations
    # with the biggest 5-order interaction effects should be retained.
    #
    # Each of them should be setted in [0,1], By default, it is setted as 1.
    #
    # However, this is not an exact control. Parameters 'RatioThreshold' and 'NumberThreshold'
    # are used together to control them.
    #
    # That is,
    # min(NumberThreshold[1],RatioThreshold[1]*nrow(SingleEffect))
    # min(NumberThreshold[2],RatioThreshold[2]*nrow(TwoEffect))
    # min(NumberThreshold[3],RatioThreshold[3]*nrow(ThreeEffect)) 
    # min(NumberThreshold[4],RatioThreshold[4]*nrow(FourEffect))
    # min(NumberThreshold[5],RatioThreshold[5]*nrow(FiveEffect))
    #
    # Using above settings, the numbers of virtual Vertexes denoting high order epistatic 
    # interactions in the graph are clear. Nevertheless, the number of Real Vertexes (i.e., the
    # number of SNPs) is unclear, which is obviously more than the above set, since some of them
    # must be included to connect virtual vertexes.
    
    cat(" Ratio thresholds control the numbers of retained top SNPs or SNP-combinations.\n")
    cat(" For each order, there is a ratio threshold.\n")
    cat(" Each of them should be setted in [0,1], By default, it is setted as 1.\n\n")
    cat(" For example:\n")
    cat("    RatioThreshold[1] means that top (RatioThreshold[1]*nrow(SingleEffect)) SNPs
        with the biggest main effects should be retained.\n")
    cat("    RatioThreshold[3] means that top (RatioThreshold[3]*nrow(ThreeEffect)) SNP-combinations
        with the biggest 3-order interaction effects should be retained.\n\n")
        
    cat(" Please input ",MaxOrder," ratio thresholds.\n\n")
    
    RatioThreshold <- character(length=MaxOrder)
    
    # input ratio thresholds
    for (i in 1:MaxOrder){
      cat(" Please input the ",i," ratio threshold.\n")
      RatioThreshold[i] <- readline()
    }
    TestRatioThreshold(MaxOrder,RatioThreshold)
    RatioThreshold <- as.numeric(RatioThreshold)
    
    ########################
    # Number Thresholds.
    #
    # Number thresholds control the numbers of retained top SNPs or SNP-combinations.
    # For each order, from 1 to MaxOrder (User Specified), there is a number threshold.
    #
    # for example
    # NumberThreshold[1] means that top NumberThreshold[1] SNPs with the biggest main
    # effects should be retained.
    # NumberThreshold[2] means that top NumberThreshold[2] SNP-combinations with the 
    # biggest 2-order interaction effects should be retained.
    # NumberThreshold[3] means that top NumberThreshold[3] SNP-combinations with the 
    # biggest 3-order interaction effects should be retained.
    # NumberThreshold[4] means that top NumberThreshold[4] SNP-combinations with the 
    # biggest 4-order interaction effects should be retained.
    # NumberThreshold[5] means that top NumberThreshold[5] SNP-combinations with the 
    # biggest 5-order interaction effects should be retained.
    #
    # Each of them should be setted as an integer, By default, it is setted as 10.
    #
    # However, this is not an exact control. Parameters 'RatioThreshold' and 'NumberThreshold'
    # are used together to control them.
    #
    # That is,
    # min(NumberThreshold[1],RatioThreshold[1]*nrow(SingleEffect))
    # min(NumberThreshold[2],RatioThreshold[2]*nrow(TwoEffect))
    # min(NumberThreshold[3],RatioThreshold[3]*nrow(ThreeEffect)) 
    # min(NumberThreshold[4],RatioThreshold[4]*nrow(FourEffect))
    # min(NumberThreshold[5],RatioThreshold[5]*nrow(FiveEffect))
    #
    # Using above settings, the numbers of virtual Vertexes denoting high order epistatic 
    # interactions in the graph are clear. Nevertheless, the number of Real Vertexes (i.e., the
    # number of SNPs) is unclear, which is obviously more than the above set, since some of them
    # must be included to connect virtual vertexes.
    
    cat(" Number thresholds control the numbers of retained top SNPs or SNP-combinations.\n")
    cat(" For each order, there is a number threshold.\n")
    cat(" Each of them should be setted as an integer, By default, it is setted as 10.\n\n")
    cat(" For example:\n")
    cat("    NumberThreshold[1] means that top NumberThreshold[1] SNPs with the biggest main
        effects should be retained.\n")
    cat("    NumberThreshold[3] means that top NumberThreshold[3] SNP-combinations with the
        biggest 3-order interaction effects should be retained.\n\n")
    
    cat(" Please input ",MaxOrder," number thresholds.\n\n")
    
    NumberThreshold <- character(length=MaxOrder)
    
    # input number thresholds
    for (i in 1:MaxOrder){
      cat(" Please input the ",i," number threshold.\n")
      NumberThreshold[i] <- readline()
    }
    TestNumberThreshold(MaxOrder,NumberThreshold)
    NumberThreshold <- as.numeric(NumberThreshold)

    ########################
    # measure: the label of current used evaluation measure
    #         1 -> The classic co-information measure
    #         2 -> The Normalized co-information measure
    #         3 -> TingHu's Co-Information Measure
    #         others -> The classic co-information measure
    
    cat(" Please select the evaluation measure (1/2/3), and 1 is the Recommendation Option.\n")
    cat("    1: The Classic Co-Information Measure\n")
    cat("    2: The Normalized Co-Information Measure\n")
    cat("    3: TingHu's Co-Information Measure\n")
    cat("    Others: Also the Classic Co-Information Measure\n")
    
    measure <- readline()
    if (measure %in% c("1","2","3")){
      measure <- as.numeric(measure)
    }else
    {
      measure <- 1
    }
  
  ########################
  # Strategy: Search Strategies
  #         1 -> The exhaustive strategy
  #         2 -> The PSO based strategy
  #         others -> The classic co-information measure
  
  cat(" Please select the search strategy (1/2), and 1 is the Recommendation Option.\n")
  cat("    1: The exhaustive strategy\n")
  cat("    2: The PSO strategy\n")
  cat("    Others: Also the exhaustive strategy\n")
  
  Strategy <- readline()
  if (Strategy %in% c("1","2")){
    Strategy <- as.numeric(Strategy)
  }else
  {
    Strategy <- 1
  }
  
  ########################
  # parameters of PSO based strategy, such as, Population, Iteration 
  if (Strategy==2){
    # Population: numeric. The number of particles.
    
    cat("   Please input the number of particles, and 1000 is the Recommendation option.\n")
    Population <- readline()
    
    # Iteration: numeric. The number of iterations.
    cat("   Please input the number of iterations, and 100 is the Recommendation option.\n")
    Iteration <- readline()
    
    TestPSOParameters(Population,Iteration)
    
    Population <- as.numeric(Population)
    Iteration <- as.numeric(Iteration)
  }
  
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
    SaveFileName <- paste(substr(FileName,1,nchar(FileName)-4),MaxOrder,
                          RatioN,NumberN,Strategy,sep="_")
  }
  if (Strategy==2) {
    SaveFileName <- paste(substr(FileName,1,nchar(FileName)-4),MaxOrder,
                          RatioN,NumberN,Strategy,Population,Iteration,sep="_")
  }
  
  #############################
  # Search Strategies
  #############################
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
    # print(Degrees)
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
    if (is.na(SNPNameFileName)){
      
      SNPCombination <- CombinationEffect[length(CombinationEffect)]
      SNPCombination <- names(SNPCombination)
      Title <- SNPCombination
      SNPCombination <- strsplit(SNPCombination,split=":")
      factor <- as.numeric(c(SNPCombination, recursive=T))
      
      HeatMapFactors <- HeatMapFactor(pts,class,factor,SaveFileName,Title)
      
    }else
    {
      SNPCombination <- CombinationEffect[length(CombinationEffect)]
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
  parameters <- list(FileName=FileName,MaxOrder=MaxOrder,RatioThreshold=
                       RatioThreshold,NumberThreshold=NumberThreshold,
                     SNPNameFileName=SNPNameFileName,measure=measure)
  
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
  
  cat("############################\n")
  cat("#####      Finish      #####\n")
  cat("############################\n")
  
  list(parameters=parameters,data=data,results=results,graphs=graphs)
}
