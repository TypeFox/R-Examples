Sel.Features <- function(ES, Y, K = "Min", Verbose = FALSE){
  Y <- as.factor(Y)
  levels(Y) = c(1,2)
  Core   <- CI.emprical(ES, Y)
  G.Mask <- GMask(ES, Core, Y)
  Pos    <- POS(ES, Core, Y)
  Min.S  <- Min.Set(G.Mask, Pos)
  if (K == "Min"){
    MinSub <- cbind(Feature = Min.S$Min.Subset, Pos = Pos[Min.S$Min.Subset])
    return(list(Features = MinSub, Covered.Obs = Min.S$Covered.Obs))
  }
  else{
    R.DC         <- RDC(G.Mask, Y)
    Net.Features <- cbind(Feature = c(1:length(Pos)), Pos = Pos, RDC = R.DC)[-Min.S$Min.Subset,]
    # dividing feauters into two ordered groups (one for each class)
    Cl1.Features <- Net.Features[Net.Features[,"RDC"] == 1,][order(Net.Features[Net.Features[,"RDC"] == 1,][,"Pos"]), ]
    Cl2.Features <- Net.Features[Net.Features[,"RDC"] == 2,][order(Net.Features[Net.Features[,"RDC"] == 2,][,"Pos"]), ]
    # craete an empty vector to save list of ranked features within
    G.Rank <- vector(mode = "integer", length = nrow(Cl1.Features) + nrow(Cl2.Features))
    # get the start position for class whose features' number is less - 1st index in "G.Rank" is for less Pos measure
    Start <- order(c(Cl1.Features[1,"Pos"], Cl2.Features[1,"Pos"]))[which.min(c(nrow(Cl1.Features), nrow(Cl2.Features)))]
    # generate an indices sequence for class whose features' number is less (those indices are used to apply the alternating fashion)
    Squns <- seq(from = Start, by = 2, length.out = min(nrow(Cl1.Features), nrow(Cl2.Features)))
    
    if(which.min(c(nrow(Cl1.Features), nrow(Cl2.Features))) == 1){
      G.Rank[Squns]  <- Cl1.Features[,"Feature"]
      G.Rank[-Squns] <- Cl2.Features[,"Feature"]
    }
    else{
      G.Rank[Squns]  <- Cl2.Features[,"Feature"]
      G.Rank[-Squns] <- Cl1.Features[,"Feature"]
    }
    
    Final.Rank <- c(Min.S$Min.Subset, G.Rank)[1:K]
    nMin.Genes <- length(Min.S$Min.Subset)
    if (Verbose){
      Status  <- c(rep("Min.Set", nMin.Genes), R.DC[G.Rank])[1:K]
      Details <- data.frame(Features = Final.Rank, Pos = Pos[Final.Rank], Status = Status)
      return(list(Features = Final.Rank, nMin.Features = nMin.Genes, Measures = Details))
    }
    else{return(list(Features = Final.Rank, nMin.Features = nMin.Genes))}
  }
}