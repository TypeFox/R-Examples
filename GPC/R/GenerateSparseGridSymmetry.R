GenerateSparseGridSymmetry <- function(d,L,NominalList,NominalSize,ExpPoly,ParamDistrib) {
    # alpha <- ParamDistrib$alpha
    # beta <- ParamDistrib$beta
    # Getting the combinations of 'k' and 'q' needed        
    M           <- getM(d,L-1);
    GridIndex   <- indexCardinal(d,L-1)+1;    
  
    QIndexMaster = GridIndex; # re replace QIndexMaster with GridIndex
    
    # Determining the size of the compact and expanded list
    SizeCompact = 0;
    SizeExpand  = 0;
    for (m in 1:M){
        SizeCompact <- SizeCompact + prod(NominalSize$UniqueSymmetric[GridIndex[,m]]);
        SizeExpand  <- SizeExpand  + prod(NominalSize$Unique[GridIndex[,m]]);
    }

    # Some declaration of variables        
    ZnList      <- rep(1,d)
    
    KList       <- matrix(data=NA,nrow=SizeCompact,ncol=d) 
    JList       <- matrix(data=NA,nrow=SizeCompact,ncol=d)
    
    NodeCoord   = matrix(data=NA,nrow=SizeCompact,ncol=d)
    NodeWeight  = matrix(data=0,nrow=L,ncol=SizeExpand)
    

    SparseGridSym <- list(NodeCoord=NodeCoord, NodeWeight=NodeWeight, JList=JList,n=0)
    # need to change the name from SparseGridSym to SparseGridSymm

    # Generating the compact hypercube, [0 1]^d
    tempSize <- rep(0,M)
    SparseGridSym$n = 0
    for (m in 1:M){
        Ind = GridIndex[,m]
        ThetaNodes <- matrix(data=NA,nrow=d,ncol=NominalSize$UniqueSymmetric[max(Ind)])
        for (dd in 1:d){ 
            nSymm = NominalSize$UniqueSymmetric[Ind[dd]]
            ThetaNodes[dd,1:nSymm] = NominalList$Theta[Ind[dd],1:nSymm]; # should change Theta here to ThetaSymm
        } 
        PreCount <- SparseGridSym$n
        CoordDummy <- matrix(data=NA,nrow=1,ncol=d)
        JListDummy <- matrix(data=NA,nrow=1,ncol=d)
        SparseGridSym <- SmolGridGenFast(d,ThetaNodes,SparseGridSym,CoordDummy,JListDummy)   # this function changes only [NodeCoord n JList] 
        
        KList[(PreCount+1):SparseGridSym$n,] <- matrix(Ind,nrow=SparseGridSym$n-PreCount,ncol=length(Ind),byrow=TRUE)
        tempSize[m+1] = tempSize[m] + prod(NominalSize$Unique[Ind]);        
    }
    
    # Expanding [0 1]^d to [-1 1]^d by rotating the existing elements about all the d axis 
    n = 0;
    NodeNodes <- matrix(data=NA,nrow=SizeExpand,ncol=d) 
    tmpKList <- matrix(data=NA,nrow=SizeExpand,ncol=d) 
    JList <- matrix(data=NA,nrow=SizeExpand,ncol=d)
    
    NodeCoord  <- SparseGridSym$NodeCoord[1,]
    tmpKList  <- KList[1,]
    JList  <- SparseGridSym$JList[1,]

    
    for (i in 2:SparseGridSym$n){
        n <- n + 1;
        RotatedGrid <- SparseGridRotate(SparseGridSym$NodeCoord[i,],SparseGridSym$JList[i,],1);        
        if (RotatedGrid$counter != 0){
            NodeCoord <- rbind(NodeCoord,SparseGridSym$NodeCoord[i,],RotatedGrid$RotatedNodes)
            tmpKList <- rbind(tmpKList,matrix(KList[i,],nrow=RotatedGrid$counter+1,ncol=length(KList[i,]),byrow=TRUE))
            # cat(i,RotatedGrid$RotatedJList,'\n')
            JList <- rbind(JList,SparseGridSym$JList[i,],RotatedGrid$RotatedJList)       
        }
        n <- n + RotatedGrid$counter;
    }
    KList = tmpKList;
    JList <- abs(JList);
    
    for (m in 1:M){
      Ind = GridIndex[,m]
      for (K in 1:L) {
        tmpP = K+2*d-sum(Ind)-d-1
        if (tmpP >= 0) {
          QM = getM(d,tmpP)
          for (NodeInd in (tempSize[m]+1):(tempSize[m+1])) {
            QIndex <- matrix(data=NA,nrow=QM,ncol=d) 
            tempWeight = 0
            for (qm in 1:QM) {
              QIndex[qm,] <- QIndexMaster[,qm]
              mWeight = 1
              for (n in 1:d) {
                if (QIndex[qm,n]==1){mWeight<-mWeight*NominalList$Weights[KList[NodeInd,n],JList[NodeInd,n],KList[NodeInd,n]]
                } else {
                  mWeight <- mWeight*
                    (NominalList$Weights[KList[NodeInd,n],JList[NodeInd,n],KList[NodeInd,n]+QIndex[qm,n]-1] - 
                    NominalList$Weights[KList[NodeInd,n],JList[NodeInd,n],KList[NodeInd,n]+QIndex[qm,n]-2])                  
                }
              }
              tempWeight <- tempWeight + mWeight;  
            }
            NodeWeight[K,NodeInd] <- tempWeight;            
          }
        }
      }
    }
    StageSize   <- rep(0,L)
    
    for (i in 0:L-1){ StageSize[i+1] = tempSize[getM(d,i)+1]}

    for (i in 1:StageSize[L]){
      for (j in 1:d){
        if(JList[i,j]<0){ JList[i,j] = NominalSize$Unique[KList[i,j]] + 1 + JList[i,j]}
      }
    }
    # cat(d,M,StageSize[L],ExpPoly)
    # cat(GridIndex)

    Index = indexCardinal(d,L-1) # generate the basis of expansion
    PolyNodes = populateNDPhZn(d,M,StageSize[L],ExpPoly,t(NodeCoord),Index,ParamDistrib)
    
    res <- list(Node=NodeCoord,Weight=NodeWeight/2^d,Size=StageSize,NodeLevel=KList,NodePosition=JList,PolyNodes=PolyNodes)
    
    return(res)
}