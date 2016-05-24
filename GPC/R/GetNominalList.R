GetNominalList <- function(L, NominalSize, Rule){
    # rm(list=ls())
    # d = 3; L = 5; Growth = 'ExpOdd'; Rule = 'ClenshawCurtis';
    # d = 3; L = 5; Growth = 'ExpOdd'; Rule = 'Fejer';
    # NominalSize  <- GetNominalSize(L,Growth,Rule);

    Theta       <- matrix(data=NA,nrow=L,ncol=NominalSize$Unique[L])
    ThetaSymm   <- matrix(data=NA,nrow=L,ncol=NominalSize$UniqueSymmetric[L])
    Weights     <- array(0,dim=c(L,NominalSize$UniqueSymmetric[L],L))
    
    if (Rule == 'ClenshawCurtis'){
      for (i in 1:L){
        # cat(NominalSize$Total[i],"\n")
        NominalQuad <- WeightsRootsClenshawCurtis(NominalSize$Total[i])
        if (i == 1){
          Theta[i,1]                      <- NominalQuad$Xn; 
          ThetaSymm[i,1]                  <- NominalQuad$Xn;
          Weights[1,1,1]                  <- NominalQuad$Wn;        
        } else if (i == 2) {
          iFullSeq                        = seq(from=1,to=NominalSize$Total[i],by=2)
          iSymSeq                         = seq(from=1,to=NominalSize$UniqueSymmetric[i])
          Theta[i,(1:length(iFullSeq))]   = NominalQuad$Xn[iFullSeq]; 
          ThetaSymm[i,1]                  = NominalQuad$Xn[iSymSeq];
          Weights[i,1,2]                  = NominalQuad$Wn[1];
          Weights[1,1,2]                  = NominalQuad$Wn[2];        
        } else {
          iTotalSeq                       <- seq(from=2,to=NominalSize$Total[i]-1,by=2)
          iUniqueSeq                      <- seq(from=1,to=NominalSize$Unique[i])
          iSymSeq                         <- seq(from=1,to=NominalSize$UniqueSymmetric[i])
          Wn                              <- NominalQuad$Wn 
          Zn                              <- NominalQuad$Xn[iTotalSeq]
          Theta[i,iUniqueSeq]             <- Zn 
          ThetaSymm[i,iSymSeq]            <- Zn[iSymSeq]
          Weights[2,1,i]                  <- Wn[1];
          Wn                              <- Wn[seq(from=2,to=length(Wn)-1)]
          for (l in seq(from=i,to=3,by=-1)){
            iUniqueSeq                      <- seq(from=1,to=NominalSize$Unique[l],by=2)
            iSymSeq                         <- seq(from=1,to=NominalSize$UniqueSymmetric[i])
            Weights[l,iSymSeq,i]            <- Wn[iUniqueSeq];
            Wn                              <- Wn[seq(from=2,to=length(Wn),by=2)]
          }
        Weights[1,1,i] <- Wn;                    
        }
      }
    } else if (Rule == 'Fejer'){
      # Nested construction for Fejer
      for (i in 1:L){
        NominalQuad <- WeightsRootsFejer(NominalSize$Total[i])
        if (i == 1){
          Theta[i,1] <- NominalQuad$Xn; 
          ThetaSymm[i,1] <- NominalQuad$Xn;
          Weights[1,1,1] <- NominalQuad$Wn;        
        } else {
          Theta[i,1:NominalSize$Unique[i]] <- NominalQuad$Xn[seq(from=1,to=NominalSize$Total[i],by=2)]; 
          ThetaSymm[i,1:NominalSize$UniqueSymmetric[i]] <- NominalQuad$Xn[seq(from=1,by=2,to=NominalSize$Unique[i])];
          Wn <- NominalQuad$Wn          
          for (j in seq(from=i,by=-1,to=2)){
            Weights[j,1:NominalSize$UniqueSymmetric[j],i] <- Wn[seq(from=1,by=2,to=NominalSize$Unique[j])];
            Wn <- Wn[seq(from=2,by=2,to=NominalSize$Total[j])]         
          }
          Weights[1,1,i] <- Wn; 
        }           
      }
    } 
    res <- list(Theta=Theta,ThetaSymm=ThetaSymm,Weights=Weights) 
  return(res)
}