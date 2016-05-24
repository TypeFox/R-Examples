NetworkData <-
function(SingleEffect,TwoEffect,ThreeEffect,FourEffect,
                        FiveEffect,RatioThreshold,NumberThreshold){
  # Collect vertices and edges for network construction
  #
  # input
  #    SingleEffect: There are 2 columns. The first column saves all SNPs, and the 
  #                  second column saves their corresponding effects. Descending 
  #                  save according to their effects.                   
  #    TwoEffect: There are 3 columns. The first 2 columns save all 2-SNP combinations,
  #               and the last column saves their corresponding effects. Descending 
  #               save according to their interaction effects.
  #    ThreeEffect: There are 4 columns. The first 3 columns save all 3-SNP combinations,
  #                 and the last column saves their corresponding effects. Descending 
  #                 save according to their interaction effects.
  #    FourEffect: There are 5 columns. The first 4 columns save all 4-SNP combinations,
  #                and the last column saves their corresponding effects. Descending 
  #                save according to their interaction effects.
  #    FiveEffect: There are 6 columns. The first 5 columns save all 5-SNP combinations,
  #                and the last column saves their corresponding effects. Descending
  #                save according to their interaction effects.
  #    RatioThreshold: The length of RatioThreshold is equal to the parameter 'MaxOrder'.
  #                    Each element is a decimal in [0,1] with a character format.
  #    NumberThreshold: The length of NumberThreshold is equal to the parameter 'MaxOrder'.
  #                     Each element is a integer with a character format.
  #
  # output
  #    vertices: Vertices for network construction
  #    edges: Vertices for network construction
  #
  # Junliang Shang
  # 3.26/2014
  
  BackupSingleEffect <- SingleEffect
  
  vertices <- data.frame(0,0,0)
  edges <- data.frame(0,0)
  
  #############################
  # Return top n SNPs as vertices from SingleEffect
  #############################
  if (ncol(SingleEffect)==2){
    SingleEffectNum <- min(min(floor(nrow(SingleEffect)*RatioThreshold[1])+1
                               ,nrow(SingleEffect)),NumberThreshold[1])
    SingleEffect <- SingleEffect[1:SingleEffectNum,]
    dim(SingleEffect) <- c(length(SingleEffect)/2,2)
    vertices <- cbind(SingleEffect,"1")
  }

  #############################
  # Return vertices and edges from TwoEffect, where vertices are SNPs in top n
  # lines, and edges represent Two-SNP effects
  #############################
  if (ncol(TwoEffect)==3){
    TwoEffectNum <- min(NumberThreshold[2],
                        min(nrow(TwoEffect),floor(nrow(TwoEffect)
                                                  *RatioThreshold[2])+1))
    TwoEffect <- TwoEffect[1:TwoEffectNum,]
    dim(TwoEffect) <- c(length(TwoEffect)/3,3)
    TwoEffect <- cbind(paste(TwoEffect[,1],":",TwoEffect[,2],sep=""),TwoEffect)
    vertices <- cbind(BackupSingleEffect[
      match(unique(c(TwoEffect[,c(2,3)],vertices[,1])),
            BackupSingleEffect[,1]),],"1")
    edges <- TwoEffect[,c(1,2)]
    edges <- rbind(edges,TwoEffect[,c(1,3)])
  }
  
  #############################
  # Return vertices and edges from ThreeEffect, where vertices are SNPs in top n
  # lines, and edges represent Three-SNP effects
  #############################
  if (ncol(ThreeEffect)==4){
    ThreeEffectNum <- min(NumberThreshold[3],
                          min(nrow(ThreeEffect),floor(nrow(ThreeEffect)
                                                      *RatioThreshold[3])+1))  
    ThreeEffect <- ThreeEffect[1:ThreeEffectNum,]
    dim(ThreeEffect) <- c(length(ThreeEffect)/4,4)
    ThreeEffect <- cbind(paste(ThreeEffect[,1],":",ThreeEffect[,2],
                               ":",ThreeEffect[,3],sep=""),ThreeEffect)
    vertices <- cbind(BackupSingleEffect[
      match(unique(c(ThreeEffect[,c(2,3,4)],vertices[,1])),
            BackupSingleEffect[,1]),],"1")
    edges <- rbind(edges,ThreeEffect[,c(1,2)])
    edges <- rbind(edges,ThreeEffect[,c(1,3)])
    edges <- rbind(edges,ThreeEffect[,c(1,4)])
  }
  
  #############################
  # Return vertices and edges from FourEffect, where vertices are SNPs in top n
  # lines, and edges represent Four-SNP effects
  #############################
  if (ncol(FourEffect)==5){
    FourEffectNum <- min(NumberThreshold[4],
                          min(nrow(FourEffect),floor(nrow(FourEffect)
                                                      *RatioThreshold[4])+1))  
    FourEffect <- FourEffect[1:FourEffectNum,]
    dim(FourEffect) <- c(length(FourEffect)/5,5)
    FourEffect <- cbind(paste(FourEffect[,1],":",FourEffect[,2],
                               ":",FourEffect[,3],":",FourEffect[,4],sep=""),FourEffect)
    vertices <- cbind(BackupSingleEffect[
      match(unique(c(FourEffect[,c(2,3,4,5)],vertices[,1])),
            BackupSingleEffect[,1]),],"1")
    edges <- rbind(edges,FourEffect[,c(1,2)])
    edges <- rbind(edges,FourEffect[,c(1,3)])
    edges <- rbind(edges,FourEffect[,c(1,4)])
    edges <- rbind(edges,FourEffect[,c(1,5)])
  }
  
  #############################
  # Return vertices and edges from FiveEffect, where vertices are SNPs in top n
  # lines, and edges represent Five-SNP effects
  #############################
  if (ncol(FiveEffect)==6){
    FiveEffectNum <- min(NumberThreshold[5],
                         min(nrow(FiveEffect),floor(nrow(FiveEffect)
                                                    *RatioThreshold[5])+1))  
    FiveEffect <- FiveEffect[1:FiveEffectNum,]
    dim(FiveEffect) <- c(length(FiveEffect)/6,6)
    FiveEffect <- cbind(paste(FiveEffect[,1],":",FiveEffect[,2],
                              ":",FiveEffect[,3],":",FiveEffect[,4],
                              ":",FiveEffect[,5],sep=""),FiveEffect)
    vertices <- cbind(BackupSingleEffect[
      match(unique(c(FiveEffect[,c(2,3,4,5,6)],vertices[,1])),
            BackupSingleEffect[,1]),],"1")
    edges <- rbind(edges,FiveEffect[,c(1,2)])
    edges <- rbind(edges,FiveEffect[,c(1,3)])
    edges <- rbind(edges,FiveEffect[,c(1,4)])
    edges <- rbind(edges,FiveEffect[,c(1,5)])
    edges <- rbind(edges,FiveEffect[,c(1,6)])
  }
  
  
  #############################
  # Return two-SNP effects as virtual vertices
  #############################
  if (ncol(TwoEffect)==4){
    vertices <- rbind(vertices,cbind(TwoEffect[,c(1,4)],"2"))
  }
  
  #############################
  # Return three-SNP effects as virtual vertices
  #############################
  if (ncol(ThreeEffect)==5){
    vertices <- rbind(vertices,cbind(ThreeEffect[,c(1,5)],"3"))
  }
  
  #############################
  # Return Four-SNP effects as virtual vertices
  #############################
  if (ncol(FourEffect)==6){
    vertices <- rbind(vertices,cbind(FourEffect[,c(1,6)],"4"))
  }
  
  #############################
  # Return Five-SNP effects as virtual vertices
  #############################
  if (ncol(FiveEffect)==7){
    vertices <- rbind(vertices,cbind(FiveEffect[,c(1,7)],"5"))
  }
  
  colnames(vertices)=c("id","value","label")
  colnames(edges)=c("From","To")
  
  list(vertices=vertices,edges=edges) 
}
