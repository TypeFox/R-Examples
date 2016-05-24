turnover <-
function(lab.mat, dist.mat, type =c("quart", "octal")){
     ## whether the input "lab.mat" is matrix
     if(!is.matrix(lab.mat)){
         stop("The \"grid.lab.mat\" object must be matrix\n
              represent the names for each plot!")
     }
     if(!is.matrix(dist.mat)){
         dist.mat <- as.matrix(dist.mat)
     }
     plotnames.lab <- as.vector(lab.mat)
     plotnames.dist <- rownames(dist.mat)
     if(!all(plotnames.lab %in% plotnames.dist)){
         stop("the distance matrix does not represents all\n 
              the distances amonge the grids.")
     }
     type <- match.arg(type)
     res <- rep(NA, nrow(lab.mat)*ncol(lab.mat))
     dim(res) <- c(nrow(lab.mat),ncol(lab.mat))
     
     for( i in 2:(nrow(lab.mat)-1)){
         for(j in 2:(ncol(lab.mat)-1)){  
          ## focal plot
          target    <- lab.mat[i,j]
          plotname1 <- lab.mat[i-1,j-1]
          plotname2 <- lab.mat[i-1,j]
          plotname3 <- lab.mat[i-1,j+1]
          plotname4 <- lab.mat[i,j-1]
          plotname5 <- lab.mat[i,j+1]
          plotname6 <- lab.mat[i+1,j-1]
          plotname7 <- lab.mat[i+1,j]
          plotname8 <- lab.mat[i+1,j+1]
          rowse <- dist.mat[which(rownames(dist.mat) == target),]
     
          dist1 <- rowse[which(names(rowse)== plotname1)]
          dist2 <- rowse[which(names(rowse)== plotname2)]
          dist3 <- rowse[which(names(rowse)== plotname3)]
          dist4 <- rowse[which(names(rowse)== plotname4)]
          dist5 <- rowse[which(names(rowse)== plotname5)]
          dist6 <- rowse[which(names(rowse)== plotname6)]
          dist7 <- rowse[which(names(rowse)== plotname7)]
          dist8 <- rowse[which(names(rowse)== plotname8)]
     
          if(type == "octal"){
              res[i,j] <- mean(c(dist1,dist2,dist3,dist4,dist5,dist6,dist7,dist8))
            }
          else{   
               if(type == "quart"){
                   res[i,j] <- mean(c(dist2,dist4,dist5,dist7))
                 }
              } 
          }
      }
return(res)
}

