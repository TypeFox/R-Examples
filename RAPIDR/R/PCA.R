#############################################################
## Internal functions - not exported 
#############################################################

correct.counts.with.PCA<- function ( PCA_output, ratios.mat ) { 
  
  message("correcting counts with PCA components")
  PCA.rotations  <- PCA_output[["PCA.rotations"]]
  PCA.centers    <- PCA_output[["PCA.centers"]] 
  coeff.mat      <- PCA_output[["coeff.mat"]] 
  numPC          <- PCA_output[["numPC"]] 
 
  n.bins <- ncol(ratios.mat)
  #ratios.mat<-matrix(data=NA,nrow=nrow(binned.counts),ncol=ncol(binned.counts))
  
  for (i in 1:n.bins){
    ratios.mat[,i] <- ratios.mat[,i]-PCA.centers[i]
  }
  
  message("Making PCA.scores")
  PCA.scores <- t(t(PCA.rotations) %*% t(ratios.mat))
  
  PCA.var<-data.frame(intercept=rep(1)) 
  for(i in 1:numPC) { 
      pcname = paste("PCA", i, sep = "") 
      PCA.var <- cbind(PCA.var, PCA.scores[,i]) 
      names(PCA.var)[1+i] <- pcname 
  }
 
  #PCA.var<-data.frame(intercept=rep(1), PCA1=c(PCA.scores[,1]), PCA2=c(PCA.scores[,2]), 
  #                    PCA3=c(PCA.scores[,3]), PCA4=c(PCA.scores[,4]),
  #                    PCA5 = c(PCA.scores[,5]), PCA6 = c(PCA.scores[,6]), 
  #                    PCA7 = c(PCA.scores[,7]), PCA8 = c(PCA.scores[,8]), 
  #                    PCA9 = c(PCA.scores[,9]), PCA10 = c(PCA.scores[,10]))
  
  rm(PCA.scores)
  rm(PCA.rotations)
  PCA.var <- data.matrix(PCA.var)
  #ratios.clean<-matrix(data=NA,nrow=nrow(binned.counts),ncol=ncol(binned.counts))
  
  message("calc responses")
  responses <- PCA.var %*% coeff.mat
  
  for (i in 1:n.bins){
    #message(i,  ' out of ', nrow(ratios.mat))
    residuals <- ratios.mat[,i] - responses[,i]
    ratios.mat[,i] <- residuals + coeff.mat[1,i]
  }
  
  # Add back the center values 
  for (i in 1:n.bins){
    ratios.mat[,i] <- ratios.mat[,i]+PCA.centers[i]
  }
  
  return(ratios.mat)
}

doPCA <- function (data.mat.clean.for.PCA, numPC = 10 ) {
  
  message("doing PCA...")
  my.pca        <- prcomp(data.mat.clean.for.PCA)  
  PCA.scores    <- my.pca$x 
  
  #save(PCA.scores, PCA.centers, PCA.rotations, binned.ratios, outcomes, sampleIDs, total.counts, file="~/UCL/PhaseI_all/P1_733_gccounts.Rdata")    
  #load("~/UCL/PhaseI_all/P1_733_gccounts.Rdata")
  
  message("doing the linear regression modelling...")

  PCA.var<-data.frame(intercept=rep(1)) 
  for(i in 1:numPC) { 
      pcname = paste("PCA", i, sep = "") 
      PCA.var <- cbind(PCA.var, PCA.scores[,i]) 
      names(PCA.var)[1+i] <- pcname 
  } 
  
  #PCA.var<-data.frame(intercept=rep(1), PCA1=c(PCA.scores[,1]), PCA2=c(PCA.scores[,2]), 
  #                    PCA3=c(PCA.scores[,3]), PCA4=c(PCA.scores[,4]),
  #                    PCA5 = c(PCA.scores[,5]), PCA6 = c(PCA.scores[,6]), 
  #                    PCA7 = c(PCA.scores[,7]), PCA8 = c(PCA.scores[,8]), 
  #                    PCA9 = c(PCA.scores[,9]), PCA10 = c(PCA.scores[,10]))
  
  PCA.mat <- data.matrix(PCA.var)
  reg.mat<-solve(t(PCA.mat) %*% PCA.mat)%*%t(PCA.mat)
  
  coeff.mat<-reg.mat %*% data.mat.clean.for.PCA    
  #save(coeff.mat, file= "~/UCL/PhaseI_all/P1_733_coeff.Rdata")
  
  PCA_output <- list() 
  PCA_output[["PCA.scores"]] <- my.pca$x 
  PCA_output[["PCA.centers"]] <- my.pca$center
  PCA_output[["PCA.rotations"]] <- my.pca$rotation
  PCA_output[["coeff.mat"]] <- coeff.mat 
  PCA_output[["numPC"]]     <- numPC 
 
  return(PCA_output)
  
}
