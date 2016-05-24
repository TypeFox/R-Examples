mgVarPart <- function (genD, vectorsMEM, coords, perm=1000) {

    X <- as.matrix(vectorsMEM)
    X <- apply(X, 2, scale)
    W <- as.matrix(coords)
    W <- apply(W, 2, scale)
    genD <- as.matrix(genD)
    Number_Predictors_X <- ncol(X)
    Number_Predictors_W <- ncol(W)
    n<-nrow(X)
    
    XW <- as.matrix(cbind(X,W));
    result <- mgRDA(genD, XW, full=FALSE)
    abc <- result$RsqAdj
    Fabc <- result$F
    
    result <- mgRDA(genD, X, full=FALSE)
    ab <- result$RsqAdj
    Fab <- result$F
    residuals_X <- result$res_matrix 
    predicted_X <- result$pred_matrix
    
    result <- mgRDA(genD, W, full=FALSE)
    bc <- result$RsqAdj
    Fbc <- result$F
    residuals_W <- result$res_matrix 
    predicted_W <- result$pred_matrix
    
    # unique fraction of contribution related to X
    a <- abc - bc
    # unique fraction of contribution related to W
    c <- abc- ab
    # common fraction of contribution between X and W
    b <- abc - a - c
    # residual fraction
    d <- 1-abc
    
    Fa=(a/Number_Predictors_X)/(d/(n-Number_Predictors_X-Number_Predictors_W));
    Fc=(c/Number_Predictors_W)/(d/(n-Number_Predictors_X-Number_Predictors_W));
    
    Prob_abc=1/perm; Prob_ab=1/perm; Prob_bc=1/perm; Prob_a=1/perm; Prob_c=1/perm;
    
    # permutations test
    for (i in 1:(perm-1)) {
      # testing fraction a; notice that we permute the residual values in W and not in X
      permuted_rows=sample(n,replace=FALSE)
      # permuting the residual matrix, which is from the distance, and hence the need to permute
      # rows and columns in the same way, hence the use of permuted_rows for columns and rows below
      # testing fraction a
      # Yperm=predicted_W+residuals_W[permuted_rows,permuted_rows] # implement permutation of residuals in the future
      
      # testing fraction a
      result <- mgRDA(genD, XW[permuted_rows,], full=FALSE)
      abcRnd <- result$RsqAdj
      FabcRnd <- result$F
      result <- mgRDA(genD, W[permuted_rows,], full=FALSE)
      bcRnd <- result$RsqAdj
      FbcRnd <- result$F
      aRnd=abcRnd-bcRnd;
      dRnd=1-abcRnd;
      FaRnd=(aRnd/Number_Predictors_X)/(dRnd/(n-Number_Predictors_X-Number_Predictors_W));
      if (FaRnd >= Fa) {Prob_a<-Prob_a+1/perm}
      
      # testing fraction c
      result <- mgRDA(genD, X[permuted_rows,], full=FALSE)
      abRnd <- result$RsqAdj
      FabRnd <- result$F
      cRnd=abcRnd-abRnd;
      FcRnd=(cRnd/Number_Predictors_W)/(dRnd/(n-Number_Predictors_X-Number_Predictors_W));
      if (FcRnd >= Fc) {Prob_c<-Prob_c+1/perm}
      
      # testing abc
      if (FabcRnd >= Fabc) {Prob_abc<-Prob_abc+1/perm}
      # testing ab
      if (FabRnd >= Fab) {Prob_ab<-Prob_ab+1/perm}
      # testing bc
      if (FbcRnd >= Fbc) {Prob_bc<-Prob_bc+1/perm}
      
    }
    
    result <- mat.or.vec(7,2)
    
    result[1,1] <- abc
    result[2,1] <- ab
    result[3,1] <- bc
    result[4,1] <- a
    result[5,1] <- c
    result[6,1] <- b
    result[7,1] <- d
    
    result[1,2] <- Prob_abc
    result[2,2] <- Prob_ab
    result[3,2] <- Prob_bc
    result[4,2] <- Prob_a
    result[5,2] <- Prob_c
    result[6,2] <- NA
    result[7,2] <- NA
    
    colnames(result) <- c("Estimate","p-value")
    rownames(result) <- c("[abc]","[ab]","[bc]","[a]","[c]","[b]","[d]")
    return(result)
}