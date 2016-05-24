
generate.normal.statistics <- function(resp.var, null, raw.geno, pathway, family, lambda){
  
  pathway <- pathway[pathway$SNP %in% colnames(raw.geno), ]
  
  chr <- sort(unique(pathway$Chr))
  G <- list()
  for(i in 1:length(chr)){
    snps <- pathway$SNP[pathway$Chr == chr[i]]
    id <- which(colnames(raw.geno) %in% snps)
    G[[i]] <- as.matrix(raw.geno[, id, drop = FALSE])
    setDT(raw.geno)
    rs <- colnames(raw.geno)[id]
    for(r in rs){
      raw.geno[, r := NULL, with = FALSE]
    }
    class(raw.geno) <- "data.frame"
    #raw.geno <- raw.geno[, -id, drop = FALSE]
    gc()
  }
  rm(raw.geno)
  gc()
  
  formula <- paste(resp.var, "~ . -1")
  
  mdl.null <- glm(formula, data = null, family = family)
  nchr <- length(G)
  
  X <- as.matrix(null[, -1, drop = FALSE])
  V <- list()
  score0 <- list()
  if(family == "binomial"){
    y.hat <- mdl.null$fitted.values
    r <- null[, 1] - y.hat
    A <- y.hat * (1 - y.hat)
    for(i in 1:nchr){
      tmp <- try(V[[i]] <- t(G[[i]])%*% (A * G[[i]]) - t(G[[i]])%*% (A * X) %*% solve(t(X) %*% (A * X)) %*% t(X) %*% (A * G[[i]]), silent = TRUE)
      
      if(error.try(tmp)){
        msg <- "Potential existence of multicollinearity detected and ARTP2 cannot automatically deal with it right now. Please check your covariates specified in formula"
        stop(msg)
      }
      
      score0[[i]] <- as.vector(t(G[[i]]) %*% r)
      G[[i]] <- NA
      gc()
      V[[i]] <- V[[i]] / nrow(X)
      score0[[i]] <- score0[[i]] / sqrt(nrow(X)) / sqrt(lambda)
      names(score0[[i]]) <- colnames(V[[i]])
    }
  }else{
    r <- mdl.null$residuals
    s2 <- sum(r^2)/(length(r)-ncol(X))
    for(i in 1:nchr){
      tmp <- try(V[[i]] <- (t(G[[i]]) %*% G[[i]] - t(G[[i]]) %*% X %*% solve(t(X) %*% X) %*% t(X) %*% G[[i]]) / s2, silent = TRUE)
      
      if(error.try(tmp)){
        msg <- "Potential existence of multicollinearity detected and ARTP2 cannot automatically deal with it right now. Please check your covariates specified in formula"
        stop(msg)
      }
      
      score0[[i]] <- as.vector(t(G[[i]]) %*% r / s2) / sqrt(lambda)
      G[[i]] <- NA
      gc()
      V[[i]] <- V[[i]] / nrow(X)
      score0[[i]] <- score0[[i]] / sqrt(nrow(X))
      names(score0[[i]]) <- colnames(V[[i]])
    }
  }
  
  for(i in 1:length(V)){
    rs <- sort(names(score0[[i]]))
    score0[[i]] <- score0[[i]][rs]
    V[[i]] <- V[[i]][rs, rs, drop = FALSE]
  }
  
  names(V) <- as.character(chr)
  names(score0) <- as.character(chr)
  
  list(V = V, score0 = score0)
  
}

