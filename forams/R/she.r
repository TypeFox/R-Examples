`she` <-
  function (df, method = "abundance") {
    
    SHE <- new("she")
    
    if (!is.na(pmatch(method, "abundance")))
      method <- "abundance"
    METHODS <- c("abundance", "frequency")
    method <- pmatch(method, METHODS)
    if (is.na(method))
      stop("Invalid method, type 'abundance' or 'frequency'.")
    if (method == 1) {
      S <- NA
      H <- NA
      E <- NA
      N <- NA
      S[1] <- as.matrix(specnumber(t(df)))[1]
      H[1] <- as.matrix(diversity(t(df), 'shannon'))[1]
      E[1] <- exp(H[1]) / S[1]
      N[1] <- as.matrix(colSums(df))[1]
      for(i in 2:ncol(df)) {
        S[i] <- as.matrix(specnumber(t(cbind(matrix(rowSums(df[1:i]),
                                                    dimnames = list(c(rownames(df)), paste('AC', i, sep = ''))),
                                             df[-1:-i]))))[1]
        H[i] <- as.matrix(diversity(t(cbind(matrix(rowSums(df[1:i]),
                                                   dimnames = list(c(rownames(df)), paste('AC', i, sep = ''))),
                                            df[-1:-i])), 'shannon'))[1]
        E[i] <- exp(H[i]) / S[i]
        N[i] <- as.matrix(colSums(cbind(matrix(rowSums(df[1:i]),
                                               dimnames = list(c(rownames(df)), paste('AC', i, sep=''))),
                                        df[-1:-i])))[1]  
      }             
      SHE@bi <- as.data.frame(cbind(S, H, E, N))
      colnames(SHE@bi) <- c('S', 'H', 'E', 'N')
      rownames(SHE@bi) <- colnames(df)
      
      return(SHE)
    }
    if (method == 2) {
      XF <- as.data.frame(matrix(nrow=nrow(df), ncol=ncol(df)))
      rownames(XF) <- rownames(df)
      for(i in 1:ncol(df)) {
        XF[, i] <- as.matrix(df[1:nrow(df), i] / colSums(df)[i])
      }
      S <- NA
      H <- NA
      E <- NA
      L <- NA
      S[1] <- as.matrix(specnumber(t(XF)))[1]
      H[1] <- as.matrix(diversity(t(XF), 'shannon'))[1]
      E[1] <- exp(H[1]) / S[1]
      L[1] <- 1
      for(i in 2:ncol(XF)) {
        S[i] <- as.matrix(specnumber(t(cbind(matrix(rowSums(XF[1:i]),
                                                    dimnames = list(c(rownames(XF)), paste('AC', i, sep=''))),
                                             XF[-1:-i]))))[1]
        H[i] <- as.matrix(diversity(t(cbind(matrix(rowSums(XF[1:i]),
                                                   dimnames = list(c(rownames(XF)), paste('AC', i, sep=''))),
                                            XF[-1:-i])), 'shannon'))[1]
        E[i] <- exp(H[i]) / S[i]
        L[i] <- (i+1)  
      } 
      
      SHE@bi <- as.data.frame(cbind(S, H, E, L))
      colnames(SHE@bi) <- c('S', 'H', 'E', 'L')
      rownames(SHE@bi) <- colnames(df)
      
      return(SHE)
    }
  }
