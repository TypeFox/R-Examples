compareModels <- 
  function (fm0, fm1, scaled = FALSE, fm0.scaling = 1, fm1.scaling = 1) {
    if(!scaled) {
      delta.chi <- as.numeric(fm0[1] - fm1[1])
      delta.df <- as.numeric(fm0[2] - fm1[2])
      delta.p.value <- (1 - stats::pchisq(delta.chi, delta.df))
      delta.cfi <- (fm1[3] - fm0[3])
      delta.rmsea <- (fm1[4] - fm0[4])
      delta.bic <- (fm1[5] - fm0[5])
      result <- c(delta.chi, delta.df, delta.p.value, delta.cfi, delta.rmsea, delta.bic)
      names(result) <- c("delta.chisq", "delta.df", "delta.p.value", "delta.cfi", "delta.rmsea", "delta.bic")
      result  
    
    } else {
      
      delta.df <- as.numeric(fm0[2] - fm1[2])
      naive.chi <- as.numeric(fm0[1] - fm1[1])
      delta.rmsea <- (fm1[4] - fm0[4])
      delta.bic <- (fm1[5] - fm0[5])
      
      # use formula from mplus web note (www.statmodel.com)
      cd <- as.numeric((fm0[2] * fm0.scaling - fm1[2] * fm1.scaling)/delta.df)
      delta.chi <- naive.chi/cd
      delta.p.value <- (1 - stats::pchisq(delta.chi, delta.df))
      delta.cfi <- (fm1[3] - fm0[3])
      result <- c(delta.chi, delta.df, delta.p.value, delta.cfi, delta.rmsea, delta.bic)
      names(result) <- c("delta.chisq.scaled", "delta.df.scaled", "p.value.scaled", "delta.cfi.scaled", "delta.rmsea.scaled", "delta.bic")
      result
    }
    
}