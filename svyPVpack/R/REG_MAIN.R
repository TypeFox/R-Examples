########################################################################
# REG MAIN FUNCTIONS
########################################################################
svyPVglm <- function(formula, 
                     design, 
                     placeholder = 1:10, 
                     family = gaussian()){
  if(length(grep("\\.\\.",formula)) == 0){stop("NO PVS")}
  else{
    as.char.forms <- lapply(placeholder,
                            function(x) gsub("\\.\\.", x, formula))
    char.calls <-  lapply(as.char.forms, 
                          function(x)(paste(x[2],x[1], x[3])))
    formulas <- lapply(char.calls, as.formula)
      mods <- lapply(lapply(formulas,as.formula), 
                     glm.dummy.no.std, 
                     design, family)
      COEF.TAB <- lapply(mods, function(x)(x$coef))
      C.est <- do.call(rbind, 
                       lapply(COEF.TAB, 
                              function(x)
                                (t(as.data.frame(x)[1]))
                       )
      )
      C.sd <- do.call(rbind, 
                      lapply(COEF.TAB, 
                             function(x)
                               (t(as.data.frame(x)[2]))
                      )
      )
      C.var <- C.sd^2
      Means <- apply(C.est, 2, mean)
      sampling.v <- apply(C.var, 2, mean)
      imputation.v <- apply(C.est,2, var) 
      tot.se <- sqrt(sampling.v + (1.1 * imputation.v))
      t.value <- Means/tot.se
      t.value.4.test <- abs(t.value)
      Pr.t <- pt(t.value.4.test , design$degf, lower.tail = FALSE)*2
      
      COEF.wo.Pr.t <- round(rbind(Means, tot.se, t.value), 
                            digits = getOption("digits"))
      COEF <- rbind(COEF.wo.Pr.t, Pr.t)
      rownames(COEF) <- c("mean", "se", "t.value", "Pr.t") 
      FIT.TAB <- lapply(mods, function(x)(x$mod.fit))
      FIT.est <- t(do.call(cbind, 
                           lapply(FIT.TAB, 
                                  function(x)
                                    (t(as.data.frame(x)))
                           )
      ))
      FIT <- t(as.data.frame(apply(FIT.est, 2, mean)))
      Working.2.logLR <- FIT[1]/mean(FIT[2:(length(FIT)-1)])
      p.dummy <- pFsum(FIT[1],
                       rep(1,length(FIT[2:(length(FIT)-1)])),
                       FIT[2:(length(FIT)-1)],ddf=FIT[length(FIT)],
                       method="sad",
                       lower.tail=FALSE)
      p <- format.pval(p.dummy)
      FIT.return <- data.frame(Working.2.logLR, p)
      rownames(FIT.return) <- ""    
      t.COEF <- as.data.frame(t(COEF))
      t.COEF$Pr.t <- format.pval(t.COEF$Pr.t)
      FINAL.RES <- list(coef=t.COEF, mod.fit=FIT.return)
      return(FINAL.RES)}
}
