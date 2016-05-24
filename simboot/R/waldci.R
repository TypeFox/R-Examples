waldci <- function(cmat, estp, varp, varcor, alternative = "two.sided", conf.level = 0.95)
  {
    k=ncol(cmat)
    m=nrow(cmat)
    if(any( c(length(estp), length(varp), length(varcor) )!=k ))
      { stop("lengths of estp, nadj, varp, varcor, and ncol(cmat) must be the same")}
    estC <- cmat %*% estp
    CorrMat<-corrmatgen(CM=cmat, varp=varcor)
    varC <- (cmat^2) %*% (varp) 
    switch(alternative,
           "two.sided" =
           { quanti <- qmvnorm(p=conf.level, sigma=CorrMat, tail="both.tails")$quantile
             stderr <- sqrt(varC)
             lCI <- estC - quanti * stderr
             uCI <- estC + quanti * stderr
           },
           "less" =
           { quanti <- qmvnorm(p=conf.level, sigma=CorrMat, tail="lower.tail")$quantile
             stderr <- sqrt(varC)
             lCI <- rep(-Inf, m)
             uCI <- estC + quanti * stderr
           },
           "greater" =
           { quanti <- qmvnorm(p=conf.level, sigma=CorrMat, tail="upper.tail")$quantile
             stderr <- sqrt(varC)
             lCI <- estC + quanti * stderr
             uCI <- rep(Inf, m)
           })
    conf.int <- cbind(estC, lCI, uCI)
    colnames(conf.int) <- c("estimate", "lower", "upper")
    return(list(conf.int = conf.int, conf.level = conf.level, alternative = alternative))
  }

