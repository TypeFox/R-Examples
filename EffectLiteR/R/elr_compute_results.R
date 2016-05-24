
computeResults <- function(obj){
  
  sem.call <- call("sem", model=obj@lavaansyntax@model,
                   group="cell", missing=obj@input@missing,
                   se=obj@input@se, bootstrap=obj@input@bootstrap,
                   mimic=obj@input@mimic,
                   group.label=obj@input@vlevels$cell, data=obj@input@data,
                   fixed.x=obj@input@fixed.z, group.w.free = !obj@input@fixed.cell)
  
  m1 <- eval(sem.call)
  
  ## lavaan.survey -- complex survey designs
  ids <- obj@input@complexsurvey$ids
  weights <- obj@input@complexsurvey$weights
  
  if((ids != ~0) | (!is.null(weights))){
    
    if(!obj@input@fixed.cell){## currently only works for fixed cell sizes
      stop("EffectLiteR error: The complex survey functionality currently only works for fixed cell sizes. Please use fixed.cell=TRUE.")
    } 
    survey.design <- survey::svydesign(ids=ids, weights=weights, 
                                       data=obj@input@data)
    m1 <- lavaan.survey::lavaan.survey(lavaan.fit=m1, 
                                       survey.design=survey.design)    
  }
  
  
  est <- parameterEstimates(m1)$est ## parameter estimates
  se <- parameterEstimates(m1)$se ## standard errors
  names(se) <- names(est) <- parameterEstimates(m1)$label 
  ## Yves: what would be the best way to get (default) parameter names 
  tval <- est/se
  pval <- 2*(1-pnorm(abs(tval)))
  
  ng <- obj@input@ng
  nz <- obj@input@nz
  nk <- obj@input@nk  
  
  # any(partable(m1)$op %in% c("==",">","<"))
  ## main hypotheses
  if(obj@input@se != "standard" | obj@input@interactions != "all" |
     any(grepl("==", obj@input@add)) | any(grepl(">", obj@input@add)) |
     any(grepl("<", obj@input@add))){ 
    ## no Wald Test for robust, bootstrapped SE...
    ## no Wald Test for models with equality constraints (ask Yves to adjust...)
    ## maybe we could come up with something similar
    hypotheses <- data.frame()
  }else{
    if(nz==0 & nk==1){
      hypotheses <- data.frame(
        lavTestWald(m1, constraints=obj@lavaansyntax@hypotheses$hypothesis1)[1:3])    
      row.names(hypotheses) <- "No average effects"    
    }else{
      hypotheses <- data.frame(rbind(
        lavTestWald(m1, constraints=obj@lavaansyntax@hypotheses$hypothesis1)[1:3],
        lavTestWald(m1, constraints=obj@lavaansyntax@hypotheses$hypothesis2)[1:3],
        lavTestWald(m1, constraints=obj@lavaansyntax@hypotheses$hypothesis3)[1:3],
        lavTestWald(m1, constraints=obj@lavaansyntax@hypotheses$hypothesis4)[1:3]    
      ))
      row.names(hypotheses) <- c("No average effects",
                                 "No covariate effects in control group",
                                 "No treatment*covariate interaction",
                                 "No treatment effects")    
    }    
  }
  
  
  
  ## average total effects
  sdyx0 <- NA
  
  mm <- obj@input@measurement
  Pkgx <- est[obj@parnames@Pkgx][1:nk]
  vary <- meany <- numeric(nk)  
  
  ## manifest outcome variable  
  if(obj@input@vnames$y %in% names(obj@input@data)){
    
    fv <- fitted.values(m1)[1:nk] ## means and variances given X=0, K=k
    
    for(i in 1:nk){
      tmp <- fv[[i]]
      vary[i] <- tmp$cov[[obj@input@vnames$y,obj@input@vnames$y]]
      meany[i] <- tmp$mean[obj@input@vnames$y]
    }
    
  }else{ ## latent outcome variable
    
    fv.cov <- inspect(m1, what="cov.lv")[1:nk]
    fv.mean <- inspect(m1, what="mean.lv")[1:nk]    
    
    for(i in 1:nk){
      tmp.cov <- fv.cov[[i]]
      tmp.mean <- fv.mean[[i]]
      
      vary[i] <- tmp.cov[[obj@input@vnames$y,obj@input@vnames$y]]
      meany[i] <- tmp.mean[[obj@input@vnames$y]]
    }
    
  }
  
  meanyx0 <- sum(meany*Pkgx) ##TODO: compute parameter in model constraint
  sdyx0 <- sqrt(sum(vary*Pkgx) + sum(Pkgx*(meany-meanyx0)^2)) ## law of total variance  
  
  Egx <- data.frame(est[obj@parnames@Egx],
                    se[obj@parnames@Egx],
                    tval[obj@parnames@Egx],
                    pval[obj@parnames@Egx],
                    est[obj@parnames@Egx]/sdyx0)
  
  names(Egx) <- c("Estimate", "SE", "Est./SE", "p-value", "Effect Size")
  
  ## Effects given a treatment condition
  Egxgx <- data.frame(est[obj@parnames@Egxgx],
                      se[obj@parnames@Egxgx],
                      tval[obj@parnames@Egxgx],
                      pval[obj@parnames@Egxgx],
                      est[obj@parnames@Egxgx]/sdyx0)
  names(Egxgx) <- c("Estimate", "SE", "Est./SE", "p-value", "Effect Size")
  
  
  ## Effects given a a value K=k
  Egxgk <- data.frame()
  if(nk>1){
    Egxgk <- data.frame(est[obj@parnames@Egxgk],
                        se[obj@parnames@Egxgk],
                        tval[obj@parnames@Egxgk],
                        pval[obj@parnames@Egxgk],
                        est[obj@parnames@Egxgk]/sdyx0)
    names(Egxgk) <- c("Estimate", "SE", "Est./SE", "p-value", "Effect Size")    
  }
  
  ## Effects given X=x and K=k
  Egxgxk <- data.frame()
  if(nk>1 & nz>0){
    Egxgxk <- data.frame(est[obj@parnames@Egxgxk],
                         se[obj@parnames@Egxgxk],
                         tval[obj@parnames@Egxgxk],
                         pval[obj@parnames@Egxgxk],
                         est[obj@parnames@Egxgxk]/sdyx0)
    names(Egxgxk) <- c("Estimate", "SE", "Est./SE", "p-value", "Effect Size")    
  }
  
  
  ## g functions
  gammas <- matrix(obj@parnames@gammas, ncol=ng)
  gammalabels <- matrix(obj@parnames@gammalabels, ncol=ng)
  gx <- vector("list",ng)
  
  for(i in 1:ng){
    tmp <- data.frame(gammas[,i],
                      est[gammas[,i]],
                      se[gammas[,i]],
                      tval[gammas[,i]],
                      pval[gammas[,i]])
    names(tmp) <- c("Coefficient", "Estimate", "SE", "Est./SE", "p-value")
    rownames(tmp) <- gammalabels[,i]
    gx[[i]] <- tmp 
  }
  
  ## adjusted means
  adjmeans <- data.frame(est[obj@parnames@adjmeans],
                         se[obj@parnames@adjmeans],
                         tval[obj@parnames@adjmeans])
  names(adjmeans) <- c("Estimate", "SE", "Est./SE")
  
  ## conditional effects
  vcov <- lavInspect(m1, "vcov.def", add.class = FALSE)
  condeffects <- computeConditionalEffects(obj, est, vcov, m1)
  
  res <- new("results",
             lavresults=m1,
             hypotheses=hypotheses,
             Egx=Egx,
             Egxgx=Egxgx,
             Egxgk=Egxgk,
             Egxgxk=Egxgxk,
             gx=gx,
             adjmeans=adjmeans,
             condeffects=condeffects
  )
  
  return(res)
}

