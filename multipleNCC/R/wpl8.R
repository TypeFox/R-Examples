################################################################################
##########################  WPL  ###############################################
################################################################################


##Depend on survival, mgcv

wplEst = function(x,data, brukt, m, weight.method, no.intervals,
                  variance, no.intervals.left, match.var, match.int)  {
  
  if(dim(x[[2]])[2] == 3)  {
    left.time = x[[2]][,1]
    survtime = x[[2]][,2]
  }
  
  if(dim(x[[2]])[2] == 2)  {
    left.time = 0
    survtime = x[[2]][,1]
  }
  
  x = x[[1]]
  
  ##No left truncation
  if(sum(left.time) == 0)  {
    left.time = 0
  }
  
  status = brukt-1
  status[status<0] = 0
  samplestat = brukt
  samplestat[brukt!=0] = 1
  
  
  forventet.n.cohort = length(samplestat)
  
  if(sum(left.time>=survtime) > 0)  {
    stop("Stop time must be > start time")
  }
  
  if(any(is.na(match.int)))  {
    stop("The matching interval cannot include NA")
  }
  
  if(variance == "Modelbased" & weight.method != "KM")  {
    stop("Estimated variance (Modelbased) can only be applied together with KM-weighs")
  }
  
  if(variance == "Poststrat" & weight.method != "Chen")  {
    stop("Estimated variance (Poststrat) can only be applied together with Chen-weighs")
  }
  
  if(length(survtime) < forventet.n.cohort | length(status) < forventet.n.cohort)  {
    stop("The number of elements in survtime, status, left.time or covariate matrix does not match the number of subjects in the cohort. It might be due to NA-values, these should be removed.")
  }
  
  if(any(is.na(match.var)))  {
    stop("The matching variable(s) cannot be NA")
  }
  
  if(any(is.na(samplestat)))  {
    stop("samplestat cannot be NA")
  }
  
  if(any(is.na(m)))  {
    stop("m cannot be NA")
  }
  
  if(length(m) > 1 & length(m)!=sum(status!=0))  {
    stop("m must either be a scalar or a vector of the same length as number of cases")
  }
  
  if(match.int!=0 && length(match.int) < dim(match.var)[2]*2)  {
    stop("Too few elements in match.int")
  }
  
  if(match.int != 0 && length(match.int) > dim(match.var)[2]*2)  {
    stop("Too many elements in match.int")
  }
  
  
  #cat("Note that the data are sorted by survtime")
  #cat( "\n")
  
  
  endpoints = length(unique(status))-1-1*any(is.na(status))
  ind.no = 1:forventet.n.cohort
  
  n.cohort = length(survtime)
  x = as.matrix(x[which(samplestat==1),])
  n = dim(x)[1]
  
  
  ##Sorting vectors by survtime
  o = order(survtime)
  oo = order(survtime[samplestat!=0])
  ooo = order(survtime[status > 0])
  status = status[o]
  samplestat = samplestat[o]
  x = x[oo,]
  if(length(m) > 1)  {
    m = m[ooo]
  }
  
  samplestat = as.numeric(samplestat!=0)
  
  survtime = survtime[o]
  if(length(left.time)==n.cohort)  {
    left.time = left.time[o]
  }
  
  if(length(match.var)==n.cohort)  {
    match.var = match.var[o]
  }
  
  if(length(match.var) > n.cohort)  {
    match.var = match.var[o,]
  }
  
  ind.no.ncc = ind.no[samplestat==1]
  status.ncc = status[samplestat==1]
  survtime.ncc = survtime[samplestat==1]
  left.time.ncc = left.time[samplestat==1]
  
  if(weight.method != "KM" & weight.method!= "gam" & weight.method!= "glm" &
       weight.method!= "Chen")  {
    stop("Invalid weight method. It should be either KM, gam, glm or Chen.")
  }
  
  if(weight.method == "KM")  {
    p = pKM(status,survtime,m,n.cohort,left.time,match.var,match.int)[samplestat==1]
    p[status[samplestat==1]!=0] = 1
    if(sum(p==0) > 0)  {
      who = which(p==0)
      cat("Warning: Some sampling probabilities are estimated as zero:", "[",who,"]" , "\n")
      cat("Sampled controls with zero sampling probability is given weight=1","\n")
      cat("\n")
      
      p[which(p==0)] = 1
    }
  }
  
  if(weight.method == "gam")  {
    pp = rep(1,n.cohort)
    p = pGAM(status,survtime,samplestat,n.cohort,left.time,match.var,match.int)
    pp[status==0] = p
    p = pp[samplestat==1]
    
    if(sum(p==0) > 0)  {
      who = which(p==0)
      cat("Warning: Some sampling probabilities are estimated as zero:", "[",who,"]" , "\n")
      cat("Sampled controls with zero sampling probability is given weight=1","\n")
      cat("\n")
      
      p[which(p==0)] = 1
    }
  }
  
  if(weight.method == "glm")  {
    pp = rep(1,n.cohort)
    
    p = pGLM(status,survtime,samplestat,n.cohort,left.time,match.var,match.int)
    pp[status==0] = p
    p = pp[samplestat==1]
    
    if(sum(p==0) > 0)  {
      who = which(p==0)
      cat("Warning: Some sampling probabilities are estimated as zero:", "[",who,"]" , "\n")
      cat("Sampled controls with zero sampling probability is given weight=1","\n")
      cat("\n")
      
      p[which(p==0)] = 1
    }
  }
  
  if(weight.method == "Chen")  {
    pp = pChen(status,survtime,samplestat,ind.no,n.cohort,no.intervals, left.time,
               no.intervals.left)
    pp[status!=0] = 1
    p = pp[samplestat==1]
    
    if(sum(p==0) > 0)  {
      who = which(p==0)
      cat("Warning: Some sampling probabilities are estimated as zero:", "[",who,"]" , "\n")
      cat("Sampled controls with zero sampling probability is given weight=1","\n")
      cat("\n")
      
      p[which(p==0)] = 1
    }
  }
  
  res = list()
  if(length(left.time)==1)  {
    strat = which(substr(colnames(x),1,7)=="strata(")
    if(length(strat) != 0) {
      xx = x
      x = xx[,-strat]
      stratum = xx[,strat]
      temp2 = stratum*1:dim(stratum)[2]
      stratum = (apply(temp2,1,sum))
    }
    for(i in 1:(endpoints))  {
      if(length(strat)==0){
        fit = coxph(Surv(survtime.ncc,status.ncc==i)~x+cluster(ind.no.ncc),weights=1/p)
      }
      if(length(strat)>0)  {
        fit = coxph(Surv(survtime.ncc,status.ncc==i)~x+strata(stratum)+cluster(ind.no.ncc),weights=1/p)
        fit = fit[-21]
      }
      fit$est.var=F
      if(variance == "Modelbased")   {
        left = rep(0,n.cohort)
        fit$var = ModelbasedVar(fit,status,survtime,left,samplestat,i,m,p,match.var,match.int)
        fit$est.var = T
      }
      if(variance == "Poststrat")  {
        left = rep(0,n)
        fit$var = PoststratVar(fit,status,samplestat,survtime,left,i,m,no.intervals)
        fit$est.var=T
      }
      names(fit)[3] = "weighted.loglik"
      fit = fit[-c(4,9)]
      res = c(res,fit)
    }
  }
  
  if(length(left.time)==n.cohort)  {
    strat = which(substr(colnames(x),1,7)=="strata(")
    if(length(strat) != 0) {
      xx = x
      x = xx[,-strat]
      stratum = xx[,strat]
      temp2 = stratum*1:dim(stratum)[2]
      stratum = (apply(temp2,1,sum))
    }
    
    for(i in 1:(endpoints))  {
      if(length(strat) == 0) {
        fit = coxph(Surv(left.time.ncc,survtime.ncc,status.ncc==i)~x+cluster(ind.no.ncc),weights=1/p)
      }
      if(length(strat)>0) {
        fit = coxph(Surv(left.time.ncc,survtime.ncc,status.ncc==i)~x+strata(stratum)+cluster(ind.no.ncc),weights=1/p)
        fit = fit[-21]
      }
      fit$est.var=F
      if(variance == "Modelbased")   {
        fit$var = ModelbasedVar(fit,status,survtime,left.time,samplestat,i,m,p,match.var,match.int)
        fit$est.var=T
      }
      if(variance == "Poststrat")  {
        fit$var = PoststratVar(fit,status,samplestat,survtime,left.time,i,m,no.intervals.left)
        fit$est.var=T
      }
      names(fit)[3] = "weighted.loglik"
      fit = fit[-c(4,9)]
      res = c(res,fit)
    }
  }
  res
  
}


pKM = function(status,survtime,m,n.cohort,left.time, match.var, match.int)  {
  
  failuretimes = survtime[which(status != 0 )]
  if(length(m) == 1)  {
    m = rep(m,length(failuretimes))
  }
  
  ##Without left truncation, without matching
  if(length(left.time)==1 & length(match.var)==1)  {
    #Number at risk at each event time
    nfail = rep(0,length(failuretimes))
    psample = rep(0,n.cohort)
    qsample = rep(0,length(failuretimes))
    
    for(k in 1:length(failuretimes))  {
      nfail[k] =  length(survtime[which(survtime >= failuretimes[k])])
      if (nfail[k] >= m[k])  {
        qsample[k] = (1-m[k]/(nfail[k]-1))
      }
      if(nfail[k] < m[k])  {
        cat("Warning: Fewer subjects at risk than number of sampled controls for case ", k , "\n")
        cat("Controls for case ", k , "are given weight=1", "\n")
        cat("\n")
      }
    }
    
    for(k in 1:n.cohort) {
      indeks = 1:length(failuretimes)
      indeks = indeks[which(failuretimes <= survtime[k])]
      psample[k] = 1-prod(qsample[indeks])
    }
  }
  
  ##With left truncation, without matching
  if(length(left.time) == n.cohort & length(match.var)==1)  {
    ##Number at riks at each event time
    nfail = rep(0,length(failuretimes))
    psample = rep(0,n.cohort)
    qsample = rep(0,sum(status!=0))
    for(k in 1:length(failuretimes))  {
      nfail[k] =  length(survtime[which(survtime >= failuretimes[k] & left.time
                                        < failuretimes[k])])
      
      if (nfail[k] >= m[k])  {
        qsample[k] = (1-m[k]/(nfail[k]-1))
      }
      if(nfail[k] < m[k])  {
        cat("Warning: Fewer subjects at risk than number of sampled controls for case ", k , "\n")
        cat("Controls for case ", k , "are given weight=1", "\n")
        cat("\n")
      }
    }
    
    for(k in 1:n.cohort) {
      indeks = left.time[k] < survtime[status!=0] & survtime[k] >
        survtime[status!=0]
      psample[k] = 1-prod(qsample[indeks])
    }
  }
  
  #Without left truncation, with matching
  
  if(length(left.time) == 1 & length(match.var) > 1)  {
    failuretimes = survtime[which(status != 0 )]
    match.var = as.matrix(match.var)
    fail.match = as.matrix(match.var[which(status != 0 ),])
    
    #One matching variable
    if(dim(match.var)[2] == 1 | length(match.var) == n.cohort)  {
      
      #Number at risk at each event time
      nfail = rep(0,length(failuretimes))
      psample = rep(0,n.cohort)
      qsample = rep(0,length(failuretimes))
      
      
      for(k in 1:length(failuretimes))  {
        nfail[k] = length(survtime[which(survtime >= failuretimes[k] &
                                           match.var[,1] >= fail.match[k,1]+match.int[1] & match.var[,1] <= fail.match[k,1]+match.int[2])])
        if (nfail[k] >= m[k])  {
          qsample[k] = (1-m[k]/(nfail[k]-1))
        }
        if(nfail[k] < m[k])  {
          cat("Warning: Fewer subjects at risk than number of sampled controls for case ", k , "\n")
          cat("Controls for case ", k , "are given weight=1", "\n")
          cat("\n")
        }
      }
      
      for(k in 1:n.cohort) {
        ind = survtime[k] > failuretimes  &
          match.var[k,1] >= fail.match[,1]+match.int[1] & match.var[k,1] <= fail.match[,1]+match.int[2]
        psample[k] = 1-prod(qsample[ind])
      }
    }
    
    #Two matching variables
    if(dim(match.var)[2] == 2)  {
      
      #Number at risk at each event time
      nfail = rep(0,length(failuretimes))
      psample = rep(0,n.cohort)
      qsample = rep(0,length(failuretimes))
      
      
      for(k in 1:length(failuretimes))  {
        nfail[k] = length(survtime[which(survtime >= failuretimes[k] &
                                           match.var[,1] >= fail.match[k,1]+match.int[1] & match.var[,1] <= fail.match[k,1]+match.int[2] &
                                           match.var[,2] >= fail.match[k,2]+match.int[3] & match.var[,2] <= fail.match[k,2]+match.int[4])])
        if (nfail[k] >= m[k])  {
          qsample[k] = (1-m[k]/(nfail[k]-1))
        }
        if(nfail[k] < m[k])  {
          cat("Warning: Fewer subjects at risk than number of sampled controls for case ", k , "\n")
          cat("Controls for case ", k , "are given weight=1", "\n")
          cat("\n")
        }
      }
      
      for(k in 1:n.cohort) {
        ind = survtime[k] > failuretimes  &
          match.var[k,1] >= fail.match[,1]+match.int[1] & match.var[k,1] <= fail.match[,1]+match.int[2] &
          match.var[k,2] >= fail.match[,2]+match.int[3] & match.var[k,2] <= fail.match[,2]+match.int[4]
        psample[k] = 1-prod(qsample[ind])
      }
    }
    
    #Three matching variables
    if(dim(match.var)[2] == 3)  {
      
      #Number at risk at each event time
      nfail = rep(0,length(failuretimes))
      psample = rep(0,n.cohort)
      qsample = rep(0,length(failuretimes))
      
      
      for(k in 1:length(failuretimes))  {
        nfail[k] = length(survtime[which(survtime >= failuretimes[k] &
                                           match.var[,1] >= fail.match[k,1]+match.int[1] & match.var[,1] <= fail.match[k,1]+match.int[2] &
                                           match.var[,2] >= fail.match[k,2]+match.int[3] & match.var[,2] <= fail.match[k,2]+match.int[4] &
                                           match.var[,3] >= fail.match[k,3]+match.int[5] & match.var[,3] <= fail.match[k,3]+match.int[6])])
        if (nfail[k] >= m[k])  {
          qsample[k] = (1-m[k]/(nfail[k]-1))
        }
        if(nfail[k] < m[k])  {
          cat("Warning: Fewer subjects at risk than number of sampled controls for case ", k , "\n")
          cat("Controls for case ", k , "are given weight=1", "\n")
          cat("\n")
        }
      }
      
      for(k in 1:n.cohort) {
        ind = survtime[k] > failuretimes &
          match.var[k,1] >= fail.match[,1]+match.int[1] & match.var[k,1] <= fail.match[,1]+match.int[2] &
          match.var[k,2] >= fail.match[,2]+match.int[3] & match.var[k,2] <= fail.match[,2]+match.int[4] &
          match.var[k,3] >= fail.match[,3]+match.int[5] & match.var[k,3] <= fail.match[,3]+match.int[6]
        psample[k] = 1-prod(qsample[ind])
      }
    }
    
    #More than three matching variables
    if(dim(match.var)[2] > 3)  {
      tt = function(rad)  {
        sum(rad==1)==length(rad)
      }
      
      ncont1 = function(M,T)  {
        mat = matrix(ncol = length(M)+1,nrow=dim(match.var)[1],0)
        mat[which(survtime >= T),1] = 1
        for(i in 1:length(M))  {
          mat[which(match.var[,i] >= M[i]+match.int[(2*i-1)] & match.var[,i] <= M[i]+match.int[(2*i)]),(i+1)] = 1
        }
        sum(apply(mat,1,tt))
      }
      #Number at risk at each event time
      nfail = rep(0,length(failuretimes))
      psample = rep(0,n.cohort)
      qsample = rep(0,length(failuretimes))
      
      for(k in 1:length(failuretimes))  {
        nfail[k] =  ncont1(fail.match[k,],failuretimes[k])
        if (nfail[k] >= m[k])  {
          qsample[k] = (1-m[k]/(nfail[k]-1))
        }
        if(nfail[k] < m[k])  {
          cat("Warning: Fewer subjects at risk than number of sampled controls for case ", k , "\n")
          cat("Controls for case ", k , "are given weight=1", "\n")
          cat("\n")
        }
      }
      
      indeks1 = function(M,T)  {
        mat = matrix(ncol = length(M)+1,nrow=dim(fail.match)[1],0)
        mat[which(failuretimes <T),1] = 1
        for(i in 1:length(M))  {
          mat[which(M[i] >= fail.match[,i]+match.int[(2*i-1)] & M[i] <= fail.match[,i]+match.int[(2*i)]),(i+1)] = 1
        }
        apply(mat,1,tt)
      }
      
      for(k in 1:n.cohort) {
        ind = indeks1(match.var[k,],survtime[k])
        psample[k] = 1-prod(qsample[ind])
      }
    }
  }
  
  #With left truncation, with matching
  if(length(left.time) == n.cohort & length(match.var) > 1)  {
    failuretimes = survtime[which(status != 0 )]
    match.var = as.matrix(match.var)
    fail.match = as.matrix(match.var[which(status != 0 ),])
    
    #One matching variabel
    if(dim(match.var)[2] == 1 | length(match.var) == n.cohort)  {
      
      #Number at risk at each event time
      nfail = rep(0,length(failuretimes))
      psample = rep(0,n.cohort)
      qsample = rep(0,length(failuretimes))
      
      for(k in 1:length(failuretimes))  {
        nfail[k] = length(survtime[which(survtime >= failuretimes[k] & left.time < failuretimes[k] &
                                           match.var[,1] >= fail.match[k,1]+match.int[1] & match.var[,1] <= fail.match[k,1]+match.int[2])])
        if (nfail[k] >= m[k])  {
          qsample[k] = (1-m[k]/(nfail[k]-1))
        }
        if(nfail[k] < m[k])  {
          cat("Warning: Fewer subjects at risk than number of sampled controls for case ", k , "\n")
          cat("Controls for case ", k , "are given weight=1", "\n")
          cat("\n")
        }
      }
      
      for(k in 1:n.cohort) {
        ind = survtime[k] > failuretimes & left.time[k] < failuretimes &
          match.var[k,1] >= fail.match[,1]+match.int[1] & match.var[k,1] <= fail.match[,1]+match.int[2]
        psample[k] = 1-prod(qsample[ind])
      }
    }
    
    #Two matching variables
    if(dim(match.var)[2] == 2)  {
      
      #Number at risk at each event time
      nfail = rep(0,length(failuretimes))
      psample = rep(0,n.cohort)
      qsample = rep(0,length(failuretimes))
      
      for(k in 1:length(failuretimes))  {
        nfail[k] = length(survtime[which(survtime >= failuretimes[k] & left.time < failuretimes[k] &
                                           match.var[,1] >= fail.match[k,1]+match.int[1] & match.var[,1] <= fail.match[k,1]+match.int[2] &
                                           match.var[,2] >= fail.match[k,2]+match.int[3] & match.var[,2] <= fail.match[k,2]+match.int[4])])
        if (nfail[k] >= m[k])  {
          qsample[k] = (1-m[k]/(nfail[k]-1))
        }
        if(nfail[k] < m[k])  {
          cat("Warning: Fewer subjects at risk than number of sampled controls for case ", k , "\n")
          cat("Controls for case ", k , "are given weight=1", "\n")
          cat("\n")
        }
      }
      
      for(k in 1:n.cohort) {
        ind = survtime[k] > failuretimes & left.time[k] < failuretimes &
          match.var[k,1] >= fail.match[,1]+match.int[1] & match.var[k,1] <= fail.match[,1]+match.int[2] &
          match.var[k,2] >= fail.match[,2]+match.int[3] & match.var[k,2] <= fail.match[,2]+match.int[4]
        psample[k] = 1-prod(qsample[ind])
      }
    }
    
    #Three matching variables
    if(dim(match.var)[2] == 3)  {
      
      ##Number at risk at each event time
      nfail = rep(0,length(failuretimes))
      psample = rep(0,n.cohort)
      qsample = rep(0,length(failuretimes))
      
      for(k in 1:length(failuretimes))  {
        nfail[k] = length(survtime[which(survtime >= failuretimes[k] & left.time < failuretimes[k] &
                                           match.var[,1] >= fail.match[k,1]+match.int[1] & match.var[,1] <= fail.match[k,1]+match.int[2] &
                                           match.var[,2] >= fail.match[k,2]+match.int[3] & match.var[,2] <= fail.match[k,2]+match.int[4] &
                                           match.var[,3] >= fail.match[k,3]+match.int[5] & match.var[,3] <= fail.match[k,3]+match.int[6])])
        if (nfail[k] >= m[k])  {
          qsample[k] = (1-m[k]/(nfail[k]-1))
        }
        if(nfail[k] < m[k])  {
          cat("Warning: Fewer subjects at risk than number of sampled controls for case ", k , "\n")
          cat("Controls for case ", k , "are given weight=1", "\n")
          cat("\n")
        }
      }
      
      for(k in 1:n.cohort) {
        ind = survtime[k] > failuretimes & left.time[k] < failuretimes &
          match.var[k,1] >= fail.match[,1]+match.int[1] & match.var[k,1] <= fail.match[,1]+match.int[2] &
          match.var[k,2] >= fail.match[,2]+match.int[3] & match.var[k,2] <= fail.match[,2]+match.int[4] &
          match.var[k,3] >= fail.match[,3]+match.int[5] & match.var[k,3] <= fail.match[,3]+match.int[6]
        psample[k] = 1-prod(qsample[ind])
      }
    }
    
    #More than three matching variables
    if(dim(match.var)[2] > 3)  {
      
      tt = function(rad)  {
        sum(rad==1)==length(rad)
      }
      
      ncont2 = function(M,T,V)  {
        mat = matrix(ncol = length(M)+2,nrow=dim(match.var)[1],0)
        mat[which(survtime >= T),1] = 1
        mat[which(left.time < T),2] = 1
        for(i in 1:length(M))  {
          mat[which(match.var[,i] >= M[i]+match.int[(2*i-1)] & match.var[,i] <= M[i]+match.int[(2*i)]),(i+2)] = 1
        }
        sum(apply(mat,1,tt))
      }
      #Number at risk at each event time
      nfail = rep(0,length(failuretimes))
      psample = rep(0,n.cohort)
      qsample = rep(0,length(failuretimes))
      
      for(k in 1:length(failuretimes))  {
        nfail[k] =  ncont2(fail.match[k,],failuretimes[k],survtime[status!=0][k])
        if (nfail[k] >= m[k])  {
          qsample[k] = (1-m[k]/(nfail[k]-1))
        }
        if(nfail[k] < m[k])  {
          cat("Warning: Fewer subjects at risk than number of sampled controls for case ", k , "\n")
          cat("Controls for case ", k , "are given weight=1", "\n")
          cat("\n")
        }
      }
      
      indeks2 = function(M,T,V)  {
        mat = matrix(ncol = length(M)+2,nrow=dim(fail.match)[1],0)
        mat[which(failuretimes <T),1] = 1
        mat[which(failuretimes > V),2] = 1
        for(i in 1:length(M))  {
          mat[which(M[i] >= fail.match[,i]+match.int[(2*i-1)] & M[i] <= fail.match[,i]+match.int[(2*i)]),(i+2)] = 1
        }
        apply(mat,1,tt)
      }
      
      for(k in 1:n.cohort) {
        ind = indeks2(match.var[k,],survtime[k],left.time[k])
        psample[k] = 1-prod(qsample[ind])
      }
    }
  }
  psample
}

pGAM = function(status,survtime,samplestat,n.cohort,left.time,match.var,match.int)  {
  
  if(length(match.int)>1)  {
    kat = 0 
    kont = 0
    temp1 = which(match.int == 0)
    temp2 = which(match.int != 0)
    if(length(temp1)>1) {
      indtemp = seq(2,length(temp1),2)
      kat = temp1[indtemp]/2  
    }
    if(length(temp2)>1) {
      indtemp = seq(2,length(temp2),2)
      kont = temp2[indtemp]/2  
    }
  }
  
  if(length(left.time)==1 & length(match.var) == 1)  {
    pgam = gam(samplestat~s(survtime),family=binomial,subset=status==0)
    pgam = pgam$fitted
  }
  if(length(left.time)==n.cohort & length(match.var) == 1)  {
    pgam = gam(samplestat~s(survtime)+s(left.time),family=binomial,subset=status==0)
    pgam = pgam$fitted
  }
  
  if(length(left.time)==1 & length(match.var) > 1)  {
    if(length(kat) == 1 && kat == 0) {
      pgam = gam(samplestat~s(survtime)+s(match.var),family=binomial,subset=status==0)
    }
    else if(length(kont) == 1 && kont == 0) {
      pgam = gam(samplestat~s(survtime)+factor(match.var),family=binomial,subset=status==0)
    }
    else{
      pgam = gam(samplestat~s(survtime)+s(match.var[,kont])+factor(match.var[,kat]),family=binomial,subset=status==0)
    }
    pgam = pgam$fitted
  }
  
  if(length(left.time)==n.cohort & length(match.var) > 1)  {
    if(length(kat) == 1 && kat == 0) {
      pgam = gam(samplestat~s(survtime)+s(left.time)+s(match.var),family=binomial,subset=status==0)
    }
    else if(length(kont) == 1 && kont == 0) {
      pgam = gam(samplestat~s(survtime)+s(left.time)+factor(match.var),family=binomial,subset=status==0)
    }
    else{
      pgam = gam(samplestat~s(survtime)+s(left.time)+s(match.var[,kont])+factor(match.var[,kat]),family=binomial,subset=status==0)
    }
    pgam = pgam$fitted
  }
  pgam
}

pGLM = function(status,survtime,samplestat,n.cohort,left.time,match.var,match.int)   {
  if(length(match.int)>1)  {
    kat = 0 
    kont = 0
    temp1 = which(match.int == 0)
    temp2 = which(match.int != 0)
  
    if(length(temp1)>1) {
      indtemp = seq(2,length(temp1),2)
      kat = temp1[indtemp]/2  
    }
    if(length(temp2)>1) {
      indtemp = seq(2,length(temp2),2)
      kont = temp2[indtemp]/2  
    }
  }
  
  if(length(left.time)==1 & length(match.var) == 1)  {
    pglm = glm(samplestat~survtime,family=binomial,subset=status==0)
    pglm = pglm$fitted
  }
  if(length(left.time)==1 & length(match.var) > 1)  {
    if(kat == 0) {
      pglm = glm(samplestat~survtime+match.var,family=binomial,subset=status==0)
    }
    if(kont == 0) {
      pglm = glm(samplestat~survtime+factor(match.var),family=binomial,subset=status==0)
    }
    if(kont != 0 & kat != 0) {
      pglm = glm(samplestat~survtime+match.var[,kont]+factor(match.var[,kat]),family=binomial,subset=status==0)
    }
    pglm = pglm$fitted
  }
  
  if(length(left.time)==n.cohort & length(match.var) == 1)  {
    pglm = glm(samplestat~survtime+left.time,family=binomial,subset=status==0)
    pglm = pglm$fitted
  }
  
  if(length(left.time)==n.cohort & length(match.var) > 1)  {
    if(length(kat) == 1 && kat == 0) {
      pglm = glm(samplestat~survtime+left.time+match.var,family=binomial,subset=status==0)
    }
    else if(length(kont)==1 && kont==0) {
      pglm = glm(samplestat~survtime+left.time+factor(match.var),family=binomial,subset=status==0)
    }
    else{#(kont != 0 & kat != 0) {
      pglm = glm(samplestat~survtime+left.time+match.var[,kont]+factor(match.var[,kat]),family=binomial,subset=status==0)
    }  
    pglm = pglm$fitted
  }
  pglm
}

pChen = function(status,survtime,samplestat,ind.no,n.cohort,no.intervals,left.time,
                  no.intervals.left)  {
  
  if(length(left.time)==1)  {
    pchen = 1:no.intervals
    partT = seq(min(survtime-0.001),max(survtime),length=(no.intervals+1))
    for(i in 1:(length(partT)))   {
      ne = sum(status == 0 & survtime > partT[i] &
                 survtime <= partT[i+1])
      te = sum(status == 0 & samplestat != 0 &
                 survtime > partT[i] & survtime <= partT[i+1])
      pchen[i] = te/ne
    }
    
    #controls
    pp = rep(0,n.cohort)
    for(i in 1:(length(partT)))  {
      test = ind.no[(survtime > partT[i] & survtime <= partT[i+1]
                     & status == 0)]
      pp[test] = pchen[i]
    }
  }
  
  if(length(left.time)==n.cohort)  {
    pchen = matrix(0,nrow = (no.intervals.left[2]+1), ncol = (no.intervals.left[1]+1))
    partT = seq(min(survtime-0.001),max(survtime),length=(no.intervals.left[2]+1))
    partV = seq(min(left.time-0.001),max(left.time),length=(no.intervals.left[1]+1))
    for(i in 1:(length(partT)))   {
      for(ii in 1:(length(partV))) {
        testT = survtime > partT[i]  & survtime <= partT[i+1]
        testV = left.time > partV[ii]  & left.time <= partV[ii+1]
        ne = sum(status == 0 & testT & testV)
        te = sum(status == 0 & samplestat != 0 & testT & testV)
        pchen[i,ii] = te/ne
      }
    }
    
    #Controls
    pp = rep(0,n.cohort)
    for(i in 1:(length(partT)))  {
      for(ii in 1:(length(partV)))  {
        testT = survtime > partT[i]  & survtime <= partT[i+1]
        testV = left.time > partV[ii]  & left.time <= partV[ii+1]
        test = ind.no[(testT & testV & samplestat != 0 & status == 0)]
        pp[test] = pchen[i,ii]
      }
    }
  }
  pp
}


ModelbasedVar = function(fit,status,survtime,left.time,samplestat,endpoint,m,psample,match.var,match.int)  {
  samplestat = samplestat!=0
  renkont = sum(samplestat!=0 & status==0)
  failuretimes = survtime[which(status != 0)]
  tidcase<-survtime[status==endpoint]
  leftcase = left.time[status==endpoint]
  tidkont<-survtime[status==0 & samplestat!=0]
  leftkont = left.time[status==0 & samplestat!=0]
  
  #With matching without left truncation
  
  if(sum(left.time) == 0 & length(match.var) > 1)  {
    
    match.var = as.matrix(match.var)
    fail.match = as.matrix(match.var[which(status != 0 ),])
    kont.match = as.matrix(match.var[which(status == 0 & samplestat!=0),])
    nfail = rep(0,length(failuretimes))
    
    if(dim(match.var)[2] == 1)  {
      
      for(k in 1:length(failuretimes))  {
        nfail[k] = length(survtime[which(survtime >= failuretimes[k]&
                                           match.var[,1] >= fail.match[k,1]+match.int[1] & match.var[,1] <= fail.match[k,1]+match.int[2])])
      }
      
      H = 1-2*m/(nfail-1)+m*(m-1)/((nfail-1)*(nfail-2))
      H = H/(1-m/(nfail-1))^2
      
      Rho<-matrix(rep(0,renkont^2),ncol=renkont)
      for(i in 1:(renkont-1)) {
        for(j in ((i+1):renkont))  {
          ind = survtime[i] > failuretimes & survtime[j] > failuretimes &
            match.var[i,1] >= fail.match[,1]+match.int[1] & match.var[i,1] <= fail.match[,1]+match.int[2] &
            match.var[j,1] >= fail.match[,1]+match.int[1] & match.var[j,1] <= fail.match[,1]+match.int[2]
          Rho[i,j] = prod(H[ind])-1
        }
      }
    }
    
    
    if(dim(match.var)[2] == 2)  {
      
      for(k in 1:length(failuretimes))  {
        nfail[k] = length(survtime[which(survtime >= failuretimes[k] &
                                           match.var[,1] >= fail.match[k,1]+match.int[1] & match.var[,1] <= fail.match[k,1]+match.int[2] &
                                           match.var[,2] >= fail.match[k,2]+match.int[3] & match.var[,2] <= fail.match[k,2]+match.int[4])])
      }
      
      H = 1-2*m/(nfail-1)+m*(m-1)/((nfail-1)*(nfail-2))
      H = H/(1-m/(nfail-1))^2
      
      Rho<-matrix(rep(0,renkont^2),ncol=renkont)
      for(i in 1:(renkont-1)) {
        for(j in ((i+1):renkont))  {
          ind = survtime[i] > failuretimes & survtime[j] > failuretimes &
            match.var[i,1] >= fail.match[,1]+match.int[1] & match.var[i,1] <= fail.match[,1]+match.int[2] &
            match.var[j,1] >= fail.match[,1]+match.int[1] & match.var[j,1] <= fail.match[,1]+match.int[2] &
            match.var[i,2] >= fail.match[,2]+match.int[3] & match.var[i,2] <= fail.match[,2]+match.int[4] &
            match.var[j,2] >= fail.match[,2]+match.int[3] & match.var[j,2] <= fail.match[,2]+match.int[4]
          Rho[i,j] = prod(H[ind])-1
        }
      }
    }
    
    
    if(dim(match.var)[2] == 3)  {
      
      for(k in 1:length(failuretimes))  {
        nfail[k] = length(survtime[which(survtime >= failuretimes[k] &
                                           match.var[,1] >= fail.match[k,1]+match.int[1] & match.var[,1] <= fail.match[k,1]+match.int[2] &
                                           match.var[,2] >= fail.match[k,2]+match.int[3] & match.var[,2] <= fail.match[k,2]+match.int[4] &
                                           match.var[,3] >= fail.match[k,3]+match.int[5] & match.var[,3] <= fail.match[k,3]+match.int[6])])
      }
      
      H = 1-2*m/(nfail-1)+m*(m-1)/((nfail-1)*(nfail-2))
      H = H/(1-m/(nfail-1))^2
      
      Rho<-matrix(rep(0,renkont^2),ncol=renkont)
      for(i in 1:(renkont-1)) {
        for(j in ((i+1):renkont))  {
          ind = survtime[i] > failuretimes & survtime[j] > failuretimes &
            match.var[i,1] >= fail.match[,1]+match.int[1] & match.var[i,1] <= fail.match[,1]+match.int[2] &
            match.var[j,1] >= fail.match[,1]+match.int[1] & match.var[j,1] <= fail.match[,1]+match.int[2] &
            match.var[i,2] >= fail.match[,2]+match.int[3] & match.var[i,2] <= fail.match[,2]+match.int[4] &
            match.var[j,2] >= fail.match[,2]+match.int[3] & match.var[j,2] <= fail.match[,2]+match.int[4] &
            match.var[i,3] >= fail.match[,3]+match.int[5] & match.var[i,3] <= fail.match[,3]+match.int[6] &
            match.var[j,3] >= fail.match[,3]+match.int[5] & match.var[j,3] <= fail.match[,3]+match.int[6]
          Rho[i,j] = prod(H[ind])-1
        }
      }
    }
    
    
    if(dim(match.var)[2] > 3)  {
      
      tt = function(rad)  {
        sum(rad==1)==length(rad)
      }
      
      nfail = rep(0,length(failuretimes))
      ncont3 = function(M,T)  {
        mat = matrix(ncol = length(M)+1,nrow=dim(match.var)[1],0)
        mat[which(survtime >= T),1] = 1
        for(i in 1:length(M))  {
          mat[which(match.var[,i] >= M[i]+match.int[(2*i-1)] & match.var[,i] <= M[i]+match.int[(2*i)]),(i+1)] = 1
        }
        sum(apply(mat,1,tt))
      }
      
      for(k in 1:length(failuretimes))  {
        nfail[k] =  ncont3(fail.match[k,],failuretimes[k])
      }
      
      H = 1-2*m/(nfail-1)+m*(m-1)/((nfail-1)*(nfail-2))
      H = H/(1-m/(nfail-1))^2
      
      indeks3 = function(M1,T1,M2,T2)  {
        mat = matrix(ncol = length(M1)*2+2,nrow=dim(fail.match)[1],0)
        mat[which(tidcase <T1),1] = 1
        mat[which(tidcase <T2),2] = 1
        for(i in 1:length(M1))  {
          mat[which(M1[i] >= fail.match[,i]+match.int[(2*i-1)] & M1[i] <= fail.match[,i]+match.int[(2*i)]),(i+2)] = 1
          mat[which(M2[i] >= fail.match[,i]+match.int[(2*i-1)] & M2[i] <= fail.match[,i]+match.int[(2*i)]),(i+dim(match.var[2])+2)] = 1
        }
        apply(mat,1,tt)
      }
      
      
      Rho<-matrix(rep(0,renkont^2),ncol=renkont)
      for(i in 1:(renkont-1)) {
        for(j in ((i+1):renkont))  {
          ind = indeks3(kont.match[i,],tidkont[i],kont.match[j,],tidkont[j])
          Rho[i,j] = prod(H[ind])-1
        }
      }
    }
  }
  
  
  #With matching and left truncation
  if(sum(left.time) > 0 & length(match.var) > 1)  {
    
    match.var = as.matrix(match.var)
    fail.match = as.matrix(match.var[which(status != 0 ),])
    kont.match = as.matrix(match.var[which(status == 0 & samplestat!=0),])
    nfail = rep(0,length(failuretimes))
    if(dim(match.var)[2] == 1)  {
      
      for(k in 1:length(failuretimes))  {
        nfail[k] = length(survtime[which(survtime >= failuretimes[k] & left.time < failuretimes[k] &
                                           match.var[,1] >= fail.match[k,1]+match.int[1] & match.var[,1] <= fail.match[k,1]+match.int[2])])
      }
      
      H = 1-2*m/(nfail-1)+m*(m-1)/((nfail-1)*(nfail-2))
      H = H/(1-m/(nfail-1))^2
      
      Rho<-matrix(rep(0,renkont^2),ncol=renkont)
      for(i in 1:(renkont-1)) {
        for(j in ((i+1):renkont))  {
          ind = survtime[i] > failuretimes & survtime[j] > failuretimes & left.time[i] < failuretimes & left.time[j] < failuretimes &
            match.var[i,1] >= fail.match[,1]+match.int[1] & match.var[i,1] <= fail.match[,1]+match.int[2] &
            match.var[j,1] >= fail.match[,1]+match.int[1] & match.var[j,1] <= fail.match[,1]+match.int[2]
          Rho[i,j] = prod(H[ind])-1
        }
      }
    }
    
    
    if(dim(match.var)[2] == 2)  {
      
      for(k in 1:length(failuretimes))  {
        nfail[k] = length(survtime[which(survtime >= failuretimes[k] & left.time < failuretimes[k] &
                                           match.var[,1] >= fail.match[k,1]+match.int[1] & match.var[,1] <= fail.match[k,1]+match.int[2] &
                                           match.var[,2] >= fail.match[k,2]+match.int[3] & match.var[,2] <= fail.match[k,2]+match.int[4])])
      }
      
      H = 1-2*m/(nfail-1)+m*(m-1)/((nfail-1)*(nfail-2))
      H = H/(1-m/(nfail-1))^2
      
      Rho<-matrix(rep(0,renkont^2),ncol=renkont)
      for(i in 1:(renkont-1)) {
        for(j in ((i+1):renkont))  {
          ind = survtime[i] > failuretimes & survtime[j] > failuretimes & left.time[i] < failuretimes & left.time[j] < failuretimes &
            match.var[i,1] >= fail.match[,1]+match.int[1] & match.var[i,1] <= fail.match[,1]+match.int[2] &
            match.var[j,1] >= fail.match[,1]+match.int[1] & match.var[j,1] <= fail.match[,1]+match.int[2] &
            match.var[i,2] >= fail.match[,2]+match.int[3] & match.var[i,2] <= fail.match[,2]+match.int[4] &
            match.var[j,2] >= fail.match[,2]+match.int[3] & match.var[j,2] <= fail.match[,2]+match.int[4]
          Rho[i,j] = prod(H[ind])-1
        }
      }
    }
    
    
    if(dim(match.var)[2] == 3)  {
      
      for(k in 1:length(failuretimes))  {
        nfail[k] = length(survtime[which(survtime >= failuretimes[k] & left.time < failuretimes[k] &
                                           match.var[,1] >= fail.match[k,1]+match.int[1] & match.var[,1] <= fail.match[k,1]+match.int[2] &
                                           match.var[,2] >= fail.match[k,2]+match.int[3] & match.var[,2] <= fail.match[k,2]+match.int[4] &
                                           match.var[,3] >= fail.match[k,3]+match.int[5] & match.var[,3] <= fail.match[k,3]+match.int[6])])
      }
      
      H = 1-2*m/(nfail-1)+m*(m-1)/((nfail-1)*(nfail-2))
      H = H/(1-m/(nfail-1))^2
      
      Rho<-matrix(rep(0,renkont^2),ncol=renkont)
      for(i in 1:(renkont-1)) {
        for(j in ((i+1):renkont))  {
          ind = survtime[i] > failuretimes & survtime[j] > failuretimes & left.time[i] < failuretimes & left.time[j] < failuretimes &
            match.var[i,1] >= fail.match[,1]+match.int[1] & match.var[i,1] <= fail.match[,1]+match.int[2] &
            match.var[j,1] >= fail.match[,1]+match.int[1] & match.var[j,1] <= fail.match[,1]+match.int[2] &
            match.var[i,2] >= fail.match[,2]+match.int[3] & match.var[i,2] <= fail.match[,2]+match.int[4] &
            match.var[j,2] >= fail.match[,2]+match.int[3] & match.var[j,2] <= fail.match[,2]+match.int[4] &
            match.var[i,3] >= fail.match[,3]+match.int[5] & match.var[i,3] <= fail.match[,3]+match.int[6] &
            match.var[j,3] >= fail.match[,3]+match.int[5] & match.var[j,3] <= fail.match[,3]+match.int[6]
          Rho[i,j] = prod(H[ind])-1
        }
      }
    }
    
    if(dim(match.var)[2] > 3)  {
      tt = function(rad)  {
        sum(rad==1)==length(rad)
      }
      
      nfail = rep(0,length(failuretimes))
      ncont4 = function(M,T,V)  {
        mat = matrix(ncol = length(M)+2,nrow=dim(match.var)[1],0)
        mat[which(survtime >= T),1] = 1
        mat[which(left.time < T),2] = 1
        for(i in 1:length(M))  {
          mat[which(match.var[,i] >= M[i]+match.int[(2*i-1)] & match.var[,i] <= M[i]+match.int[(2*i)]),(i+2)] = 1
        }
        sum(apply(mat,1,tt))
      }
      for(k in 1:length(failuretimes))  {
        nfail[k] =  ncont4(fail.match[k,],failuretimes[k],survtime[status!=0][k])
      }
      
      H = 1-2*m/(nfail-1)+m*(m-1)/((nfail-1)*(nfail-2))
      H = H/(1-m/(nfail-1))^2
      
      indeks4 = function(M1,T1,V1,M2,T2,V2)  {
        mat = matrix(ncol = length(M1)*2+4,nrow=dim(fail.match)[1],0)
        mat[which(tidcase <T1),1] = 1
        mat[which(tidcase <T2),2] = 1
        mat[which(leftcase < T1),3] = 1
        mat[which(leftcase < T2),4] = 1
        for(k in 1:length(M1))  {
          mat[which(M1[k] >= fail.match[,k]+match.int[(2*k-1)] & M1[k] <= fail.match[,k]+match.int[(2*k)]),(k+4)] = 1
          mat[which(M2[k] >= fail.match[,k]+match.int[(2*k-1)] & M2[k] <= fail.match[,k]+match.int[(2*k)]),(k+dim(match.var[2])+4)] = 1
        }
        apply(mat,1,tt)
      }
      
      Rho<-matrix(rep(0,renkont^2),ncol=renkont)
      for(i in 1:(renkont-1)) {
        for(j in ((i+1):renkont))  {
          ind = indeks4(kont.match[i,],tidkont[i],left.time[i],kont.match[j,],tidkont[j],left.time[i])
          Rho[i,j] = prod(H[ind])-1
        }
      }
    }
  }
  
  #With left truncation without matching
  if(sum(left.time) > 0 & length(match.var) == 1)  {
    
    nfail = rep(0,length(failuretimes))
    for(i in 1:length(failuretimes))  {
      nfail[i] =  length(survtime[which(survtime >= failuretimes[i] & left.time < failuretimes[i])])
    }
    
    H<-1-2*m/(nfail-1)+m*(m-1)/((nfail-1)*(nfail-2))
    H<-H/(1-m/(nfail-1))^2
    
    rho<-cumprod(H)-1
    indeksT = order(tidkont)
    indeksL = order(leftkont)
    
    Rho<-matrix(rep(0,renkont^2),ncol=renkont)
    for(i in (1:renkont-1)) {
      for(j in ((i+1):renkont))  {
        
        ind = 1:length(failuretimes)
        ind = ind[which(failuretimes >= leftkont[i] & failuretimes >= leftkont[j] & failuretimes<tidkont[i] & failuretimes<tidkont[j])]
        
        Rho[i,j] = prod(H[ind])-1
        
      }
    }
  }
  
  
  
  #Without left truncation and without matching
  if(sum(left.time) == 0 & length(match.var) == 1)  {
    nfail = rep(0,length(failuretimes))
    for(i in 1:length(failuretimes))  {
      nfail[i] =  length(survtime[which(survtime >= failuretimes[i])])
    }
    
    H<-1-2*m/(nfail-1)+m*(m-1)/((nfail-1)*(nfail-2))
    H<-H/(1-m/(nfail-1))^2
    
    rho = rep(0,length(H))
    rho<-cumprod(H)-1
    
    Rho<-matrix(rep(0,renkont^2),ncol=renkont)
    
    for (i in 1:(renkont-1)){
      index<-sum(tidkont[i]>tidcase)
      Rho[i,(i+1):renkont] = rho[index]
    }
  }
  Rhony = Rho + t(Rho)
  p0 = psample[status[samplestat==1]==0]
  Q = matrix((1-p0)/p0^2,ncol=1)%*%matrix((1-p0)/p0^2,nrow=1)
  R<-Rhony*Q+diag((1-p0)/p0^2)
  
  scoreres<-residuals(fit,type="score")
  scoreres =as.matrix(scoreres)
  w<-scoreres[status[samplestat==1]==0,]
  w = as.matrix(w)
  
  Deltavar = t(w)%*%R%*%w
  adjvar<-fit$naive.var+fit$naive.var%*%Deltavar%*%fit$naive.var
}

PoststratVar = function(fit,status,samplestat,survtime,left.time,i,m,no.intervals)  {
  dfb = resid(fit,type="dfbeta")
  if(length(dfb)==sum(samplestat))  {
    dfb = matrix(dfb)
  }
  if(sum(left.time)==0)  {
    partT = seq(min(survtime-0.001),max(survtime),length=(no.intervals+1))
    gamma = 0
    for(i in 1:no.intervals)  {
      index = which(survtime >= partT[i] & survtime <= partT[i+1] & status==0)
      index.ncc = which(survtime[samplestat==1] >= partT[i] & survtime[samplestat==1] <= partT[i+1] & status[samplestat==1]==0)
      m = length(index[samplestat[index]==1])
      n = length(index)
      
      if(m > 1 & n > 0)  {
        gamma = gamma + (1-m/n)*m*var(dfb[index.ncc,])
      }
    }
  }
  else  {
    partT = seq(min(survtime-0.001),max(survtime),length=(no.intervals[2]+1))
    partV = seq(min(left.time-0.001),max(left.time),length=(no.intervals[1]+1))
    gamma = 0
    for(i in 1:no.intervals[2])  {
      for(j in 1:no.intervals[1])  {
        testT = survtime > partT[i]  & survtime <= partT[i+1] & status==0
        testV = left.time > partV[j]  & left.time <= partV[j+1] & status==0
        index = which(testT & testV)
        
        testT.ncc = survtime[samplestat==1] > partT[i]  & survtime[samplestat==1] <= partT[i+1] & status[samplestat==1]==0
        testV.ncc = left.time[samplestat==1] > partV[j]  & left.time[samplestat==1] <= partV[j+1] & status[samplestat==1]==0
        index.ncc = which(testT.ncc & testV.ncc)
        
        m = length(index[samplestat[index]==1])
        n = length(index)
        
        if(m > 1 & n > 0)  {
          gamma = gamma + (1-m/n)*m*var(dfb[index.ncc,])
        }
      }
    }
  }
  adjvar = fit$naive.var + gamma
}


##wrap-functions for probability estimation

KMprob = function(survtime,samplestat,m,left.time=0,match.var=0,match.int=0)  {
  n.cohort = length(survtime)
  status = rep(0,n.cohort)
  status[samplestat > 1] = 1
  o = order(survtime)
  status = status[o]
  survtime = survtime[o]
  if(length(left.time)==n.cohort)  {
    left.time = left.time[o]
  }
  if(length(match.var)==n.cohort)  {
    match.var = match.var[o]
  }
  if(length(match.var) > n.cohort)  {
    match.var = match.var[o,]
  }
  tilbakestill = (1:n.cohort)[o]
  p = pKM(status,survtime,m,n.cohort,left.time,match.var,match.int)
  p[status > 0] = 1
  p = p[order(tilbakestill)]
  p
}

GAMprob = function(survtime,samplestat,left.time=0, match.var=0,match.int=0)  {
  #library(mgcv)
  n.cohort = length(survtime)
  status = rep(0,n.cohort)
  status[samplestat > 1] = 1
  samplestat[samplestat > 1] = 1
  pgam = pGAM(status,survtime,samplestat,n.cohort,left.time,match.var,match.int)
  p = rep(1,n.cohort)
  p[status==0] = pgam
  p
}

GLMprob = function(survtime,samplestat,left.time=0, match.var=0,match.int=0)  {
  n.cohort = length(survtime)
  status = rep(0,n.cohort)
  status[samplestat > 1] = 1
  samplestat[samplestat > 1] = 1
  pglm = pGLM(status,survtime,samplestat,n.cohort,left.time,match.var,match.int)
  p = rep(1,n.cohort)
  p[status==0] = pglm
  p
}


Chenprob = function(survtime,samplestat,no.intervals=10,left.time=0,no.intervals.left = c(3,4))  {
  n.cohort = length(survtime)
  
  status = rep(0,n.cohort)
  status[samplestat > 1] = 1
  samplestat[samplestat > 1] = 1
  ind.no = 1:length(samplestat)
  p = pChen(status,survtime,samplestat,ind.no,n.cohort,no.intervals,left.time,no.intervals.left)
  
  p[status==1] = 1
  p
}

wpl = function(x,data, samplestat, m=1, weight.method="KM", no.intervals=10,
               variance="robust",no.intervals.left = c(3,4), match.var=0, match.int=0)  {
  
  UseMethod("wpl")
}


wpl.default = function(x,data, samplestat, m=1, weight.method="KM",
                       no.intervals=10, variance = "robust", no.intervals.left = c(3,4),
                       match.var=0, match.int=0,...)  {
  
  
  est = wplEst(x,data, samplestat, m, weight.method, no.intervals,
               variance, no.intervals.left, match.var, match.int)
  
  
  class(est) = "wpl"
  est
}

print.wpl = function(x,...)  {
  print.coxph = function (x, digits = max(options()$digits - 4, 3), ...)  {
    if (!is.null(cl <- x$call)) {
      cat("Call:\n")
      dput(cl)
      cat("\n")
    }
    if (!is.null(x$fail)) {
      cat(" Coxph failed.", x$fail, "\n")
      return()
    }
    savedig <- options(digits = digits)
    on.exit(options(savedig))
    coef <- x$coefficients
    se <- sqrt(diag(x$var))
    if (is.null(coef) | is.null(se))
      stop("Input is not valid")
    if (is.null(x$naive.var)) {
      tmp <- cbind(coef, exp(coef), se, coef/se, signif(1 -
                                                          pchisq((coef/se)^2, 1), digits - 1))
      dimnames(tmp) <- list(c(substring(names(coef)[1:length(coef)],2,30000)), c("coef", "exp(coef)",
                                                                                 "se(coef)", "z", "p"))
    }
    else {
      nse <- sqrt(diag(x$naive.var))
      tmp <- cbind(coef, exp(coef), nse, se, coef/se, signif(1 -
                                                               pchisq((coef/se)^2, 1), digits - 1))
      dimnames(tmp) <- list(c(substring(names(coef)[1:length(coef)],2,30000)), c("coef", "exp(coef)",
                                                                                 "se(coef)", "robust se", "z", "p"))
    }
    if(x$est.var==T)  {
      nse <- sqrt(diag(x$naive.var))
      tmp <- cbind(coef, exp(coef), nse, se, coef/se, signif(1 -
                                                               pchisq((coef/se)^2, 1), digits - 1))
      dimnames(tmp) <- list(c(substring(names(coef)[1:length(coef)],2,30000)), c("coef", "exp(coef)",
                                                                                 "se(coef)", "est.se(coef)", "z", "p"))
    }
    cat("\n")
    prmatrix(tmp)
    logtest <- -2 * (x$loglik[1] - x$loglik[2])
    if (is.null(x$df))
      df <- sum(!is.na(coef))
    else df <- round(sum(x$df), 2)
    cat("\n")
    omit <- x$na.action
    cat("  n=", x$n)
    if (!is.null(x$nevent))
      cat(", number of events=", x$nevent, "\n")
    else cat("\n")
    if (length(omit))
      cat("   (", naprint(omit), ")\n", sep = "")
    invisible(x)
  }
  
  #endpoints = length(x)/20
  elements = length(unique(names(x)))
  endpoints = length(x)/elements
  for(i in 1:endpoints)  {
    fit = x[(1+(i-1)*elements):(i*elements)]
    class(fit) = "coxph"
    cat("Endpoint", i,":","\n")
    print.coxph(fit)
    cat("\n")
  }
}

summary.wpl = function(object,...)  {
  summary.coxph = function (object, conf.int = 0.95, scale = 1, ...)  {
    cox <- object
    beta <- cox$coefficients
    names(beta) = c(substring(names(beta)[1:length(beta)],2,30000))
    
    if (is.null(cox$coefficients)) {
      return(object)
    }
    nabeta <- !(is.na(beta))
    beta2 <- beta[nabeta]
    if (is.null(beta) | is.null(cox$var))
      stop("Input is not valid")
    se <- sqrt(diag(cox$var))
    if (!is.null(cox$naive.var))
      nse <- sqrt(diag(cox$naive.var))
    rval <- list(call = cox$call, fail = cox$fail, na.action = cox$na.action,
                 n = cox$n, loglik = cox$loglik)
    if (!is.null(cox$nevent))
      rval$nevent <- cox$nevent
    if (is.null(cox$naive.var)) {
      tmp <- cbind(beta, exp(beta), se, beta/se, 1 - pchisq((beta/se)^2,
                                                            1))
      dimnames(tmp) <- list(names(beta), c("coef", "exp(coef)",
                                           "se(coef)", "z", "Pr(>|z|)"))
    }
    else {
      tmp <- cbind(beta, exp(beta), nse, se, beta/se, 1 - pchisq((beta/se)^2,
                                                                 1))
      dimnames(tmp) <- list(names(beta), c("coef", "exp(coef)",
                                           "se(coef)", "robust se", "z", "Pr(>|z|)"))
    }
    if(cox$est.var==T)  {
      tmp <- cbind(beta, exp(beta), nse, se, beta/se, 1 - pchisq((beta/se)^2,
                                                                 1))
      dimnames(tmp) <- list(names(beta), c("coef", "exp(coef)",
                                           "se(coef)", "est. se(coef)", "z", "Pr(>|z|)"))
    }
    rval$coefficients <- tmp
    if (conf.int) {
      z <- qnorm((1 + conf.int)/2, 0, 1)
      beta <- beta * scale
      se <- se * scale
      tmp <- cbind(exp(beta), exp(-beta), exp(beta - z * se),
                   exp(beta + z * se))
      dimnames(tmp) <- list(names(beta), c("exp(coef)", "exp(-coef)",
                                           paste("lower .", round(100 * conf.int, 2), sep = ""),
                                           paste("upper .", round(100 * conf.int, 2), sep = "")))
      rval$conf.int <- tmp
    }
    if (is.R())
      class(rval) <- "summary.coxph"
    else oldClass(rval) <- "summary.coxph"
    rval
  }
  
  print.summary.coxph = function (x, digits = max(getOption("digits") - 3, 3), signif.stars = getOption("show.signif.stars"), ...)  {
    if (!is.null(x$call)) {
      cat("Call:\n")
      dput(x$call)
      cat("\n")
    }
    if (!is.null(x$fail)) {
      cat(" Coxreg failed.", x$fail, "\n")
      return()
    }
    savedig <- options(digits = digits)
    on.exit(options(savedig))
    omit <- x$na.action
    cat("  n=", x$n)
    if (!is.null(x$nevent))
      cat(", number of events=", x$nevent, "\n")
    else cat("\n")
    if (length(omit))
      cat("   (", naprint(omit), ")\n", sep = "")
    if (nrow(x$coef) == 0) {
      cat("   Null model\n")
      return()
    }
    if (!is.null(x$coefficients)) {
      cat("\n")
      if (is.R())
        printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars,
                     ...)
      else prmatrix(x$coefficients)
    }
    if (!is.null(x$conf.int)) {
      cat("\n")
      prmatrix(x$conf.int)
    }
    cat("\n")
    invisible()
  }
  
  
  #endpoints = length(object)/20
  elements = length(unique(names(object)))
  endpoints = length(object)/elements
  for(i in 1:endpoints)  {
    fit = object[(1+(i-1)*elements):(i*elements)]
    class(fit) = "coxph"
    cat("Endpoint", i,":","\n")
    print.summary.coxph(summary.coxph(fit))
    cat("\n")
  }
}


wpl.formula = function(formula,data=data, samplestat,...)  {
  
  mf = model.frame(formula=formula,data=data)
  
  x.coef = model.matrix(attr(mf,"terms"),data=mf)
  x.coef = as.matrix(x.coef[,-1])
  ant.coef = dim(x.coef)[2]
  y = model.response(mf)
  
  if(dim(y)[2] == 2)  {
    left.time = 0
    survtime = y[,1]
    status = y[,2]
  }
  
  if(dim(y)[2] == 3)  {
    left.time = y[,1]
    survtime = y[,2]
    status = y[,3]
  }
  
  
  x = array(data=list(x.coef,y))
  
  est = wpl.default(x, data, samplestat,...)
  est$call = match.call()
  est$formula = formula
  
  est
}
################################################################################
################################################################################