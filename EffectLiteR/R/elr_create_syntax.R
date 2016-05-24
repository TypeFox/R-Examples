


create_syntax_intercepts <- function(y, alphas){

  res <- paste0(y," ~ c(", paste(alphas[1,,],collapse=","), ")*1")
  return(res)
  
}


create_syntax_regcoef <- function(y, z, nz, alphas){
  
  res <- NULL
  if (nz>0) {
    for (i in 1:nz) {
      tmp <- paste0(y," ~ c(", paste(alphas[i+1,,],collapse=","), ")*",z[i])
      res <- paste0(res, "\n", tmp)
    }
  }
  return(res)
  
}


create_syntax_cellmeanz <- function(z, nz, fixed.z, cellmeanz, sampmeanz){
  
  res <- NULL
  if(!fixed.z){
    ## stochastic z
    if (nz>0) {
      cellmeanz <- matrix(cellmeanz, nrow=nz)    
      for (i in 1:nz) {
        tmp <- paste0(z[i]," ~ c(", paste(cellmeanz[i,],collapse=","), ")*1")
        res <- paste0(res, "\n", tmp)
      }
    }
    
  }else if(fixed.z){
    ## fixed z
    if (nz>0) {
      res <- paste0(res, "\n\n## Fixed Means of Z")
      cellmeanz <- matrix(cellmeanz, nrow=nz)
      tmp <- paste0(cellmeanz, " := ",  sampmeanz, collapse="\n")
      res <- paste0(res, "\n", tmp)
    }
  }

  return(res)
}


create_syntax_covz <- function(z, nz, fixed.z){

  ## covariances between stochastic z
  res <- NULL
  
  ## syntax covariances between z in each cell
  if(!fixed.z){
    if(nz > 1){
      tmp <- combn(z,2)
      res <- paste0("\n", paste0(tmp[1,], " ~~ ", tmp[2,], collapse="\n"))
    }
  }
  
  return(res)
}
  

create_syntax_homoscedasticity <- function(y, ng, nk, homoscedasticity){

  res <- NULL
  if(homoscedasticity){
    tmp <- paste0(y, " ~~ c(", 
                  paste(rep("veps", times=ng*nk),collapse=","),
                  ")*", y)
    res <- paste0(res, "\n", tmp)
  }
  
  return(res)
}


create_syntax_group_freq <- function(fixed.cell, relfreq, observed.freq, groupw){
  
  
  res <- "\n\n## Relative Group Frequencies \n"
  
  if(fixed.cell){
    tmp <- paste(paste0(relfreq, " := ", observed.freq), collapse="\n")
    res <- paste0(res, tmp)       
    
  }else if(!fixed.cell){
    
    ## syntax group weights
    tmp <- paste0("group % c(", paste(groupw, collapse=","), ")*w")
    res <- paste0(res, tmp)
    
    tmp <- paste(paste0("exp(", groupw, ")"), collapse=" + ")
    tmp <- paste0("N := ",tmp)
    res <- paste0(res, "\n", tmp)
    
    tmp <- paste(paste0(relfreq, " := exp(", groupw, ")/N"), collapse="\n")
    res <- paste0(res, "\n", tmp)    
  }
  
  return(res)
}



create_syntax_betas <- function(betas, alphas, ng, nk, nz){
  
  res <- "\n\n## beta Coefficients"
  
  ## create temporary beta array (to be overwritten in next step)
  beta <- betas
  
  ## compute betas based on alphas
  for(q in 1:(nz+1)){
    for(x in 1:ng){
      beta[q,1,x] <- paste(betas[q,1,x], alphas[q,1,x], sep=" := ")
      if(nk>1){
        for(k in 2:nk){
          beta[q,k,x] <- paste0(betas[q,k,x], " := ", 
                                alphas[q,k,x], "-", alphas[q,1,x])
        }        
      }
    }
  }
  res <- paste0(res, "\n", paste(beta, collapse="\n"))
  
  return(res)
}



create_syntax_gammas <- function(gammas, betas, ng, nk, nz){

  res <- "\n\n## gamma Coefficients"
  
  ## create temporary gamma array (to be overwritten in next step)
  gamma <- gammas
  
  ## compute gammas based on betas
  for(q in 1:(nz+1)){
    for(k in 1:nk){
      gamma[q,k,1] <- paste(gammas[q,k,1], betas[q,k,1], sep=" := ")
      for(x in 2:ng){
        gamma[q,k,x] <- paste0(gammas[q,k,x], " := ", 
                               betas[q,k,x], "-", betas[q,k,1])
      }
    }
  }
  res <- paste0(res, "\n", paste(gamma, collapse="\n"))
  
  return(res)
}


create_syntax_ez <- function(nz, meanz, cellmeanz, relfreq){
  
  res <- NULL
  if (nz>0) {
    res <- paste0(res, "\n\n## Unconditional Expectations E(Z)")
    cellmeanz <- matrix(cellmeanz, nrow=nz)    
    for (i in 1:nz) {
      tmp <- paste0(meanz[i]," := ", paste(cellmeanz[i,], relfreq, sep="*", collapse=" + "))
      res <- paste0(res, "\n", tmp)
    }
  }
  
  return(res)
}  


create_syntax_pk <- function(nk, pk, relfreq){
  
  res <- NULL
  if (nk>1) {
    res <- paste0(res, "\n\n## Unconditional Probabilities P(K=k)")
    relfreq <- matrix(relfreq, nrow=nk)
    for (i in 1:nk) {
      tmp <- paste0(pk[i], " := ", paste(relfreq[i,], collapse=" + "))
      res <- paste0(res, "\n", tmp)
    }
  }
  
  return(res)
}


create_syntax_ex <- function(px, ng, nk, relfreq){
  
  res <- "\n\n## Unconditional Probabilities P(X=x)"
  
  relfreq <- matrix(relfreq, nrow=nk)
  for (i in 1:ng) {
    tmp <- paste0(px[i], " := ", paste(relfreq[,i], collapse=" + "))
    res <- paste0(res, "\n", tmp)
  }
  
  return(res)
}


create_syntax_Ezk <- function(ng, nk, nz, Ezk, cellmeanz, relfreq){
  
  res <- NULL
  if (nk>1 & nz>0) {
    
    res <- paste0(res, "\n\n## Unconditional Expectations E(Z*I_K=k)")
    Ezk <- array(Ezk, dim=c(nk,nz))
    cellmeanz <- array(cellmeanz, dim=c(nz,nk,ng))
    relfreq <- array(relfreq, dim=c(nk,ng))
    
    for(q in 1:nz){
      for(k in 1:nk){
        tmp <- paste0(Ezk[k,q], " := ", 
                      paste(cellmeanz[q,k,1:ng], relfreq[k,], sep="*", collapse=" + "))
        res <- paste0(res, "\n", tmp)
      }
    }
    
  }
  
  return(res)
}


create_syntax_Egx <- function(ng,nk,nz,pk,meanz,Ezk,Egx,gammas){
  
  res <- "\n\n## Average Effects"
  
  ## create vector of unconditional expectations of Z, K, Z*K
  pk_tmp <- pk[-1]
  meanz_tmp <- Ezk_tmp <- pkEzk_tmp <- character()
  
  if(nz>0){meanz_tmp <- meanz}  
  if(nk>1 & nz>0){Ezk_tmp <- c(matrix(Ezk, nrow=nk)[-1,])}  
  if(nk>1){pkEzk_tmp <- c(matrix(c(pk_tmp,Ezk_tmp), ncol=nk-1, byrow=T))}
  
  expectations <- c(1,meanz_tmp,pkEzk_tmp)
  
  ## average total effects  
  for(i in 2:ng){
    tmp <- paste0(Egx[i-1]," := ",
                  paste(gammas[,,i],expectations, sep="*", collapse=" + "))
    res <- paste0(res, "\n", tmp)
  }
  
  return(res)
}


create_syntax_adjmeans <- function(ng,nk,nz,pk,meanz,Ezk,adjmeans,betas){
  
  ##TODO compute adjusted means based on gammas so that it can be used in
  ## single group models as well
  
  res <- "\n\n## Adjusted Means"
  
  ## create vector of unconditional expectations of Z, K, Z*K
  pk_tmp <- pk[-1]
  meanz_tmp <- Ezk_tmp <- pkEzk_tmp <- character()
  
  if(nz>0){meanz_tmp <- meanz}  
  if(nk>1 & nz>0){Ezk_tmp <- c(matrix(Ezk, nrow=nk)[-1,])}  
  if(nk>1){pkEzk_tmp <- c(matrix(c(pk_tmp,Ezk_tmp), ncol=nk-1, byrow=T))}
  
  expectations <- c(1,meanz_tmp,pkEzk_tmp)
  
  ## adjusted means
  for(i in 1:ng){
    tmp <- paste0(adjmeans[i]," := ",
                  paste(betas[,,i],expectations, sep="*", collapse=" + "))
    res <- paste0(res, "\n", tmp)
  }
  
  
  return(res)
}


create_syntax_Pkgx <- function(ng, nk, relfreq, Pkgx, px){
  
  res <- "\n\n## Conditional Probabilities P(K=k|X=x)"
  
  relfreq <- matrix(relfreq, nrow=nk)
  Pkgx <- matrix(Pkgx, nrow=nk)
  
  for(i in 1:ng){
    for(k in 1:nk){
      Pkgx[k,i] <- paste0(Pkgx[k,i], " := ", relfreq[k,i], "/", px[i])
    }
  }   
  res <- paste0(res, "\n", paste(Pkgx, collapse="\n"))  
  
  return(res)
}


create_syntax_Pxgk <- function(ng, nk, relfreq, Pxgk, pk){
  
  res <- NULL
  if(nk>1){
    res <- paste0(res, "\n\n## Conditional Probabilities P(X=x|K=k)")
    relfreq <- matrix(relfreq, nrow=nk)
    Pxgk <- matrix(Pxgk, nrow=ng)

    for(i in 1:ng){
      for(k in 1:nk){
        Pxgk[i,k] <- paste0(Pxgk[i,k], " := ", relfreq[k,i], "/", pk[k])
      }
    }   
    res <- paste0(res, "\n", paste(Pxgk, collapse="\n"))      
  }
  
  return(res)
}


create_syntax_Ezgx <- function(ng,nk,nz,Ezgx,Pkgx,cellmeanz){
  
  res <- NULL
  
  if(nz!=0){
    res <- paste0(res, "\n\n## Conditional Expectations E(Z|X=x)")
    cellmeanz <- array(cellmeanz, dim=c(nz,nk,ng))
    Ezgx <- matrix(Ezgx, nrow=nz)
    Pkgx <- matrix(Pkgx, nrow=nk)
    for(i in 1:ng){
      for(q in 1:nz){
        Ezgx[q,i] <- paste0(Ezgx[q,i], " := ", 
                            paste(cellmeanz[q,1:nk,i], Pkgx[,i], sep="*", collapse=" + "))
      }
    }
    res <- paste0(res, "\n", paste(Ezgx, collapse="\n"))
  }
  
  return(res)
}


create_syntax_Ezgk <- function(ng,nk,nz,cellmeanz,Ezgk,Pxgk){
  
  res <- NULL

  if(nz>0 & nk>1){
    res <- paste0(res, "\n\n## Conditional Expectations E(Z|K=k)")
    cellmeanz <- array(cellmeanz, dim=c(nz,nk,ng))
    Ezgk <- matrix(Ezgk, nrow=nz)
    Pxgk <- matrix(Pxgk, nrow=ng)
    
    for(i in 1:nk){
      for(q in 1:nz){
        Ezgk[q,i] <- paste0(Ezgk[q,i], " := ", 
                            paste(cellmeanz[q,i,1:ng], Pxgk[,i], sep="*", collapse=" + "))
      }
    }
    res <- paste0(res, "\n", paste(Ezgk, collapse="\n"))
  }  
  
  return(res)
}


create_syntax_Ezkgx <- function(ng,nk,nz,Ezkgx,cellmeanz,Pkgx){
  
  res <- NULL

  if(nz>0 & nk>1){
    res <- paste0(res, "\n\n## Conditional Expectations E(Z*I_K=k|X=x)")
    Ezkgx <- array(Ezkgx, dim=c(nz,nk,ng))
    cellmeanz <- array(cellmeanz, dim=c(nz,nk,ng))
    Pkgx <- array(Pkgx, dim=c(nk,ng))
    
    for(i in 1:ng){
      for(k in 1:nk){
        for(q in 1:nz){
          Ezkgx[q,k,i] <- paste0(Ezkgx[q,k,i], " := ", cellmeanz[[q,k,i]],
                                 "*", Pkgx[k,i])
        }
      }
    }
    res <- paste0(res, "\n", paste(Ezkgx, collapse="\n"))  
  }
  
  return(res)
}


create_syntax_Egxgx <- function(ng,nk,nz,Pkgx,Ezgx,Ezkgx,Egxgx,gammas){
  
  ## create matrix of conditional expectations of Z, K, Z*K given X=x
  ##TODO: maybe we can get rid of the if conditions and only use last one
  expectationsgx <- character()
  if(nz==0 & nk==1){
    expectationsgx <- matrix("1", nrow=ng)
  }
  if(nz>0 & nk==1){
    Ezgx <- matrix(Ezgx, nrow=ng, byrow=T)
    expectationsgx <- cbind("1", Ezgx)
  }
  if(nz==0 & nk>1){    
    Pkgx <- matrix(Pkgx, nrow=ng, byrow=T)[,-1]
    expectationsgx <- cbind("1", Pkgx)
  }
  if(nz>0 & nk>1){
    Ezgx <- matrix(Ezgx, nrow=ng, byrow=T)
    
    select <- c(matrix(Pkgx, nrow=ng, byrow=T)[,-1])
    Pkgx <- matrix(select, nrow=ng)
    
    select <- c(array(Ezkgx, dim=c(nz,nk,ng))[,-1,])
    Ezkgx <- array(select, dim=c(nz,nk-1,ng))
    
    tmp1 <- matrix(Pkgx[1,], nrow=nk-1, ncol=1)
    tmp2 <- matrix(Ezkgx[,,1], nrow=nk-1, ncol=nz, byrow=TRUE)
    expectationsgx <- c(1, Ezgx[1,], c(t(cbind(tmp1, tmp2))))
    
    for(i in 2:ng){
      tmp1 <- matrix(Pkgx[i,], nrow=nk-1, ncol=1)
      tmp2 <- matrix(Ezkgx[,,i], nrow=nk-1, ncol=nz, byrow=TRUE)      
      expectationsgx <- rbind(expectationsgx,
                              c(1, Ezgx[i,], c(t(cbind(tmp1, tmp2)))))
    }
    
  }
  
  ## Effects given a treatment condition
  res <-  "\n\n## Effects given X=x"
  Egxgx <- matrix(Egxgx, ncol=ng)
  
  for(gx in 1:ng){
    for(i in 2:ng){
      tmp <- paste0(Egxgx[i-1,gx], " := ", 
                    paste(c(gammas[,,i]),expectationsgx[gx,], sep="*", collapse=" + "))
      res <- paste0(res, "\n", tmp)
      
    }    
  }
  
  return(res)
}


create_syntax_Egxgk <- function(ng,nk,Egxgk,Ezgk,gammas){
  
  res <- NULL
  Egxgk <- matrix(Egxgk, nrow=ng-1)
  Ezgk <- matrix(Ezgk, ncol=nk)
  Ezgk <- rbind("1", Ezgk)
  
  if(nk>1){
    res <- paste0(res, "\n\n## Effects given K=k")
    for(i in 2:ng){     ## effects given K=0
      tmp <- paste0(Egxgk[i-1,1], " := ", 
                    paste(gammas[,1,i], Ezgk[,1], sep="*", collapse=" + ")) 
      res <- paste0(res, "\n", tmp)
    }
    
    for(i in 2:ng){     ## effects given K=1,...,k
      for(k in 2:nk){
        gammaselect <- c(gammas[,c(1,k),i])
        tmp <- paste0(Egxgk[i-1,k], " := ", 
                      paste(gammaselect, Ezgk[,k], sep="*", collapse=" + "))
        res <- paste0(res, "\n", tmp)
      }
    }
    
  }
  
  return(res)
}


create_syntax_Egxgxk <- function(ng,nk,nz,Egxgxk,gammas,cellmeanz){
  
  res <- NULL
  
  if(nk>1 & nz>0){
    res <- paste0(res, "\n\n## Effects given X=x, K=k")
    cellmeanz <- array(cellmeanz, dim=c(nz,nk,ng))
    Egxgxk <- array(Egxgxk, dim=c(ng-1,ng,nk))
    
    
    for(t in 2:ng){
      for(x in 1:ng){
        for(k in 1:nk){        
          if(k==1){
            rhs <- paste0(gammas[,1,t],"*", c(1,cellmeanz[,1,x]), collapse=" + ")
            tmp <- paste(Egxgxk[t-1,x,k], ":=", rhs)
            res <- paste0(res, "\n", tmp) 
          }else{
            gammaselect <- c(gammas[,1,t], gammas[,k,t])
            cellmeanzselect <- rep(c(1,cellmeanz[,k,x]),2)
            rhs <- paste0(gammaselect,"*", cellmeanzselect, collapse=" + ")
            tmp <- paste(Egxgxk[t-1,x,k], ":=", rhs)
            res <- paste0(res, "\n", tmp)
          }
        }
      }
    }
    
  }
  
  return(res)
}


create_syntax_interaction_constraints <- function(ng,nk,nz,interactions,gammas){
  
  res <- NULL
  
  stopifnot(interactions %in% c("all","none","2-way","X:K","X:Z"))
  if(interactions == "none"){
    res <- paste0(res, "\n\n## Equality Constraints")
    gammas <- matrix(c(gammas), ncol=ng)[-1,-1]
    res <- paste0(res, "\n", paste(gammas, "== 0", collapse="\n"))
  } 
  if(interactions == "2-way"){
    if(nk>1 & nz>0){
      res <- paste0(res, "\n\n## Equality Constraints")
      gammas <- gammas[2:(nz+1), 2:nk, 2:ng]
      res <- paste0(res, "\n", paste(gammas, "== 0", collapse="\n"))
    }
  }
  if(interactions == "X:K"){
    if(nz>0){
      res <- paste0(res, "\n\n## Equality Constraints")
      gammas <- gammas[2:(nz+1), , 2:ng]
      res <- paste0(res, "\n", paste(gammas, "== 0", collapse="\n"))
    }
  }
  if(interactions == "X:Z"){
    if(nk>1){
      res <- paste0(res, "\n\n## Equality Constraints")
      gammas <- gammas[, 2:nk, 2:ng]
      res <- paste0(res, "\n", paste(gammas, "== 0", collapse="\n"))
    }
  }

  return(res)  
}



create_syntax_hypotheses <- function(ng, Egx, gammas){
  
  ## Hypothesis 1: No average treatment effects
  hypothesis1 <- paste(Egx, "== 0", collapse="\n")
  
  ## Hypothesis 2: No covariate effects in control group
  gammas_tmp <- matrix(c(gammas), ncol=ng)[-1,1]
  hypothesis2 <- paste(gammas_tmp, "== 0", collapse="\n")
  
  ## Hypothesis 3: No treatment*covariate interaction
  gammas_tmp <- matrix(c(gammas), ncol=ng)[-1,-1]
  hypothesis3 <- paste(gammas_tmp, "== 0", collapse="\n")
  
  ## Hypothesis 4: No treatment effects
  gammas_tmp <- matrix(c(gammas), ncol=ng)[,-1]
  hypothesis4 <- paste(gammas_tmp, "== 0", collapse="\n")
  
  hypotheses <- list(hypothesis1=hypothesis1,
                     hypothesis2=hypothesis2,
                     hypothesis3=hypothesis3,
                     hypothesis4=hypothesis4)
  
  return(hypotheses)
}





createLavaanSyntax <- function(obj) {
  
  inp <- obj@input
  parnames <- obj@parnames
  
  ## input information
  y <- inp@vnames$y
  z <- inp@vnames$z  
  ng <- inp@ng
  nz <- inp@nz
  nk <- inp@nk
  fixed.cell <- inp@fixed.cell
  fixed.z <- inp@fixed.z
  sampmeanz <- inp@sampmeanz
  homoscedasticity <- inp@homoscedasticity
  observed.freq <- inp@observed.freq
  interactions <- inp@interactions
  
  ## parnames
  alphas <- parnames@alphas
  betas <- parnames@betas
  gammas <- parnames@gammas
  cellmeanz <- parnames@cellmeanz
  relfreq <- parnames@relfreq
  groupw <- parnames@groupw
  meanz <- parnames@meanz
  pk <- parnames@pk
  px <- parnames@px
  Ezk <- parnames@Ezk
  Egx <- parnames@Egx
  adjmeans <- parnames@adjmeans
  Pkgx <- parnames@Pkgx
  Pxgk <- parnames@Pxgk
  Ezgx <- parnames@Ezgx
  Ezgk <- parnames@Ezgk
  Ezkgx <- parnames@Ezkgx
  Egxgx <- parnames@Egxgx
  Egxgk <- parnames@Egxgk
  Egxgxk <- parnames@Egxgxk
  
  
  model <- "#### lavaan Syntax for EffectLiteR Model ####"
  
  ## measurement model
  if(length(inp@measurement) != 0){
    model <- paste0(model, "\n\n## Measurement Model \n")
    model <- paste0(model, inp@measurement)
  }
  
  ## syntax intercepts
  model <- paste0(model, "\n\n## Structural Model \n")
  model <- paste0(model, create_syntax_intercepts(y,alphas))

  ## syntax regression coefficients in each cell
  model <- paste0(model, create_syntax_regcoef(y,z,nz,alphas))
  
  ## mean z in each cell
  model <- paste0(model, create_syntax_cellmeanz(z, nz, fixed.z, cellmeanz, 
                                                 sampmeanz))
  
  ## covariances between stochastic z
  model <- paste0(model, create_syntax_covz(z, nz, fixed.z))

  ## homoscedastic residual variances
  model <- paste0(model, create_syntax_homoscedasticity(y,ng,nk,homoscedasticity))
  
  ## compute relative group frequencies
  model <- paste0(model, create_syntax_group_freq(fixed.cell, relfreq, 
                                                  observed.freq, groupw))
  
  ## compute betas based on alphas
  model <- paste0(model, create_syntax_betas(betas, alphas, ng, nk, nz))
  
  ## compute gammas based on betas
  model <- paste0(model, create_syntax_gammas(gammas, betas, ng, nk, nz))
  
  ## compute unconditional means of z
  model <- paste0(model, create_syntax_ez(nz, meanz, cellmeanz, relfreq))
  
  ## compute unconditional probabilities of K*=k
  model <- paste0(model, create_syntax_pk(nk, pk, relfreq))
  
  ## compute unconditional probabilities of X=x
  model <- paste0(model, create_syntax_ex(px, ng, nk, relfreq))
  
  ## compute unconditional means of Z*K 
  model <- paste0(model, create_syntax_Ezk(ng, nk, nz, Ezk, cellmeanz, relfreq))

  ## compute average effects
  model <- paste0(model, create_syntax_Egx(ng,nk,nz,pk,meanz,Ezk,Egx,gammas))  

  ## compute adjusted means
  model <- paste0(model, create_syntax_adjmeans(ng,nk,nz,pk,meanz,Ezk,adjmeans,betas))
  
  ## conditional probabilities of K=k given X=x (Pkgx)
  model <- paste0(model, create_syntax_Pkgx(ng, nk, relfreq, Pkgx, px))
  
  ## conditional probabilities of X=x given K=k (Pxgk)
  model <- paste0(model, create_syntax_Pxgk(ng, nk, relfreq, Pxgk, pk))
  
  ## conditional expectations of Z given X=x (Ezgx)
  model <- paste0(model, create_syntax_Ezgx(ng, nk, nz, Ezgx, Pkgx, cellmeanz))
  
  ## conditional expectations of Z given K=k (Ezgk)
  model <- paste0(model, create_syntax_Ezgk(ng,nk,nz,cellmeanz,Ezgk,Pxgk))
  
  ## conditional expectations of Z*K given X=x (Ezkgx)
  model <- paste0(model, create_syntax_Ezkgx(ng,nk,nz,Ezkgx,cellmeanz,Pkgx))
  
  ## effects given a treatment condition
  model <- paste0(model, create_syntax_Egxgx(ng,nk,nz,Pkgx,Ezgx,Ezkgx,Egxgx,gammas))
  
  ## Effects given a value k of K  
  model <- paste0(model, create_syntax_Egxgk(ng,nk,Egxgk,Ezgk,gammas))
  
  ## Effects given X=x and K=k
  model <- paste0(model, create_syntax_Egxgxk(ng,nk,nz,Egxgxk,gammas,cellmeanz))
  
  ## Constraints about 2 and 3 way interactions
  model <- paste0(model, 
              create_syntax_interaction_constraints(ng,nk,nz,interactions,gammas))
  
  ## additional syntax
  if(length(inp@add) != 0){
    model <- paste0(model, "\n\n## Additional User Defined Syntax \n")
    model <- paste0(model, inp@add)
  }
  
  ## syntax for main hypotheses
  hypotheses <- create_syntax_hypotheses(ng, Egx, gammas)
    
  
  res <- new("lavsyntax",             
             model=model,
             hypotheses=hypotheses
  )
  
  return(res)
  
}
