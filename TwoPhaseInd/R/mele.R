mele <- function(data,			## dataset (data.frame)
                 response,		## response variable (binary)
                 treatment,		## treatment variable (binary)       ### 1st phase
                 BaselineMarker,		## environment variable (continuous) ### 2nd phase
                 extra=NULL,		## extra variable(s)
                 phase,			## variable for phase indicator
                 ind=TRUE,		## independent or non-indepentent
                 maxit=2000  	## max iteration
                 ){

  ## check if the data argument is missing
  call <- match.call()
  if (missing(data)){ 
    data <- environment(formula)
  }

  ## argument validation
  if(!is.data.frame(data)){
    stop("Argument data must be a data.frame object.")
  }else{
    colNames <- colnames(data)
    if(!(response %in% colNames)){
      stop("Response variable not found in the data.")
    }else{
      if(any(levels(factor(data[, response])) != c("0", "1"))){
        warning("Response variable must be either 0 or 1 only.")
      }
    }
    if(!(treatment %in% colNames)){
      stop("Treatment variable not found in the data.")
    }else{
      if(any(levels(factor(data[, treatment])) != c("0", "1"))){
        warning("Treatment variable must be either 0 or 1 only.")
      }
    }
    if(!(BaselineMarker %in% colNames)){
      stop("BaselineMarker variable not found in the data.")
    }
    if(!is.null(extra)){
      if(!any(extra %in% colNames)){
        extraNotFound <- paste(extra[!(extra %in% colNames)], sep="", collapse=", ")
        stop(paste("Extra variable(s) not found in the data:", extraNotFound))
      }
      tmp <- remove_rarevariants(data[,extra])
      if (any(tmp))
      {
        idx <- tmp==TRUE
        toremove=NULL
        for (i in 1:length(idx))
        {
          if (idx[i]) toremove <- c(toremove,extra[idx[i]])
        }
        warnings(paste0(paste(toremove,sep=", "), " were removed due to rare vairant"))
        extra <- extra[!idx]
        if (length(extra)==0) extra <- NULL
      }
    }

    if(!(phase %in% colNames)){
      stop("Phase variable not found in the data.")
    }else{
      if(any(levels(factor(data[, phase])) != c("1", "2"))){
        stop("Phase variable must be either 1 or 2 only.")
      }
    }
  }
  
  ##check response vaiable
  idx <- data[,response]==0 | data[,response]==1
  if (any(idx==FALSE))
  {
    data <- data[idx,]
    warning(paste0(sum(!idx)," rows were removed, whose response !=0 or 1"))
  }
  ##check treatment vaiable
  idx <- data[,treatment]==0 | data[,treatment]==1
  if (any(idx==FALSE))
  {
    data <- data[idx,]
    warning(paste0(sum(!idx)," rows were removed, whose treatment !=0 or 1"))
  }
  
  ##check phase variable
  idx <- data[,phase]==1 | data[,phase]==2
  if (any(idx==FALSE))
  {
    data <- data[idx,]
    warning(paste0(sum(!idx)," rows were removed, whose phase != 1 or 2"))
  }
  
  ##check if extra variables in phase2 are missing
  if (!is.null(extra))
  {
    nextra_remove <- 0
    vars=NULL
    for (var in extra)
    {
      idx <- data[,phase]==2 & (is.na(data[,var]) | is.null(data[,var]))
      if (any(idx==TRUE))
      {
        nextra_remove <- nextra_remove + sum(idx)
        data[idx,phase]=3
        vars=c(vars,var)
      }
    }
    idx <- data[,phase]==3
    if (any(idx==TRUE))
    {
      data <- data[!idx,]
      missingcols=NULL
      for (i in 1:length(vars))
      {
        missingcols=paste0(missingcols," ",vars[i])
      }
      warning(paste0(nextra_remove, "rows in phase2 were removed due to missing data in",missingcols))
    }
  }
  
  ##check BaselineMarker
  idx <- data[,phase]==2 & (is.na(data[,BaselineMarker]) | is.null(data[,BaselineMarker]))
  if (any(idx==TRUE))
  {
    data <- data[!idx,]
    warning(paste0(sum(idx)," rows were removed, whose BaselineMarker is missing in phase2"))
  }
  #indicator of divergence
  converged=1
  if (remove_rarevariants(data[, BaselineMarker]))
  {
    warnings("BaselineMarker variable is rare variant")
    tmpResult <- data.frame(beta <- rep(NA,length(extra)+4), stder=rep(NA,length(extra)+4), pVal=rep(NA,length(extra)+4))
    rownames(tmpResult)[1:4] <- c("(Intercept)","x1","z1","x1:z1")
    rownames(tmpResult)[5:nrow(tmpResult)] <- extra
  }else
  {
    isPhase1 <- data[, phase] == 1
    isPhase2 <- data[, phase] == 2
    
    a <- data[, response]
    b <- data[, treatment]
    
    N1 <- sum(a==1)
    N0 <- sum(a==0)
    N  <- N1 + N0
    NXY1 <- sum(b==0 & a==0)
    NXY2 <- sum(b==1 & a==0)
    NXY3 <- sum(b==0 & a==1)
    NXY4 <- sum(b==1 & a==1)
    
    if (NXY1==0 | NXY2==0 | NXY3==0 |NXY4==0)
    {
      warning(paste0("(x==0 & y==0):",NXY1,
                     " (x==1 & y==0):",NXY2,
                     " (x==0 & y==1):",NXY3,
                     " (x==1 & y==1):",NXY4))
      stop("A treatment-response category has zero counts!")
    }
    y <- data[isPhase2, response]
    x <- data[isPhase2, treatment]
    z <- data[isPhase2, BaselineMarker]
    
    n0 <- sum(y==0)
    n1 <- sum(y==1)
    n  <- n0 + n1    
    
    if(!ind){
      ## using independence (start)
      n11 <- sum(y==1 & x==1)
      n10 <- sum(y==1 & x==0)
      n01 <- sum(y==0 & x==1)
      n00 <- sum(y==0 & x==0)  
      z0  <- z[x==0]
      z1  <- z[x==1]      
      nz0 <- length(z0)
      nz1 <- length(z1)
      z1  <- c(z1,z0,z1,z0)
      from.y <- c(y[x==1],y[x==0],y[x==1],y[x==0])
      from.x <- c(x[x==1],x[x==0],x[x==1],x[x==0])
      x1 <- c(rep(1,nz1),rep(0,nz0),rep(1,nz1),rep(0,nz0))
      y1 <- c(rep(1,n),rep(0,n))
      
      pzx <- rep(0,2*n)
      Nyx <- pzx
      pyx <- pzx
      Nyx[x1==0 & y1==0] <- NXY1 - sum(x==0 & y==0)
      pyx[x1==0 & y1==0] <- NXY1/(NXY1+NXY3)  
      Nyx[x1==1 & y1==0] <- NXY2 - sum(x==1 & y==0)
      pyx[x1==1 & y1==0] <- NXY2/(NXY2+NXY4)
      Nyx[x1==0 & y1==1] <- NXY3 - sum(x==0 & y==1)
      pyx[x1==0 & y1==1] <- NXY3/(NXY1+NXY3)
      Nyx[x1==1 & y1==1] <- NXY4 - sum(x==1 & y==1)
      pyx[x1==1 & y1==1] <- NXY4/(NXY2 + NXY4)
      weight <- rep(0,length(y1))
      for (i in 1:length(y1))  weight[i] <- 1/sum(x==from.x[i] & y==from.y[i])    
      pzx[from.y==1 & from.x==1] <- weight[from.y==1 & from.x==1] * unique(pyx[y1==1 & x1==1])
      pzx[from.y==1 & from.x==0] <- weight[from.y==1 & from.x==0] * unique(pyx[y1==1 & x1==0])
      pzx[from.y==0 & from.x==1] <- weight[from.y==0 & from.x==1] * unique(pyx[y1==0 & x1==1])
      pzx[from.y==0 & from.x==0] <- weight[from.y==0 & from.x==0] * unique(pyx[y1==0 & x1==0])
      ## using independence (stop)
    }else{
      ## without using independence (start)
      z.from.y <- ifelse(y==1,1,0)
      x1 <- c(rep(1,n),rep(0,n),rep(1,n),rep(0,n))
      y1 <- c(rep(1,2*n),rep(0,2*n))
      z1 <- rep(z,4)
      z.from.y<- rep(z.from.y,4)
      
      pzw <- ifelse(z.from.y==1,1/n1 *(N1/N), 1/n0 *(N0/N) )
      Nyx  <- rep(0,length(y1))
      for (i in 1:length(y1)) {
        if (x1[i]==0 & y1[i]==0) {
          Nyx[i] <- NXY1 - sum(x==0 & y==0)                     
        }
        if (x1[i]==1 & y1[i]==0) { 
          Nyx[i] <- NXY2 - sum(x==1 & y==0)                     
          
        } 
        if (x1[i]==0 & y1[i]==1) {
          Nyx[i] <- NXY3 - sum(x==0 & y==1)                     
          
        }
        if (x1[i]==1 & y1[i]==1) {
          Nyx[i] <- NXY4 - sum(x==1 & y==1)                     
        } 
      }
      ## without using independence (end)
    }
    
    X1 <- model.matrix(~ x1 + z1 + x1*z1) ## 1st phase
    X  <- model.matrix(~ x  + z  + x*z)   ## 2nd phase
    if(is.null(extra)){
      ## without adjusting the covariates 
      glmData <- data.frame(cbind(y, x, z, x*z))
    }else{
      ## adjusting the covariates 
      covar <- data[isPhase2, extra]
      for(var in extra){
        if(is.factor(covar[, var])){
          covar[, var] <- factor(covar[, var])
        }
      }    
      if(!ind){
        covar1 <- rbind(covar[x==1,],covar[x==0,],covar[x==1,],covar[x==0,])
      }else{
        covar1 <- rbind(covar, covar, covar, covar)   
      }
      
      #### example from JD's code
      ##X1 <- model.matrix(~ x1 + z1 + x1*z1 +
      ##                   covar1$age +
      ##                   factor(covar1$lmsepi)
      ##                   )
      ## generalized form TYL
      tmpX1 <- model.matrix(as.formula(paste("~", paste(extra, collapse=" + "))), data=covar1)
      toDrop <- match("(Intercept)", colnames(tmpX1))
      X1 <- cbind(X1, tmpX1[, -toDrop])
      
      
      #### example from JD's code
      ##X  <- model.matrix(~ x + z + x*z +
      ##                   covar$age +
      ##                   factor(covar$lmsepi) +
      ##                   )
      tmpX <- model.matrix(as.formula(paste("~", paste(extra, collapse=" + "))), data=covar)
      toDrop <- match("(Intercept)", colnames(tmpX))
      X <- cbind(X, tmpX[, -toDrop])
      
      glmData <- data.frame(cbind(y, x, z, x*z, data[isPhase2, extra]))
    }
    
    wgt0 <- ifelse(y==1,N1/n1,N0/n0) 
    ##sfit <- glm(y~x+z + x*z+ covar$age + factor(covar$lmsepi) + factor(covar$diabtrt) + factor(covar$hyp) + covar$syst+covar$dias, family=binomial,weights=wgt0)
    options(warn=-1) #remvove message:non-integer #successes in a binomial glm! 
    sfit <- glm(y~., data=glmData, family=binomial, weights=wgt0)
    betastart <- sfit$coef   
    options(warn=0)
    
    diff <- 1
    tol <- 1e-7
    it  <- 0
    
    beta <- betastart
    if(!ind){
      out <- scoreMat(X1,X,y1,x1,y,pzx,beta,Nyx)
      ## with IND assumption (start)
      while (diff > tol & it<=maxit & is.finite(sum(sum(out$cheese))) & !is.na(sum(sum(out$cheese)))) {
        it <- it +1 
         #cat(it,"..")
        out <- scoreMat(X1,X,y1,x1,y,pzx,beta,Nyx)
        if (is.finite(sum(sum(out$cheese))) & !is.na(sum(sum(out$cheese))))
        {
          betanew <- beta + qr.solve(out$cheese,tol=1e-21) %*% out$score
          diff <- max(abs(beta-betanew))
          beta <- betanew
        }else converged=0
      }
      if (it >maxit) {
        converged=0
      }
      if (converged==0)  warning("function diverges!")
      if (converged==1)
      {
        cheese <- out$cheese
        fit1  <- drop(exp(X1 %*% beta)/(1+ exp(X1 %*% beta)))
        pfit1 <- ifelse(y1==1,fit1,1-fit1)
        D1 <- ifelse(y1==1,pfit1*(1-pfit1),-pfit1*(1-pfit1))*X1
        DZW1 <- D1*pzx
        p <- ncol(X1)
        D2 <- matrix(0,length(y1),p)
        D2[x1==0 & y1==0,] <- matrix(rep(apply(DZW1[x1==0 & y1==0,],2,sum),nz0),nrow=nz0,byrow=TRUE)
        D2[x1==1 & y1==0,] <- matrix(rep(apply(DZW1[x1==1 & y1==0,],2,sum),nz1),nrow=nz1,byrow=TRUE)
        D2[x1==0 & y1==1,] <- matrix(rep(apply(DZW1[x1==0 & y1==1,],2,sum),nz0),nrow=nz0,byrow=TRUE)
        D2[x1==1 & y1==1,] <- matrix(rep(apply(DZW1[x1==1 & y1==1,],2,sum),nz1),nrow=nz1,byrow=TRUE)
        pyx <- rep(0,length(y1))
        pyxz<- pfit1*pzx
        pyx[x1==0 & y1==0] <- sum(pyxz[x1==0 & y1==0])
        pyx[x1==1 & y1==0] <- sum(pyxz[x1==1 & y1==0])
        pyx[x1==0 & y1==1] <- sum(pyxz[x1==0 & y1==1])
        pyx[x1==1 & y1==1] <- sum(pyxz[x1==1 & y1==1])
        
        W <- (D1/pyx - D2/(pyx)^2 *pfit1)*Nyx * weight
        
        Wbar <- matrix(0,n,p)
        num <- rep(1:n,2)
        for (i in 1:n) Wbar[i,] <- apply(W[num==i,],2,sum)         
        from.y <- c(y[x==1],y[x==0])
        from.x <- c(x[x==1],x[x==0])   
        
        m1 <- apply( Wbar[from.y==1 & from.x==1,],2,mean)*NXY4/(NXY2+NXY4) + apply( Wbar[from.y==1 & from.x==0,],2,mean)*NXY3/(NXY1+NXY3)
        nxy <- sum(from.y==1 & from.x==1)
        d1 <- NXY4/(NXY2+NXY4)*(Wbar[from.y==1 & from.x==1,] - matrix(m1,nrow=nxy,ncol=p,byrow=TRUE))
        nxy <- sum(from.y==1 & from.x==0)
        d2 <- NXY3/(NXY1+NXY3)*(Wbar[from.y==1 & from.x==0,] - matrix(m1,nrow=nxy,ncol=p,byrow=TRUE))
        I1 <- t(rbind(d1,d2)) %*% (rbind(d1,d2))
        
        m1 <- apply( Wbar[from.y==0 & from.x==1,],2,mean)*NXY2/(NXY2+NXY4) + apply( Wbar[from.y==0 & from.x==0,],2,mean)*NXY1/(NXY1+NXY3)
        nxy <- sum(from.y==0 & from.x==1)
        d3 <- NXY2/(NXY2+NXY4)*(Wbar[from.y==0 & from.x==1,] - matrix(m1,nrow=nxy,ncol=p,byrow=TRUE))
        nxy <- sum(from.y==0 & from.x==0)
        d4 <- NXY1/(NXY1+NXY3)*(Wbar[from.y==0 & from.x==0,] - matrix(m1,nrow=nxy,ncol=p,byrow=TRUE))
        I2 <- t(rbind(d3,d4)) %*% (rbind(d3,d4))
        I <- I1 + I2
        stder <- round(sqrt(diag(qr.solve(cheese,tol=1e-21) %*% (cheese+I) %*% qr.solve(cheese,tol=1e-21))),4)
        
        ##beta <- round(beta,4)
        ## with IND assumption (stop)
      }
      
    }else{
      ## without IND assumption (start)
      it <- 0
      out <- scoreMat(X1,X,y1,x1,y,pzw,beta,Nyx)
      while (diff > tol & it<=maxit & is.finite(sum(sum(out$cheese))) & !is.na(sum(sum(out$cheese)))) {
        it <- it +1 
        out <- scoreMat(X1,X,y1,x1,y,pzw,beta,Nyx)
        if (it<=maxit & is.finite(sum(sum(out$cheese))) & !is.na(sum(sum(out$cheese))))
        {
          betanew <- beta + qr.solve(out$cheese,tol=1e-21) %*% out$score
          diff <- max(abs(beta-betanew))
          beta <- betanew
        }else converged=0
      }
      if (it > maxit) {
        converged=0
      }
      if (converged==0) warning("function diverges!")
      if (converged==1)
      {
        cheese <- out$cheese   
        I3 <- VarDensity(beta,X1,y1,x1,pzw,Nyx,z.from.y,n1,n0,n, N1, N0, N)   
        CovMat <- qr.solve(cheese,tol=1e-21) %*% (cheese+I3) %*% qr.solve(cheese,tol=1e-21)
        ##stder <- round(sqrt(diag(CovMat)),3)
        ##beta <- round(beta,3)
        stder <- sqrt(diag(CovMat))
        ## without IND assumption (start)
      }
    }
    #in divergent cases,generate NAs
    if (converged==0)
    {
      tmpResult <- data.frame(beta <- rep(NA,length(extra)+4), stder=rep(NA,length(extra)+4), pVal=rep(NA,length(extra)+4))
      rownames(tmpResult)[1:4] <- c("(Intercept)","x1","z1","x1:z1")
      rownames(tmpResult)[5:nrow(tmpResult)] <- extra
    }else
    {
      pVal <- 2*(1-pnorm(abs(beta/stder)))
      ##pval1 <- round(2*(1-pnorm(abs(beta[2]/stder[2]))),4)
      ##pval2 <- round(2*(1-pnorm(abs(beta[3]/stder[3]))),4)
      ##pval3 <- round(2*(1-pnorm(abs(beta[4]/stder[4]))),4)
      ##cat(beta[2],"(",stder[2],")","&",pval1, "&", beta[3],"(",stder[3],")","&",pval2, "&",beta[4],"(",stder[4],")","&",pval3,"\n")
      
      tmpResult <- data.frame(beta=round(beta, 4), stder=round(stder, 4), pVal=pVal)
    }
  }
  
  
  rowName <- rownames(tmpResult)
  toRename <- match(c("x1", "z1", "x1:z1"), rowName)
  rowName[toRename] <- c(paste(treatment, "(Treatment)"),
                         paste(BaselineMarker, "(BaselineMarker)"),
                         paste(treatment, BaselineMarker, sep=":")
                         )
  rownames(tmpResult) <- rowName
  return(tmpResult)
}
