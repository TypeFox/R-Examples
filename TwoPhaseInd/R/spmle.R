spmle <- function(data,			## dataset (data.frame)
                  response,			## response variable (binary)
                  treatment,		## treatment variable (binary)
                  BaselineMarker,		## environment variable (continuous)
                  extra=NULL,		## extra variable(s)
                  phase,			## variable for phase indicator
                  ind=TRUE,			## independent or non-indepentent
                  difffactor=0.001,	## 
                  maxit=1000		## max iteration
                  ){
  
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
      coln=NULL
      for (i in 1:length(colNames))
      {
        coln=paste(coln,colNames[i],sep=" ")
      }
      stop(paste0("BaselineMarker variable ",BaselineMarker, " not found in the data:",coln))
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
 
  verbose <- FALSE		## TRUE to turn on debug log
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
  if (remove_rarevariants(data[, BaselineMarker]))
  {
    warnings("BaselineMarker variable is rare variant")
    tmpResult <- data.frame(beta <- rep(NA,length(extra)+4), stder=rep(NA,length(extra)+4), pVal=rep(NA,length(extra)+4))
    rownames(tmpResult)[1] <- "(Intercept)"
    rownames(tmpResult)[5:nrow(tmpResult)] <- extra
  }else
  {
    ## phase 1 
    y <- data[, response]
    x <- data[, treatment]
    z <- data[, BaselineMarker]
    N1 <- sum(y==1)
    N0 <- sum(y==0)
    
    NXYcount <- c(sum(x==0 & y==0),
                  sum(x==1 & y==0),
                  sum(x==0 & y==1),
                  sum(x==1 & y==1)
    )
    if (any(NXYcount==0)) 
    {
      warning(paste0("(x==0 & y==0):",sum(x==0 & y==0),
                     " (x==1 & y==0):",sum(x==1 & y==0),
                     " (x==0 & y==1):",sum(x==0 & y==1),
                     " (x==1 & y==1):",sum(x==1 & y==1)))
      stop("A treatment-response category has zero counts!")
    }
    
    
    ## phase 2
    phase2 <- data[, phase]==2
    phase2Data <- data[phase2, ]
    
    y <- data[phase2, response]
    x <- data[phase2, treatment]
    z <- data[phase2, BaselineMarker]
    n <- length(y)
    n0 <- sum(y==0)
    n1 <- sum(y==1)
    
    wgt0 <- ifelse(y==1,N1/n1,N0/n0)
    
    ##if(!is.null(extra)){
    ##  glmData <- data.frame(cbind(y, x, z, x*z, data[phase2, extra]))
    ##}else{
    ##  glmData <- data.frame(cbind(y, x, z, x*z))
    ##}
    ##
    ##sfit <- glm(y~., data=glmData, family=binomial, weights=wgt0)
    
    if(!is.null(extra)){
      glmFormula <- as.formula(paste(response, "~",
                                     paste(## main model
                                       paste(treatment,
                                             BaselineMarker,
                                             paste(treatment, BaselineMarker, sep=" * "), sep=" + "
                                       ),
                                       ## extra variables
                                       paste(extra, collapse=" + "),
                                       sep=" + "
                                     )
      )
      )
    }else{
      glmFormula <- as.formula(paste(response, "~",
                                     paste(## main model
                                       paste(treatment,
                                             BaselineMarker,
                                             paste(treatment, BaselineMarker, sep=" * "), sep=" + "
                                       ),
                                       sep=" + "
                                     )
      )
      )
    }
    
    options(warn=-1) #remvove message:non-integer #successes in a binomial glm! 	
    sfit <- glm(glmFormula, data=data[phase2, ], family=binomial, weights=wgt0)  
    beta <- sfit$coef
    options(warn=0)
    varmat <- rep(0,length(beta)^2)
    
    if(!is.null(extra)){
      extracov <- 1
      ## the following 2 lines don't work if some variables are factoral
      ##nextracov <- length(extra)	
      ##covvec <- c(t(as.matrix(data[, extra])))
      nextracov <- length(beta) - 4 ## 4 factors
      
      ## re-recognize levels of factors
      phase2Extra <- data.frame(data[phase2, extra])
      colnames(phase2Extra) <- extra
      for(var in extra){
        if(is.factor(phase2Extra[, var])){
          phase2Extra[, var] <- factor(phase2Extra[, var])
        }
      }
      ## dim(phase2Extra)
      
      #### FIXME:  should we constrand on phase2, i.e., data=data[phase2, extra] or data=data[, extra]?
      extraDesMat <- data.frame(model.matrix(as.formula(paste("~", paste(extra, collapse=" + "))), data=phase2Extra), check.names=FALSE)
      extraDesMat <- extraDesMat[, -1, drop=FALSE] ## remove the "(Intercept)" column
      covvec <- c(t(extraDesMat))
      
      ##setdiff(colnames(extraDesMat), names(beta))
      ##setdiff(names(beta), colnames(extraDesMat))
      
      betaName <- c("(Intercept)", treatment, BaselineMarker, paste(treatment, BaselineMarker, sep=":"),
                    colnames(extraDesMat))
      
      
    }else{
      extracov <- 0
      covvec <- 0
      nextracov <- 0
      
      betaName <- c("(Intercept)", treatment, BaselineMarker, paste(treatment, BaselineMarker, sep=":"))
      
    }
    
    ## re-order beta to fit extraDesMat
    beta <- beta[betaName]
    
    tmpResult <- data.frame(beta=rep(NA,length(beta)), stder=rep(NA,length(beta)), pVal=rep(NA,length(beta)))
    #No NA beta
    if (! any(is.na(beta)))
    {
      
      ## decide which C function to call by independent or not
      if(ind){
        nrCode <- "profile_NR_ind"
      }else{
        nrCode <- "profile_NR_noind"
      }
      converged=0
      nrResult <- .C(nrCode,
                     n_subject=as.integer(n),		## Number of subject
                     subj_id=as.integer(1:n),	## Subject ID
                     y=as.double(y),		## Y = Response variable (binary)
                     x=as.double(x),		## X = Treatment (binary)
                     z=as.double(z),		## Z = Biomarker (continuous)
                     NXYcount=as.integer(NXYcount),	## 
                     extracov=as.integer(extracov),
                     covvec=as.double(covvec),
                     n_extracov=as.integer(nextracov),
                     beta=as.double(beta),
                     varmat=as.double(varmat),
                     diff_factor=as.double(difffactor),
                     maxit=as.integer(maxit),
                     verbose=as.integer(sum(verbose)),
                     converged=as.integer(converged),
                     PACKAGE="TwoPhaseInd"
      )
      converged <- nrResult$converged
      if (converged == 0) warnings("The function didn't converge!!!")
      beta <- nrResult$beta
      ## FIXME: CovMat is not symmetric
      CovMat <- matrix(nrResult$varmat, length(beta), length(beta), byrow=TRUE)
      stder <- sqrt(diag(CovMat))
      pVal <- 2*(1-pnorm(abs(beta/stder)))
      
      tmpResult <- data.frame(beta=round(beta, 4), stder=round(stder, 4), pVal=pVal)
    }
    
  }
  
  ## add more information to the row names
  betaName[2] <- paste(treatment, "(Treatment)")
  betaName[3] <- paste(BaselineMarker, "(BaselineMarker)")
  betaName[4] <- paste(treatment, BaselineMarker, sep=":")
  rownames(tmpResult) <- betaName
    ##FIXME
  return(tmpResult)

}
