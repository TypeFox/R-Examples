phase2glm <- function(data,			## dataset (data.frame)
                  response,			## response variable (binary)
                  treatment,		## treatment variable (binary)
                  BaselineMarker,		## environment variable (continuous)
                  extra=NULL,		## extra variable(s)
                  phase			## variable for phase indicator
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
  if (remove_rarevariants(data[, BaselineMarker]))
  {
    warnings("BaselineMarker variable is rare variant")
    tmpResult <- data.frame(beta <- rep(NA,length(extra)+4), stder=rep(NA,length(extra)+4), pVal=rep(NA,length(extra)+4))
    rownames(tmpResult)[1] <- "(Intercept)"
    rownames(tmpResult)[2] <- paste(treatment, "(Treatment)")
    rownames(tmpResult)[3] <- paste(BaselineMarker, "(BaselineMarker)")
    rownames(tmpResult)[4] <- paste(treatment, BaselineMarker, sep=":")
    rownames(tmpResult)[5:nrow(tmpResult)] <- extra
  }else
  {
    ## phase 1 
    y <- data[, response]
    x <- data[, treatment]
    z <- data[, BaselineMarker]
    N1 <- sum(y==1)
    N0 <- sum(y==0)
    
    
    ## phase 2
    phase2 <- data[, phase]==2
    phase2Data <- data[phase2, ]
    
    y <- data[phase2, response]
    x <- data[phase2, treatment]
    z <- data[phase2, BaselineMarker]
    n <- length(y)
    n0 <- sum(y==0)
    n1 <- sum(y==1)
    
    #wgt0 <- ifelse(y==1,N1/n1,N0/n0)
    
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
    
    #options(warn=-1) #remvove message:non-integer #successes in a binomial glm! 	
    sfit <- glm(glmFormula, data=data[phase2, ], family=binomial)  
    beta <- summary(sfit)$coefficients
    beta <- beta[,c(1,2,4)]
    tmpResult <- beta
    
  }
  
 
  return(tmpResult)
}