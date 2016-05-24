# (Rowwise) multiple imputation via Bayesian Bootstrap
# with Predictive Mean Matching.
# Version:                           0.2
# Date:                       2015-04-27
# Author: F.M., some contributions: T.S.
# Note:   Needs MASS's stepAIC
# Further infos, references and credits:
#  See for MASS: Venables, W. N. & Ripley, B. D. (2002) Modern Applied Statistics with S. Fourth
#                Edition. New York: Springer.
# License: GPL-2 | GPL-3 (GPL >= 2)

BBPMM.row <- function(misDataPat,
                      blockImp   = length(misDataPat$blocks),
                      M          = 10,
                      outfile    = NULL,
                      manWeights = NULL,
                      stepmod    = "stepAIC",
                      verbose    = TRUE,
                      tol        = 0.25,
                      setSeed    = NULL,
                      ...)
{
  
  
  #### General warnings & stops I ####
  ### Testing data preperation
  if(class(misDataPat) != "impprep"){
    warning("Data preperation needed. Trying rowimpPrep() on given data.",immediate.=T)
    misDataPat <- rowimpPrep(misDataPat)
    blockImp   <- length(misDataPat$blocks)
  }
  
  ### Testing old arguments
  
  dot_arg <- names(sapply(substitute(list(...))[-1],deparse))
  if (any(dot_arg=="stepwise")) warning("The argument stepwise is deprecated. Please use stepmod instead.")
  
  ### Testing blocks
  if(length(setdiff(blockImp, 1:length(misDataPat$blocks))) > 0) {
    stop(paste("blockImp =",as.character(setdiff(blockImp, 1:length(misDataPat$blocks))),
               "is not a subset of the number of",
               "different missing-data patterns (blocks)!\n"))}
  
  if(length(misDataPat$blocks) > length(blockImp)){
    warning("Per default only the last block is imputed")
  }
  
  
  
  #### General declarations I ####
  fun_call   <- match.call() 
  b.order    <- order(blockImp)
  blockImp   <- sort(blockImp)
  
  
  #### General warnings & stops II ####
  if(!is.null(manWeights)) {
    if(!is.vector(manWeights, mode = "list")) {
      manWeights <- list(manWeights)
      if (length(blockImp) > 1) {
        stop(paste("Only one vector with manual weights, but more than one",
                   "block specified for imputation!\n"))}
    } else if (is.vector(manWeights, mode = "list") & (length(blockImp) < length(manWeights))) {
      stop(paste("'manWeight' contains more elements than 'blockImp'!\n"))
    }
    if (any(unlist(manWeights) < 0)) {
      stop(paste("'manWeight' contains negative value(s)!\n"))}
    if (length(b.order) > 1) manWeights <- manWeights[b.order]
  }
  
  
  #### General Declarations II ####
  
  data.set   <- misDataPat$data
  n          <- nrow(data.set)
  l          <- ncol(data.set)
  R          <- is.na(data.set)
  nmis       <- apply(R, 2, sum) 
  varnames   <- names(data.set)
  ignore       <- misDataPat$ignore
  ignored_data <- misDataPat$ignored_data
  
  key        <- misDataPat$key
  block      <- misDataPat$blocks[blockImp]
  mis.pos    <- obs.pos <- s.model <- vector(mode="list", length=length(block))
  comp.names <- misDataPat$compNames
  stepwise   <- stepmod == "stepAIC"
  weight.matrix <- model <- 
    y.hat <- impdata <- vector(mode="list", length=M)
  
  names(weight.matrix) <- names(model) <- 
    names(y.hat) <- names(impdata) <-
    paste("M",1:M,sep="")
  
  
  if (!is.null(key)) {
    pairlst <- vector(mode="list", length=M)
    names(pairlst) <- paste("M",1:M,sep="") 
  }
  
  miss <- function(x) {any(is.na(x)) }
  
  if (!is.null(key) && ncol(key) == 1) {
    Donid <- Recid <- names(key)[1]
  } else if (!is.null(key) && ncol(key) == 2) {
    Recid <- names(key)[1]
    Donid <- names(key)[2]
  } else if (is.null(key)) {
    pairlst <- NULL}
  
  if(!is.null(setSeed)) set.seed(setSeed)
  firstSeed <- .Random.seed
  
  nYC     <- "YC"
  nM      <- "M"
  
  while (any(varnames == nYC)) nYC <- paste(nYC,"C",sep="")
  while (any(varnames == nM)) nM   <- paste(nM,"M",sep="")
  
  #########################################################################
  mrow <- logical(n)
  for (j in seq(along=block)) {
    mrow <- apply(as.matrix(data.set[ ,block[[j]]]), 1, miss)
    mis.pos[[j]] <- (1:n)[mrow] 
    obs.pos[[j]] <- (1:n)[!mrow]
    ## Test for available degress of freedom in the model
    
    
    #Change for qr, T.S.
    mc.test <- data.set[obs.pos[[j]], comp.names]
	   mc.test <- matrix(as.numeric(unlist( mc.test)),NROW( mc.test))
    mc.test <- qr(mc.test,...)$rank 
    
    if((((length(obs.pos[[j]])-1) <= length(comp.names)) & (stepmod != "stepAIC")) |
         (mc.test != length(comp.names) & (stepmod != "stepAIC"))) {
      stop("Block ",j, " has insufficient rank for imputation!\n")
    }
  }
  
  
  #### General Declarations III ####
  
  xvars <- paste(comp.names,collapse= ' + ') 
  weight.matrix <- lapply(weight.matrix,
                          function(x){'<-'(x,vector(mode = "list",
                                                    length=length(block)))})
  model   <- lapply(model,
                    function(x){'<-'(x,vector(mode = "list",
                                              length=length(block)))})
  y.hat   <- lapply(y.hat,
                    function(x){'<-'(x,vector(mode = "list",
                                              length=length(block)))})
  impdata <- lapply(impdata,
                    function(x){'<-'(x,as.data.frame(matrix(nrow=n,ncol=l)))})
  
  #### General warnings & stops III ####
  if(n < 200) warning("Small data sets can reduce the quality of the predictive mean match.", immediate.=T)
  
  if(any(n == nmis)) stop(paste("Column", which(n == nmis), "contains solely missing values.\n Please remove from data set before imputation."))
  
  ##########################################################################
  ### first loop for MI----------------------------------------------------
  for (m in 1:M) {
    if(!is.null(key)) {
      pairlst[[m]] <- vector(mode = "list", length=length(block))
      names(pairlst[[m]]) <- paste("block",seq(along=block),sep="")
    }
    
    names(weight.matrix[[m]]) <- names(model[[m]]) <- names(y.hat[[m]]) <- 
      paste("block",seq(along=block),sep="")
    
    for (j in seq(along=block)) { # second loop for different blocks
      model[[m]][[j]] <- vector(mode="list", length=length(block[[j]]))
      names(model[[m]][[j]]) <- varnames[block[[j]]]
      S.xy <- NULL
      y.hat[[m]][[j]] <- matrix(nrow=n,ncol=length(block[[j]]))
      colnames(y.hat[[m]][[j]]) <- varnames[block[[j]]]
      co2 <- 0
      if(M > 1){
        BB.ind  <- BayesBoot(ind.obs = obs.pos[[j]])
        BB.data <- data.set[BB.ind, ]
      }
      for (k in block[[j]]) {
        co2 <- co2+1
        s.model <- as.formula(ifelse(co2 == 1,
                                     paste(varnames[k],'~',xvars),
                                     paste(nYC,'~',xvars)))
        if (co2 == 1) {
          if(M == 1) {
            regmod <- lm(s.model, data=data.set, subset=obs.pos[[j]])
          } else if(M > 1) {
            BB.stab <- BB.mod.stab.glm(data=data.set,BB.data=BB.data,
                                       s.model=s.model)
            regmod <- BB.stab$model
            if (any(BB.stab$mislevpos == TRUE)) {
              warning(paste("Imputation ",m,", block ",j,
                            ": Bayesian Bootstrap dropped ",
                            "at least one category of a factor variable!\n",
                            sep=""))
            }
          }
          
          if (stepwise){ regmod <- stepAIC(regmod, trace=0, k=log(n),
                                           direction = "backward")
          }
          
          var.T <- var(data.set[ ,k], na.rm=TRUE)
          var.U <- var(regmod$residuals)*length(obs.pos[[j]])/
            (length(obs.pos[[j]])-regmod$rank+1)
        } else if (co2 > 1) {
          xvars.e <- paste(varnames[block[[j]]][1:(co2-1)], collapse='+')
          s.model.e <- as.formula(paste(varnames[k],' ~ ',xvars.e))
          if(M == 1) {
            regmod.e <- lm(s.model.e, data=data.set, subset=obs.pos[[j]])
          } else if(M > 1) {
            regmod.e <- lm(s.model.e, data=BB.data)
          }
          YC <- numeric(n)
          YC[obs.pos[[j]]] <- regmod.e$residuals
          data.set <- as.data.frame(cbind(YC,data.set))
          names(data.set)[1] <- nYC
          if(M == 1) { 
            regmod <- lm(s.model, data=data.set, subset=obs.pos[[j]])
          } else if(M > 1) {
            BB.data <- as.data.frame(cbind(YC[BB.ind],BB.data))
            names(BB.data)[1] <- nYC
            BB.stab <- BB.mod.stab.glm(data=data.set,BB.data=BB.data,
                                       s.model=s.model)
            regmod <- BB.stab$model
          }
          
          
          if (stepwise){ regmod <- stepAIC(regmod, trace=0, k=log(n),
                                           direction = "backward")
          }
          
          var.T <- var(YC, na.rm=TRUE)
          var.U <- var(regmod$residuals)*length(obs.pos[[j]])/
            (length(obs.pos[[j]])-regmod$rank+1)
        }
        model[[m]][[j]][[co2]] <- regmod 
        y.hat[[m]][[j]][ ,co2] <- predict(regmod,newdata=data.set)
        ## multicollinearity among ys
        if (is.na(var.T) || var.T < 1e-16) {
          S.xy[co2] <- 1e16 
        } else {
          S.xy[co2] <- var.U}
        if (co2 > 1) { data.set <- data.set[ ,-1]
                       if (M > 1) { BB.data <- BB.data[ ,-1] }}
        
      } ## end of loop k (incomplete variables)
      suppressWarnings(rm("BB.data"))
      gc()
      
      if (length(S.xy) > 1) {
        weight.matrix[[m]][[j]] <- diag(S.xy)
      } else {
        weight.matrix[[m]][[j]] <- S.xy
      }
   
      if (!is.null(manWeights) && length(manWeights[[j]]) > 0) {
        if (length(manWeights[[j]]) > 1) {
          weight.matrix[[m]][[j]] <- diag(manWeights[[j]]^(-1)*
                                            diag(weight.matrix[[m]][[j]]))
        } else if (length(manWeights[[j]]) == 1) {
          weight.matrix[[m]][[j]] <- weight.matrix[[m]][[j]]/manWeights[[j]]
        }
      }
      
      #just solve version 
      weight.matrix[[m]][[j]] <- solve(weight.matrix[[m]][[j]])
      
      y.hat.obs <- y.hat[[m]][[j]][obs.pos[[j]], ]
      
      if (verbose) {
        cat(paste("Imputation ",m,": reciprocal weight matrix for block ",j,
                  ":\n",sep = ""))
        print(weight.matrix[[m]][[j]])
      }
      if (!is.null(key)) {
        pairlst[[m]][[j]] <- matrix(nrow=length(mis.pos[[j]]),ncol=2)}
      co3 <- 0
      for (i in mis.pos[[j]]) # third loop b) for the unobserved ys
      {
        co3 <- co3+1
        
        
        #### PMM-Distance ####
        
        pmmdist <- PMMC(y.hat[[m]][[j]][i, ],
                        as.matrix(weight.matrix[[m]][[j]]),t(y.hat.obs))
        
        pmmdist <- pmmdist == min(pmmdist)
        index   <- obs.pos[[j]][pmmdist]
        
        if (length(index) > 1) {
          index <- sample(index, 1)
        } # random selection in case of several nearest neighbours
        if (!is.null(key)) {
          pairlst[[m]][[j]][co3, ] <- c(key[i,Recid],key[index,Donid])
        }
        impdata[[m]][i, block[[j]]] <- data.set[index, block[[j]]]
      } ## end of i loop (missing values)
      
    } ## end of j loop (blocks)
    if (!is.null(key)) { impdata[[m]] <- cbind.data.frame(key, impdata[[m]]) }
    
    
    
  } ## end of m loop (MI)
  
  #### Prepare Output ####
  if (!is.null(key)){
    for( i in 1:length(impdata)){	
      for(j in 1:ncol(R)){
        impdata[[i]][,-(1:ncol(key))][!R[,j],j] <- data.set[!R[,j],j]
      }		
      names(impdata[[i]]) <- c(names(key),names(data.set))		
    }
  } else {
    for( i in 1:length(impdata)){	
      for(j in 1:ncol(R)){
        impdata[[i]][!R[,j],j] <- data.set[!R[,j],j]
      }	
      names(impdata[[i]]) <- names(data.set)
    }
  }
  
  if(!is.null(ignore)){
    for(m in 1:M){
      for(ipn in 1:length(ignore)){
        impdata[[m]] <- inbind(impdata[[m]],ignore[ipn]-1,ignored_data[ipn])
      }
    }
    if(is.null(key)){
      l          <- ncol(impdata[[1]]) - ncol(key)
    } else {
      l          <- ncol(impdata[[1]]) 
    }
  }
  
  
  if (M == 1) impdata <- impdata[[1]]
  if (!is.null(outfile)) {
    if(is.null(key)){
      outDAT <- matrix(nrow=n*M, ncol=l+1)
    } else{
      outDAT <- matrix(nrow=n*M, ncol=l+ncol(key)+1)
    }
    
    if (M == 1) {
      outDAT <- cbind(impdata,1)
    } else {
      for (i in seq(along = impdata)) {
        outDAT[((i-1)*n+1):(i*n), ] <- cbind(as.matrix(impdata[[i]]),i)
      }
    }
    outDAT <- as.data.frame(outDAT)
    if(is.null(key)){
      names(outDAT) <- c(varnames,nM) ## check for file ending
    } else {
      names(outDAT) <- c(names(key),varnames,nM) ## check for file ending
    }
    if (grep(".",outfile) > 0) {
      lastDot <- max(which(strsplit(outfile,"")[[1]]=="."))
      StrL <- length(strsplit(outfile,"")[[1]])
      if ((StrL - lastDot) > 3) {
        outfile <- paste(outfile,".dat",sep="")
      }
    } else {
      outfile <- paste(outfile,".dat",sep="")
    }
    write.table(outDAT, file = outfile, sep = "\t",
                row.names = FALSE, quote = FALSE)
  }
  
  if(!is.null(key) | !is.null(ignored_data)){
  	R <- misDataPat$indMatrix
  }
  
  
  x <- list("call"             = fun_call,
            "mis.num"          = nmis,
            "modelselection"   = stepmod,
            "seed"             = setSeed,
            "impdata"          = impdata,
            "indMatrix"        = R,
            "M"                = M,
            "weightMatrix"     = weight.matrix,
            "model"            = model,
            "pairlst"          = pairlst,
            "FirstSeed"        = firstSeed,
            "LastSeed"         = .Random.seed,
         			"ignoredvariables" = !is.null(ignore))
  class(x) <- "imp"
  return(x)
}