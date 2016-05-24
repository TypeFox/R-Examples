# (Columwise) multiple imputation via Bayesian Bootstrap
# with Predictive Mean Matching.
# Version:                                0.2
# Date:                            2015-05-17
# Author: F.M.[cre], some contributions: T.S.
# Note:   Needs MASS's stepAIC
# Further infos, references and credits:
#  See for MASS: Venables, W. N. & Ripley, B. D. (2002) Modern Applied Statistics with S. Fourth
#                Edition. New York: Springer.
# License: GPL-2 | GPL-3 (GPL >= 2)


BBPMM <- function(Data,
                  M                = 10,
                  nIter            = 10,
                  outfile          = NULL,
                  ignore           = NULL,
                  vartype          = NULL,
                  stepmod          = "stepAIC",
                  maxit.multi      = 3,  
                  maxit.glm        = 25,   
                  maxPerc          = 0.98,
                  verbose          = TRUE,
                  setSeed          = NULL,
                  chainDiagnostics = TRUE,
                  ...)
{
  #### General declarations I ####
  stepwise <- stepmod == "stepAIC"
  DAT      <- as.data.frame(Data)
  orgnames <- varnames <- attributes(DAT)$names
  n        <- nrow(DAT) 
  l        <- length(orgnames)
  onevar   <- FALSE
  lv       <- 1:l
  nv       <- 1:n
  fun_call <- match.call() 
    
  #### General warnings & stops I ####
  if(n < 200) warning("Small data sets can reduce the quality of the predictive mean match.", immediate.=T)
  
  ### Test old arguments
  dot_arg <- names(sapply(substitute(list(...))[-1],deparse))
  if(any(dot_arg=="maxit")) warning("The argument maxit is deprecated. Please use maxit.multi and maxit.glm instead.")
  if(any(dot_arg=="stepwise")) warning("The argument stepwise is deprecated. Please use stepmod instead.")
  
  ### Test for zero variance
  pos_novar <- apply((asNumericMatrix(DAT)),2,var,na.rm=T)==0
  if(sum(pos_novar) > 0){
    warning(paste("Variable",which(pos_novar),"contains no variation and will be ignored. "))
    if(!is.null(ignore)){
      if(is.character(ignore)){
        ignore <- append(ignore,names(which(pos_novar)))
      } else {
        ignore <- append(ignore,which(pos_novar))
      }
    }
  }
  
  ### Test if ignore is set
  if (!is.null(ignore)) {
    if((length(setdiff(ignore,varnames)) > 0 & is.character(ignore)) |
         (length(setdiff(ignore,lv)) > 0 & is.numeric(ignore))) {
      stop("'ignore' is no subset of either number of columns or variable names!\n")
    }
    if (is.character(ignore)) {
      origR  <- is.na(DAT)  
      ig.pos <- varnames %in% ignore
      ignore <- which(ig.pos)
    } else {
      origR    <- is.na(DAT)  
      ig.pos <- lv %in% ignore
      ignore <- which(ig.pos)
    }
    
    #### General declarations II (concerning ignored Variables) ####
    org.l     <- l
    not.inc   <- as.data.frame(DAT[ ,ignore])
    varnames  <- varnames[-ignore] 
    DAT       <- DAT[ ,-ignore]
    l         <- ncol(DAT)
    lv        <- 1:l
  }
    
  if (l <= 2) {
    ## add up to two columns with random noise to ensure variability for PMM
    DAT      <- as.data.frame(cbind(DAT,matrix(runif((3-l)*n),nrow=n)))
    varnames <- attributes(DAT)$names
    onevar   <- TRUE
  }
  
  ## take over class from data.frame or administer classes
  if (!is.null(vartype)) {
    
    #### General warnings & stops III (concerning vartype) ####
    if (length(vartype) != l) {
      stop("Error: Number of flagged variables in 'vartype'",
           "does not match number of (remaining) variables in",
           "data set!\n")
    } else if (any(vartype != "C" & vartype != "M")) {
      stop("Error: 'vartype' contains wrong character(s)!\n")
    }
    
    
    f.pos <- vartype == "C"
    if (sum(f.pos) > 0) {
      DAT[,f.pos] <- lapply(DAT[,f.pos], as.factor)}
    
    m.pos <- vartype == "M"
    if (sum(m.pos) > 0) {
      DAT[,m.pos] <- lapply(DAT[,m.pos], as.numeric)}
  }
  
  #### General declarations III ####
  M.DAT <- vector(M, mode = "list")
  M.DAT <- lapply(M.DAT,function(x){'<-'(x,data.frame(matrix(ncol=l,nrow=n)))})
  
  if (!is.null(setSeed)) set.seed(setSeed)
  firstSeed <- .Random.seed 
  nM    <- "M"
  while (any(varnames == nM)) nM <- paste0(nM,"M")
  
  ## autoselect iterations
  if(nIter=="autolin"){
    nIter <- round(19*(1-dmi(DAT))+1)
    print(paste("autolin-Iterations used. nIter set to:",nIter))
  }
    
  ## indicator matrix for missing values
  #Definitions & Declarations IV
  R          	<- matrix()
  length(R)   <- n*l
  dim(R)	     <- c(n,l)
  R           <- is.na(DAT)
  mis.num     <- colSums(R)
    
  mis.overview <- paste0("number of missing values ", attributes(DAT)$names,": ",
                         mis.num)
  if (verbose) print(mis.overview)
  
  if(sum(mis.num)==0) stop(paste("No missing values found."))
  
  ## new variable order
  n.order  <- order(mis.num)
  o.order  <- order(n.order)
  DAT      <- as.data.frame(DAT[ ,n.order])
  varnames <- varnames[n.order]
  mvar     <- mis.num[n.order] != 0
  p.impvar <- which(mvar)
  p.comp   <- which(!mvar)
  i.mis    <- i.obs <- matrix(ncol=sum(mvar),nrow=n)
  
  i.mis <- as.matrix((R[,n.order])[,mvar])
  i.obs <- as.matrix(!i.mis)
  
  imh <- integer(n/2)
  ioh <- integer(n/2)
  
  ### chainDiagnostics
  saveChain <- NULL
  
  if(chainDiagnostics){ #New
    
    ## Create a list for every variable
    saveChain <- vector("list",length(p.impvar))
    
    for(chains in 1:length(saveChain)){
      
      ## Each listelement contains a list for the imputation run
      saveChain[[chains]] <- vector("list",M)
      
      ## ... containing a matrix saving the results for the 
      ## imputed values per iteration
      for(impss in 1:M){
        saveChain[[chains]][[impss]] <-  matrix(NA,sum(i.mis[,chains]),nIter)
        if(impss == 1){
          colnames(saveChain[[chains]][[impss]]) <- c("StartSol",paste0("Iter",2:nIter))
        } else {
          colnames(saveChain[[chains]][[impss]]) <- paste0("Iter",1:nIter)
        }
        
      }
      #Give names for imputation elements
      names(saveChain[[chains]]) <- paste0("M",1:M)  
    }
    #Give names for variable elements 
    names(saveChain) <- colnames(i.mis)
  }
  
  #### General warnings & stops IV (concerning NAs) ####
  if(any(n == mis.num)) stop(paste("Column",
                                   which(n == mis.num),
                                   "contains solely missing values.\n Please remove from data set before imputation."))
  
  ## starting solution
  startSol <- TRUE
  ##+++++++++++++++++++++++++PMM++++++++++++++++++++++++++++++++++++
  ## Sequential Regression with Predictive Mean Matching
  for (m in 1:M) {
    co      <- 0
    iterate <- TRUE
    
    
    while (iterate) {
      ##--------------first loop for iterations-----------------------
      if (!startSol) co <- co + 1
      co2 <- 0
      if (verbose & !startSol) {
        cat(paste0("Imputation ", m," of ",M ,": iteration ", co), "\n") }
      ##-------------- Bayesian Bootstrap --------------------------
      if (M > 1) {
        ind1 <- BayesBoot(ind.obs = nv)
        ## Bayesian Bootstrap: draw n times with replacement as basis for
        ## imputation model parameter estimates
      }
      if (verbose & !startSol) cat("Variable:")
      ##----second loop for every variable with missing values------
      for (j in p.impvar) {
        if (verbose & !startSol) {cat("",varnames[j])
                                  if (j == rev(p.impvar)[1]) cat("\n ")
        }
        co2 <- co2 + 1
        y   <- DAT[ ,j]
        ioh <- which(i.obs[,co2])
        imh <- which(i.mis[,co2])
        if (startSol & length(p.comp) == 0) {
          ## hotdeck imputation for starting solution if no variable
          ## completely observed
          DAT[imh,j] <- sample(y[ioh], length(imh), replace = TRUE)
          p.comp     <- j
          y          <- DAT[ ,j]
        }
        
        if(!startSol){
          xvars <- paste(setdiff(varnames[-j],varnames[noV]), collapse = ' + ')
        } else {
          
          noV   <- which(sapply(DAT[ioh,],var,na.rm=TRUE) == 0) #old: which(apply(DAT[ioh,],2,var,na.rm=TRUE) == 0) 
          xvars <- paste(setdiff(c(varnames[p.comp], varnames[p.impvar[0:(co2-1)]]),
                                 varnames[noV]), collapse=' + ')
        }
        
        s.model <- as.formula(paste0(varnames[j],' ~ ',xvars))
        ## hotdeck imputation if y (almost) Dirac distributed
        if (!startSol) {usePos <- nv}  else {usePos <- ioh}
        
        if (var(y[usePos],na.rm=TRUE) == 0 | max(table(y))/length(na.omit(usePos)) > maxPerc) {
          DAT[imh,j] <- sample(y[usePos], length(imh),TRUE)
          next
        }
        if (is.numeric(y)) {
          if (M > 1 & !startSol) {
            while (var(y[ind1],na.rm=TRUE) == 0) {
              ind1 <- BayesBoot(ind.obs = nv)
            }
            
            BB.data <- DAT[ind1, ]
            
            BB.stab <- BB.mod.stab.glm(data=DAT, BB.data=BB.data,
                                       s.model=s.model, maxit.glm=maxit.glm)
            
            regmod  <- BB.stab$model
          }
          else if (M == 1 | startSol) {
            regmod <- lm(s.model, data=DAT, subset = ioh)
          } 
          
          
          if (stepwise & summary(regmod)$r.squared < 1) {
            try(regmod <- stepAIC(regmod, trace=0, k=log(n), direction = "backward"))
          } #Expand here for different selection
          
          y.pred     <- predict(regmod, newdata=DAT)
          y.pred.mis <- y.pred[imh]
          y.pred.obs <- y.pred[ioh]
          nextlist   <- numeric(length(y.pred.mis))
          for (i in seq_along(y.pred.mis)) {
            nextlist[i] <- PMMsearchMet(yHatMis = y.pred.mis[i],
                                        yHatObs = y.pred.obs)
          }
        } else if (is.factor(y) & length(table(y)) > 2) {
          if ( M > 1 & !startSol) {
            BB.data <- DAT[ind1, ]
            
            BB.stab <- BB.mod.stab.mlog(data=DAT,
                                        BB.data=BB.data,
                                        s.model=s.model,maxit.multi=maxit.multi)
            
            regmod  <- BB.stab$model
          } else if (M == 1 | startSol) {
            options(warn = -1)
            regmod <- multinom(s.model,data=DAT,trace=F,
                               subset = ioh, ...=list(maxit=maxit.multi))
            options(warn = 0)
          } 
          
          if (stepwise) try(regmod <- stepAIC(regmod, trace=0, k=log(n), direction = "backward"))
          
          y.pred                 <- predict(regmod, newdata=DAT, type="probs")
          y.pred[y.pred > 0.999] <- 0.999
          y.pred[y.pred < 0.001] <- 0.001
          l.y.pred               <- log(y.pred/(1-y.pred))
          y.pred.mis             <- l.y.pred[imh, ]
          y.pred.obs             <- l.y.pred[ioh, ]
          
          ## calculate outer product for all obs/mis columns
          
          
          nr.pred.mis <- nrow(y.pred.mis)                                                                                  
          nr.pred.obs <- nrow(y.pred.obs)                                           
          
          #Two scalars
          if(is.null(nr.pred.mis) && is.null(nr.pred.obs)){                  
            m.dist                 <- matrix(0,1,1)                         
            
            for (jj in 1:ncol(y.pred)){                                       
              m.dist <- m.dist+(y.pred.mis[jj]-y.pred.obs[jj])^2                     
            }                                                                  
            
          } else if (is.null(nr.pred.mis)){       #if mis is scalar            
            m.dist                 <- matrix(rep(0,nr.pred.obs),               
                                             nrow=1)                           
            
            for (jj in 1:ncol(y.pred)){                                        
              m.dist <- m.dist+(y.pred.mis[jj]-y.pred.obs[ ,jj])^2 
            }                                                                  
            
          } else if (is.null(nr.pred.obs)){   #if obs is scalar                
            m.dist                 <- matrix(rep(0,nr.pred.mis),              
                                             nrow=nr.pred.mis)                 
            
            for (jj in 1:ncol(y.pred)){                                        
              m.dist <- m.dist+outer(y.pred.mis[ ,jj],y.pred.obs[jj],"-")^2        
            }                                                                  
            
          } else {                          #if both are vectors               
            
            m.dist                 <- matrix(rep(0,nrow(y.pred.mis)*nrow(y.pred.obs)),  
                                             nrow=nrow(y.pred.mis))                     
            
            for (jj in 1:ncol(y.pred)){                                                 
              m.dist <- m.dist+outer(y.pred.mis[ ,jj],y.pred.obs[ ,jj],FUN="-")^2      
            }                                                                          
            
          }                                                                    
          
          
          nextlist <- max.col(as.matrix(m.dist*(-1)),
                              ties.method="random")
          rm("m.dist")
          
        } else if (is.factor(y) & length(table(y)) == 2) {
          if (M > 1 & !startSol) {
            BB.data <- DAT[ind1, ]
            
            BB.stab <- BB.mod.stab.glm(data=DAT, BB.data=BB.data,
                                       s.model=s.model, model="binomial",maxit.glm=maxit.glm)
            
            
            regmod  <- BB.stab$model
          } else if (M == 1 | startSol)  {
            
            regmod <- glm(s.model, data=DAT,
                          family = binomial(link="logit"),
                          subset = ioh,
                          control=glm.control(maxit = maxit.glm))
            
            
          } 
          
          if (stepwise) try(regmod <- stepAIC(regmod, trace=0, k=log(n), direction = "backward"))
          y.pred     <- predict(regmod, newdata=DAT)
          y.pred.mis <- y.pred[imh]
          y.pred.obs <- y.pred[ioh]
          nextlist   <- numeric(length(y.pred.mis))
          for (i in seq_along(y.pred.mis)) {
            nextlist[i] <- PMMsearchMet(yHatMis = y.pred.mis[i],
                                        yHatObs = y.pred.obs)
          }
        }
                
        DAT[imh,j] <- y[ioh][nextlist]
        
        # Save Chain when chainDiagnostics is TRUE
        if(chainDiagnostics){
          saveChain[[which(p.impvar == j)]][[m]][1:length(imh),co] <- y[ioh][nextlist]
        } 

      }
      startSol <- FALSE
      if (co == nIter) iterate <- FALSE
    } ## end of iterate cycle
    if (onevar == TRUE) {
      M.DAT[[m]] <- as.data.frame(DAT[,lv])
      attributes(M.DAT[[m]])$names <- orgnames
    }
    if (is.null(ignore)) {
      M.DAT[[m]] <- DAT[ ,o.order]
    } else {
      M.DAT[[m]]            <- as.data.frame(matrix(nrow=n,ncol=org.l))
      M.DAT[[m]][ ,!ig.pos] <- DAT[ ,o.order]
      M.DAT[[m]][ ,ig.pos]  <- not.inc
      attributes(M.DAT[[m]])$names     <- orgnames
    }
    if (M == 1) M.DAT <- M.DAT[[m]]
  }
  ###################### end of m cylce ################################
  if (!is.null(outfile)) {
    outDAT <- matrix(nrow=n*M, ncol=length(orgnames)+1)
    if (M == 1) {
      outDAT <- cbind(M.DAT,1)
    } else {
      for (i in seq_along(M.DAT)) {
        outDAT[((i-1)*n+1):(i*n), ] <- cbind(as.matrix(M.DAT[[i]]),i)
      }
    }
    outDAT <- as.data.frame(outDAT)
    names(outDAT) <- c(orgnames,nM)
    ## check for file ending
    if (grep(".",outfile) > 0) {
      lastDot <- max(which(strsplit(outfile,"")[[1]]=="."))
      StrL <- length(strsplit(outfile,"")[[1]])
      if ((StrL - lastDot) > 3) {
        outfile <- paste0(outfile,".dat")
      }
    } else {
      outfile <- paste0(outfile,".dat")
    }
    write.table(outDAT, file = outfile, sep = "\t",
                row.names = FALSE, quote = FALSE)
  }
  
  if (!is.null(ignore)) R <- origR 
  
  x <- list("call"             = fun_call,
            "mis.num"          = colSums(R),
            "modelselection"   = stepmod,
            "seed"             = setSeed,
            "impdata"          = M.DAT,
            "misOverview"      = mis.overview,
            "indMatrix"        = R,
            "M"                = M,
            "nIter"            = nIter,
            "Chains"           = saveChain,
            "FirstSeed"        = firstSeed,
            "LastSeed"         = .Random.seed,
			         "ignoredvariables" = !is.null(ignore))
  class(x) <- "imp"
  return(x)
}