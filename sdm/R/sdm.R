# Author: Babak Naimi, naimi.b@gmail.com
# Date :  April 2016
# Version 2.4
# Licence GPL v3
#--------


.methodFix <- function(n) {
  for (i in seq_along(n)) {
    nx <- .sdmMethods$whichMethod(n[i])
    if (!is.null(nx)) n[i] <- nx
    else n[i] <- NA
  }
  n
}
#----------
.replicate.methodFix <- function(n) {
  for (i in seq_along(n)) {
    nx <- .replicateMethods$whichMethod(n[i])
    if (!is.null(nx)) n[i] <- nx
    else n[i] <- NA
  }
  n
}
#----------
.getSpeciesDistribution <- function(data) {
  o <- lapply(data@species,function(x) {
    if (!is.null(x@presence)) return('binomial')
    else if (!is.null(x@abundance)) return('poisson')
    else if (!is.null(x@Multinomial)) return('multinomial')
    else return(NA)
  })
  n <- names(o)
  o <- as.character(o)
  names(o) <- n
  o
}
#-------------

.getFormula.rhs <- function(n,env=parent.frame()) {
  as.formula(paste('~',paste(n,collapse='+'),sep=''),env = env)
}

.getFormula <- function(n,env=parent.frame()) {
  as.formula(paste(n[1],'~',paste(n[-1],collapse='+'),sep=''),env = env)
}

.getFormula.gammgcv.rhs <- function(n,nFact=NULL,env=parent.frame()) {
  if (!is.null(nFact)) as.formula(paste('~',paste(c(paste(paste('s(',n,sep=''),')',sep=''),nFact),collapse='+'),sep=''),env = env)
  else as.formula(paste('~',paste(paste(paste('s(',n,sep=''),')',sep=''),collapse='+'),sep=''),env = env)
}
.getFormula.gammgcv <- function(n,nFact=NULL,env=parent.frame()) {
  if (!is.null(nFact)) as.formula(paste(n[1],'~',paste(c(paste(paste('s(',n[-1],sep=''),')',sep=''),nFact),collapse='+'),sep=''),env = env)
  else as.formula(paste(n[1],'~',paste(paste(paste('s(',n[-1],sep=''),')',sep=''),collapse='+'),sep=''),env = env)
}
#----------

.addLHSformula <- function(f,lhs,env=parent.frame()) {
  # ~ x1 + x2; add lhs to such formula
  as.formula(paste(lhs,'~',as.character(f)[2]),env = env )
} 


#----------------

.mahal <- function(d1,d2) {
  co <- solve(cov(d1))
  mahalanobis(d2,colMeans(d1,na.rm=TRUE),co,inverted=TRUE)
}
#----------
.checkFactor <- function(f1,f2) {
  f2 <- factor(f2)
  f1 <- factor(f1)
  l <- levels(f2)[!levels(f2) %in% levels(f1)]
  ww <- NULL
  if (length(l) > 0) {
    w1 <- w2 <- ww <- list()
    for (i in seq_along(l)) w1[[l[i]]] <- which(f2 == l[i])
    l2 <- levels(f2)[levels(f2) %in% levels(f1)]
    for (i in seq_along(l2)) w2[[l2[i]]] <- which(f2 == l2[i])
    ww[['p']] <- w1
    ww[['np']] <- w2
  }
  ww
}


#----------------
.eqFactLevel <- function(data1,data2) {
  nFact <- colnames(data2)[.where(is.factor,data2)]
  for (i in seq_along(nFact)) data2[,nFact[i]] <- factor(data2[,nFact[i]],levels = levels(data1[,nFact[i]]))
  data2
}
#-------------

.factorFix <- function(data1,data2,nFact,nf) {
  # assign the problematic factors to the more similar factors according to continuous variables
  # if no continus variable does exist, then it is assigned to a dominant class.
  if (missing(nFact) || is.null(nFact)) nFact <- colnames(data2)[.where(is.factor,data2)]
  if (missing(nf) || is.null(nf)) {
    nf <- colnames(data2)[!colnames(data2) %in% nFact]
    if (length(nf) == 0) nf <- NULL
  }
  
  dd <- data2
  
  for (i in seq_along(nFact)) {
    data1[,nFact[i]] <- factor(data1[,nFact[i]])
    data2[,nFact[i]] <- factor(data2[,nFact[i]])
    dd[,nFact[i]] <- as.character(dd[,nFact[i]])
    
    fc <- .checkFactor(data1[,nFact[i]],data2[,nFact[i]])
    if (!is.null(fc)) {
      p <- names(fc$p)
      np <- names(fc$np)
      for (j in seq_along(p)) {
        w <- fc[['p']][[p[[j]]]]
        if (!is.null(nf)) {
          m <- rep(NA,length(np))
          d2 <- data2[w,which(colnames(data2) %in% nf)]
          if (length(nf) > 1) {
            options(warn=-1)
            for (k in seq_along(np)) {
              ww <- fc[['np']][[np[[k]]]]
              d1 <- data2[ww,which(colnames(data2) %in% nf)]
              m[k] <- mean(try(.mahal(d1,d2),silent=TRUE),na.rm=TRUE)
            }
            options(warn=0)
          } else {
            for (k in seq_along(np)) {
              ww <- fc[['np']][[np[[k]]]]
              d1 <- data2[ww,which(colnames(data2) %in% nf)]
              m[k] <- abs(mean(d2,na.rm=TRUE) - mean(d1,na.rm=TRUE))
            }
          }
          
          ww <- which.min(m)
          if (length(ww) > 0) dd[w,nFact[i]] <- np[ww]
          else {
            dom.class <- summary(data1[,nFact[i]])
            dom.class <- attr(sort(dom.class, decreasing = TRUE)[1], "names")
            dd[w,nFact[i]] <- dom.class
          }
        } else {
          dom.class <- summary(data1[,nFact[i]])
          dom.class <- attr(sort(dom.class, decreasing = TRUE)[1], "names")
          dd[w,nFact[i]] <- dom.class
        }
        
      }
      
    }
  }
  for (i in seq_along(nFact)) dd[,nFact[i]] <- factor(dd[,nFact[i]])
  .eqFactLevel(data1,dd)
}
#--------

.factorFixW <- function(data1,data2,nFact,nf) {
  # just generates a list to report which class should be assigned to which class!
  if (missing(nFact) || is.null(nFact)) nFact <- colnames(data2)[.where(is.factor,data2)]
  if (missing(nf) || is.null(nf)) {
    nf <- colnames(data2)[!colnames(data2) %in% nFact]
    if (length(nf) == 0) nf <- NULL
  }
  
  dd <- data2
  o <- list()
  for (i in seq_along(nFact)) {
    data1[,nFact[i]] <- factor(data1[,nFact[i]])
    data2[,nFact[i]] <- factor(data2[,nFact[i]])
    dd[,nFact[i]] <- as.character(dd[,nFact[i]])
    
    fc <- .checkFactor(data1[,nFact[i]],data2[,nFact[i]])
    if (!is.null(fc)) {
      p <- names(fc$p)
      np <- names(fc$np)
      for (j in seq_along(p)) {
        w <- fc[['p']][[p[[j]]]]
        if (!is.null(nf)) {
          m <- rep(NA,length(np))
          d2 <- data2[w,which(colnames(data2) %in% nf)]
          if (length(nf) > 1) {
            options(warn=-1)
            for (k in seq_along(np)) {
              ww <- fc[['np']][[np[[k]]]]
              d1 <- data2[ww,which(colnames(data2) %in% nf)]
              m[k] <- mean(try(.mahal(d1,d2),silent=TRUE),na.rm=TRUE)
            }
            if (any(is.na(m))) {
              m <- matrix(NA,nrow=length(nf),ncol=length(np))
              for (k in seq_along(np)) {
                ww <- fc[['np']][[np[[k]]]]
                d1 <- data2[ww,which(colnames(data2) %in% nf)]
                for (nfi in seq_along(nf)) {
                  m[nfi,k] <- abs(mean(d2[[nf[nfi]]],na.rm=TRUE) - mean(d1[[nf[nfi]]],na.rm=TRUE))
                }
              }
              m <- abs(apply(t(apply(m,1,function(x) (x - mean(x)) / sd(x))),2,mean))
            }
            options(warn=0)
          } else {
            for (k in seq_along(np)) {
              ww <- fc[['np']][[np[[k]]]]
              d1 <- data2[ww,which(colnames(data2) %in% nf)]
              m[k] <- abs(mean(d2,na.rm=TRUE) - mean(d1,na.rm=TRUE))
            }
          }
          
          ww <- which.min(m)
          if (length(ww) > 0) {
            o <- c(o,list(c(field=nFact[i],old=unique(dd[w,nFact[i]]),new=np[ww])))
          } else {
            dom.class <- summary(data1[,nFact[i]])
            dom.class <- attr(sort(dom.class, decreasing = TRUE)[1], "names")
            o <- c(o,list(c(field=nFact[i],old=unique(dd[w,nFact[i]]),new=dom.class)))
          }
        } else {
          dom.class <- summary(data1[,nFact[i]])
          dom.class <- attr(sort(dom.class, decreasing = TRUE)[1], "names")
          o <- c(o,list(c(field=nFact[i],old=unique(dd[w,nFact[i]]),new=dom.class)))
        }
      }
    }
  }
  o
}

#--------

.factorFix.bm <- function(data1,data2,nFact,nf) {
  # problematic factors are moved from data2 to data1
  if (missing(nFact) || is.null(nFact)) nFact <- colnames(data2)[.where(is.factor,data2)]
  er <- FALSE
  
  for (i in seq_along(nFact)) {
    data1[,nFact[i]] <- factor(data1[,nFact[i]])
    data2[,nFact[i]] <- factor(data2[,nFact[i]])
    fc <- .checkFactor(data1[,nFact[i]],data2[,nFact[i]])
    if (!is.null(fc)) {
      er <- TRUE
      p <- names(fc$p)
      for (j in seq_along(p)) {
        ww <- fc[['p']][[p[[j]]]]
        ww <- sample(ww,1)
        data1 <- rbind(data1,data2[ww,])
        data2 <- data2[-ww,]
        #dd <- dd[-ww,]
      }
    }
  }
  if (er) {
    for (i in seq_along(nFact)) {
      data1[,nFact[i]] <- factor(data1[,nFact[i]])
      data2[,nFact[i]] <- factor(data2[,nFact[i]])
    }
    return(list(train=data1,test=data2,IDs=fc$p))
  }
}
#--------

.is.windows <- function() {
  s <- Sys.info()
  if (!is.null(s)) s[['sysname']] == 'Windows'
  else FALSE
}
#---------
.getRunID <- function(sm,sp,m) {
  # from sdmModels, extract the modelID and runID (replicates) to be assigned to the model objects generated through fitting
  w1 <- sm@run.info[,2] == sp
  w2 <- sm@run.info[,3] == m
  w1 <- w1 & w2
  list(mID=sm@run.info[w1,1],rID=sm@run.info[w1,5])
}

#-----------
.require <- function(x) {
  # based on simplifying the code of the reqiure function in the base package
  loaded <- paste("package", x, sep = ":") %in% search()
  if (!loaded) {
    value <- tryCatch(library(x,character.only = TRUE, logical.return = TRUE, warn.conflicts = FALSE, quietly = TRUE), error = function(e) e)
    if (inherits(value, "error")) {
      return(FALSE)
    }
    if (!value) return(FALSE)
  } else value <- TRUE
  value
}
#----------
.loadLib <- function(pkgs) {
  options(warn=-1)
  return(unlist(lapply(pkgs,function(x) {
    all(unlist(lapply(x,function(p) {.require(p)})))
  })))
  options(warn=0)
}
#---------
.getRecordID <- function(x,sp,id,train) {
  # x is recordID list
  # it finds the record ID of observation by the rowID used to generate train or dependent test, or in independent test
  # train = FALSE means id belongs to independent test
  if (train) {
    x[[sp]]$train$rID[x[[sp]]$train$rowID %in% id]
  } else {
    x[[sp]]$test$rID[x[[sp]]$test$rowID %in% id]
  }
}

.generateWL <- function(d,s) {
  pkgs <- .sdmMethods$getPackageNames(s@methods)
  .sdm...temp <- NULL; rm(.sdm...temp)
  pos <- 1
  
  w <- sapply(pkgs, function(x) '.temp' %in% x)
  
  if (any(w)) {
    if (".sdm...temp" %in% ls(pattern='^.sdm..',pos=1,all.names = TRUE)) {
      ww <- ls(.sdmMethods$userFunctions)
      rm(list=ww,pos=1)
      rm(.sdm...temp,pos=1)
    }
    #----
    for (i in which(w)) pkgs[[i]] <- pkgs[[i]][which(pkgs[[i]] != '.temp')]
    #----
    w <- ls(.sdmMethods$userFunctions)
    
    if (length(w) > 0) {
      assign('.sdm...temp',c(),envir = as.environment(pos))
      for (ww in w) {
        if (!exists(ww,where=1)) {
          assign(ww,.sdmMethods$userFunctions[[ww]],envir = as.environment(pos))
          .sdm...temp <<- c(.sdm...temp,ww)
        }
      }
    }
  }
  
  ww <- .loadLib(pkgs)
  if (!all(ww)) {
    if (!any(ww)) {
      cat('some methods are removed because they depend on some packages that are not installed on this machine!\n')
      cat('you can use installAll() function to simply install all the packages that may be required by some functions in the sdm package!\n')
      stop(paste('There is no installed packages rquired by the selected methods. Package names:',paste(unlist(pkgs),collapse=', ')))
    } else {
      cat('some methods are removed because they depend on some packages that are not installed on this machine!\n')
      cat('you can use installAll() function to simply install all the packages that may be required by some functions in the sdm package!\n')
      warning(paste('There is no installed packages rquired by the methods:',paste(s@methods[!ww],collapse=', '),'; These methods are excluded! The packages need to be installed for these methods:',paste(unlist(pkgs[!ww]),collapse=', ')))
      s@methods <- s@methods[ww]
    }
  }
  
  #fr <- .getFeaturetype(d,s@sdmFormula)
  w <- new('.workload',ncore=s@ncore,data=d,setting=s,frame=s@featuresFrame)
  hasTest <- 'test' %in% d@groups$training@values[,2]
  nFact <- NULL
  if (!is.null(d@factors) > 0 && any(d@factors %in% s@sdmFormula@vars)) nFact <- d@factors[d@factors %in% s@sdmFormula@vars]
  nf <- .excludeVector(d@features.name,nFact)
  nf <- nf[nf %in% s@sdmFormula@vars]
  nFact <- nFact[nFact %in% s@sdmFormula@vars]
  
  for (sp in s@sdmFormula@species) {
    dt <- as.data.frame(d,sp=sp,grp='train')
    w$recordIDs[[sp]]$train <- data.frame(rID=dt[,1],rowID=1:nrow(dt))
    f <- .getModelFrame(w$frame,dt,response=sp)
    if (!is.null(f$specis_specific)) {
      w$train[[sp]]$sdmDataFrame <- cbind(dt[,sp],f$features,f$specis_specific)
    } else w$train[[sp]]$sdmDataFrame  <- cbind(dt[,sp],f$features)
    colnames(w$train[[sp]]$sdmDataFrame)[1] <- sp
    
    if (hasTest) {
      dt <- as.data.frame(d,sp=sp,grp='test')
      w$recordIDs[[sp]]$test <- data.frame(rID=dt[,1],rowID=1:nrow(dt))
      f <- .getModelFrame(w$frame,dt,response=sp)
      if (!is.null(f$specis_specific)) {
        w$test[[sp]]$sdmDataFrame <- cbind(dt[,sp],f$features,f$specis_specific)
      } else w$test[[sp]]$sdmDataFrame <- cbind(dt[,sp],f$features)
      colnames(w$test[[sp]]$sdmDataFrame)[1] <- sp
      
      if (!is.null(nFact)) {
        
        for (nF in nFact) {
          fc <- .checkFactor(w$train[[sp]]$sdmDataFrame[,nF],w$test[[sp]]$sdmDataFrame[,nF])
          if (!is.null(fc)) {
            p <- names(fc$p)
            for (j in seq_along(p)) {
              ww <- fc[['p']][[p[[j]]]]
              if (length(ww) > 1) ww <- sample(ww,1)
              w$train[[sp]]$sdmDataFrame <- rbind(w$train[[sp]]$sdmDataFrame,w$test[[sp]]$sdmDataFrame[ww,])
              w$test[[sp]]$sdmDataFrame[ww,] <- w$test[[sp]]$sdmDataFrame[-ww,]
              w$data <- .updateGroup(w$data,c(.getRecordID(w$recordIDs,sp=sp,id = ww,train=FALSE),'test','train'))
            }
          }
        }
        
        #f <- .factorFix.bm(w$train[[sp]]$sdmDataFrame,w$test[[sp]]$sdmDataFrame,nFact)
        #if (!is.null(f)) {
        #  w$train[[sp]]$sdmDataFrame <- f$train
        #  w$test[[sp]]$sdmDataFrame <- f$test
        #}
      }
    }
  }
  
  w$funs[['fit']] <- .sdmMethods$getFitFunctions(s@methods)
  w$arguments[['fit']] <- .sdmMethods$getFitArguments(s@methods)
  w$funs[['predict']] <- .sdmMethods$getPredictFunctions(s@methods)
  w$arguments[['predict']] <- .sdmMethods$getPredictArguments(s@methods)
  #w$dataObject.names <- unique(unlist(lapply(s@methods, .sdmMethods$getDataArgumentNames)))
  mo <- s@methods
  names(mo) <- mo
  w$dataObject.names <- lapply(mo, .sdmMethods$getDataArgumentNames)
  #-----------
  
  #reserved.names <- w$getReseved.names()
  for (mo in s@methods) {
    wc <- unlist(lapply(w$arguments$fit[[mo]]$params,function(x) is.character(x)))
    if (any(!wc)) {
      if (!all(unlist(lapply(w$arguments$fit[[mo]]$params[!wc],function(x) is.function(x))))) stop(paste('parameter definition for the model',mo,'in the model container is not correctly defined!'))
      for (n in names(w$arguments$fit[[mo]]$params[!wc])) {
        #if (!all(names(formals(w$arguments$fit[[mo]]$params[[n]])) %in% reserved.names)) stop(paste('the input argument for the function generates the parameter for model',m,'is unknown (not in the reseved objects)'))
        w$params[[sp]][[n]] <- do.call(w$arguments$fit[[mo]]$params[[n]],w$generateParams(names(formals(w$arguments$fit[[mo]]$params[[n]])),sp)) 
        w$arguments$fit[[mo]]$params[[n]] <- n
      }
    }
    
    wc <- unlist(lapply(w$arguments$predict[[mo]]$params,function(x) is.character(x)))
    
    if (any(!wc)) {
      if (!all(unlist(lapply(w$arguments$predict[[mo]]$params[!wc],function(x) is.function(x))))) stop(paste('parameter definition for the model',mo,'in the model container is not correctly defined!'))
      for (n in names(w$arguments$predict[[mo]]$params[!wc])) {
        #if (!all(names(formals(w$arguments$predict[[mo]]$params[[n]])) %in% reserved.names)) stop(paste('the input argument for the function generates the parameter for model',m,'is unknown (not in the reseved objects)'))
        w$params[[sp]][[n]] <- do.call(w$arguments$predict[[mo]]$params[[n]],w$generateParams(names(formals(w$arguments$predict[[mo]]$params[[n]])),sp))
        w$arguments$predict[[mo]]$params[[n]] <- n
      }
    }
  }
  #-----------------
  
  if (!is.null(s@replicate)) {
    f <- .replicateMethods$getFunctions(s@replicate)
    for (sp in names(w$train)) {
      if (d@species[[sp]]@type == 'Presence-Absence') family <- 'binomial'
      else family <- 'xxx'
      
      # sdmDataFrame! Leter should be checked for other types of data!
      for (ff in f) {
        w$replicates[[sp]] <- c(w$replicates[[sp]],ff(x=w$train[[sp]][['sdmDataFrame']][,sp],family=family,stratify=TRUE,test.percent=s@test.percentage,nfolds=s@cv.folds,replicates=s@n.replicates))
      }
    }
  }
  
  # for each replicatios, checks whether the factor level is going to be problematic
  # if so, move the record from test to train
  if (!is.null(nFact)) {
    for (nF in nFact) {
      for (sp in names(w$replicates)) {
        for (i in 1:length(w$replicates[[sp]])) {
          #fc <- .checkFactor(w$train[[sp]]$sdmDataFrame[w$runtasks$runIndex[[sp]][[i]]$train,nF],w$train[[sp]]$sdmDataFrame[w$runtasks$runIndex[[sp]][[i]]$test,nF])
          fc <- .checkFactor(w$train[[sp]]$sdmDataFrame[w$replicates[[sp]][[i]]$train,nF],w$train[[sp]]$sdmDataFrame[w$replicates[[sp]][[i]]$test,nF])
          if (!is.null(fc)) {
            p <- names(fc$p)
            for (j in seq_along(p)) {
              ww <- fc[['p']][[p[[j]]]]
              if (length(ww) > 1) ww <- sample(ww,1)
              w$replicates[[sp]][[i]]$train <- c(w$replicates[[sp]][[i]]$train,w$replicates[[sp]][[i]]$test[ww])
              w$replicates[[sp]][[i]]$test <- w$replicates[[sp]][[i]]$test[-ww]
            }
          }
        }
      }
    }
  }
  #-----
  w
}

#----------------------------------------
if (!isGeneric("sdmSetting")) {
  setGeneric("sdmSetting", function(formula,data,methods,interaction.depth=1,n=1,replication=NULL,
                                    cv.folds=NULL,test.percent=NULL,bg=NULL,bg.n=NULL,var.importance=NULL,response.curve=TRUE,
                                    var.selection=FALSE,ncore=1L,...)
    standardGeneric("sdmSetting"))
}

setMethod('sdmSetting', signature(formula='ANY','sdmdata','character'), 
          function(formula,data,methods,interaction.depth=1,n=1,replication=NULL,
                   cv.folds=NULL,test.percent=NULL,bg=NULL,bg.n=NULL,var.importance=NULL,response.curve=TRUE,
                   var.selection=FALSE,ncore=1L,...) {
            
            if (!.sdmOptions$getOption('sdmLoaded')) .addMethods()
            
            dot <- list(...)
            sobj <- NULL
            if (length(dot) > 0) {
              ndot <- names(dot)
              if ('' %in% ndot) {
                for (i in seq_along(which(ndot == ''))) {
                  if (inherits(dot[[i]],'.sdmCorSetting')) {
                    sobj <- dot[[i]]
                    break
                  }
                }
                dot <- dot[-which(ndot == '')]
                ndot <- names(dot)
              }
              
              a <- c('interaction.depth','replication','cv.folds','test.percent','bg','bg.n','var.importance','response.curve','var.selection','ncore')
              ndot <- .pmatch(ndot,a)
              w <- !is.na(ndot)
              if (length(w) > 0) {
                dot <- dot[w]
                ndot <- ndot[w]
                names(dot) <- ndot
              }
              
              
              if ('setting' %in% names(dot) && inherits(dot[['setting']],'.sdmCorSetting')) {
                sobj <- dot[['setting']]
                dot <- dot[-which(ndot == 'setting')]
                ndot <- names(dot)
              }
              
              if (length(dot) > 0) {
                if (length(ndot) > 0) {
                  for (nd in ndot) {
                    if (nd == 'interaction.depth' && interaction.depth == 1) interaction.depth <- dot[[nd]]
                    else if (nd == 'ncore' && ncore == 1L) ncore <- dot[[nd]]
                    else if (nd == 'replication' && is.null(replication)) replication <- dot[[nd]]
                    else if (nd == 'cv.folds' && is.null(cv.folds)) cv.folds <- dot[[nd]]
                    else if (nd == 'test.percent' && is.null(test.percent)) test.percent <- dot[[nd]]
                    else if (nd == 'bg' && is.null(bg)) bg <- dot[[nd]]
                    else if (nd == 'bg.n' && is.null(bg.n)) bg.n <- dot[[nd]]
                    else if (nd == 'var.importance' && is.null(var.importance)) var.importance <- dot[[nd]]
                    else if (nd == 'response.curve' && response.curve && is.logical(dot[[nd]])) response.curve <- dot[[nd]]
                    else if (nd == 'var.selection' && !var.selection && is.logical(dot[[nd]])) var.selection <- dot[[nd]]
                  }
                }
              }
            }
            #--------
            
            m <- .methodFix(methods)
            if (any(is.na(m))) stop(paste('methods',paste(methods[is.na(m)],collapse=', '),'do not exist!'))
            m <- unique(m)
            
            s <- new('.sdmCorSetting',methods=m)
            s@distribution <- .getSpeciesDistribution(data)
            
            if (missing(formula)) {
              if (!is.null(sobj)) {
                if (all(sobj@sdmFormula@vars %in% data@features.name)) s@sdmFormula <- sobj@sdmFormula
                else s@sdmFormula <- data@sdmFormula
              } else s@sdmFormula <- data@sdmFormula
              
            } else if (inherits(formula,'sdmFormula')) s@sdmFormula <- formula
            else if (inherits(formula,'formula')) {  
              s@sdmFormula <- .exFormula(formula,as.data.frame(data)[,-1])
            } else if (inherits(formula,'.sdmCorSetting')) {
              sobj <- formula
              if (all(sobj@sdmFormula@vars %in% data@features.name)) s@sdmFormula <- sobj@sdmFormula
              else s@sdmFormula <- data@sdmFormula
            } else {
              if (!is.null(sobj)) {
                if (all(sobj@sdmFormula@vars %in% data@features.name)) s@sdmFormula <- sobj@sdmFormula
                else s@sdmFormula <- data@sdmFormula
              } else s@sdmFormula <- data@sdmFormula
            }
            
            s@featuresFrame <- .getFeaturetype(data,s@sdmFormula)  
              
            if (!is.null(test.percent)) s@test.percentage <- test.percent
            else {
              if (!is.null(sobj)) {
                if (!is.null(sobj@test.percent)) s@test.percentage <- sobj@test.percent
              }
            }
            
            s@interaction.depth <- interaction.depth
            if (interaction.depth ==1 && !is.null(sobj) && !is.null(sobj@interaction.depth)) s@interaction.depth <- sobj@interaction.depth
            
            s@ncore <- ncore
            if (ncore == 1L && !is.null(sobj) && length(sobj@ncore) == 1) s@ncore <- sobj@ncore
            
            if (!is.null(replication)) {
              nx <- .replicate.methodFix(replication)
              if (any(is.na(nx))) warning(paste(paste(replication[is.na(nx)],collapse=', '),'methods in replication are not found [They are ignored!]'))
              replication <- nx[!is.na(nx)]
              s@replicate <- replication
            } else {
              if (!is.null(sobj)) {
                if (!is.null(sobj@replicate)) s@replicate <- sobj@replicate
              }
              if (is.null(s@replicate) && !is.null(s@test.percentage)) {
                s@replicate <- "subsampling"
              }
            }
            
            s@n.replicates <- n
            if (!is.null(sobj) && !is.null(sobj@n.replicates)) s@n.replicates <- sobj@n.replicates
            
            if ("subsampling" %in% s@replicate) {
              if (is.null(s@test.percentage)) s@test.percentage <- 30
            }
            
            if (!is.null(cv.folds)) s@cv.folds <- cv.folds
            else {
              if (!is.null(sobj) && !is.null(sobj@cv.folds)) s@cv.folds <- sobj@cv.folds
              if (is.null(s@cv.folds) && "cross_validation" %in% s@replicate) s@cv.folds <- 5
            }
            
            if (!is.null(s@cv.folds) && !"cross_validation" %in% s@replicate) {
              s@replicate <- c("cross_validation",s@replicate)
            }
            
            if (!is.null(bg)) s@pseudo.absence.methods <- bg
            else {
              if (!is.null(sobj)) {
                if (!is.null(sobj@pseudo.absence.methods)) s@pseudo.absence.methods <- sobj@pseudo.absence.methods
              }
            }
            if (!is.null(bg.n)) s@n.pseudo.absence <- bg.n
            else {
              if (!is.null(sobj)) {
                if (!is.null(sobj@n.pseudo.absence)) s@n.pseudo.absence <- sobj@n.pseudo.absence
              }
              if (is.null(s@n.pseudo.absence) && !is.null(s@pseudo.absence.methods)) {
                s@n.pseudo.absence <- 1000
              }
            }
            if (!is.null(var.importance)) s@varImportance.methods <- var.importance
            else {
              if (!is.null(sobj)) {
                if (!is.null(sobj@varImportance.methods)) s@varImportance.methods <- sobj@varImportance.methods
              }
            }
            if (response.curve) s@response.curve <- TRUE
            else {
              if (!is.null(sobj)) {
                if (!is.null(sobj@response.curve) && sobj@response.curve) s@response.curve <- sobj@response.curve
              } else s@response.curve <- FALSE
            }
            
            if (var.selection) s@var.selection <- TRUE
            else {
              if (!is.null(sobj)) {
                if (!is.null(sobj@var.selection) && sobj@var.selection) s@var.selection <- sobj@var.selection
              } else s@var.selection <- FALSE
            }
            
            if (!is.null(interaction.depth)) s@interaction.depth <- interaction.depth
            else {
              if (!is.null(sobj)) {
                if (!is.null(sobj@interaction.depth)) s@interaction.depth <- sobj@interaction.depth
              }
            }
            s
          }
)
#----------------
if (!isGeneric("sdm")) {
  setGeneric("sdm", function(formula,data,methods,...)
    standardGeneric("sdm"))
}

setMethod('sdm', signature(formula='formula',data='sdmdata',methods='character'), 
          function(formula,data,methods,...) {
            a <- c('interaction.depth','n','replication','cv.folds','test.percent','bg','bg.n','var.importance','response.curve','var.selection','setting','ncore')
            .sdm...temp <- NULL; rm(.sdm...temp)
            dot <- list(...)
            ndot <- names(dot)
            if (length(ndot) > 0) {
              ndot <- .pmatch(ndot,a)
              w <- !is.na(ndot)
              ndot <- ndot[w]
              dot <- dot[w]
              names(dot) <- ndot
            }
            
            dot$data <- data
            dot$formula <- formula
            dot$methods <- methods
            
            if (!.sdmOptions$getOption('sdmLoaded')) .addMethods()
            
            s <- do.call('sdmSetting',dot)
            w <- .generateWL(data,s)
            w <- w$fit()
            #if (".sdmMethods$userFunctions" %in% search()) detach('.sdmMethods$userFunctions')
            if (".sdm...temp" %in% ls(pattern='^.sdm..',pos=1,all.names = TRUE)) {
              ww <- ls(.sdmMethods$userFunctions)
              rm(list=ww,pos=1)
              rm(.sdm...temp,pos=1)
            }
            w
          }
)

#-------------
# .getModel.info <- function(x,w,...) {
#   if (missing(w) || is.null(w)) {
#     a <- c('species','method','replication','run')
#     w1 <- w2 <- w3 <- w4 <- TRUE
#     dot <- list(...)
#     if (length(dot) > 0) {
#       ndot <- names(dot)
#       ndot <- .pmatch(ndot,a)
#       w <- !is.na(ndot)
#       ndot <- ndot[w]
#       dot <- dot[w]
#       names(dot) <- ndot
#       for (nd in ndot) {
#         if (nd == 'species' && !is.null(dot[[nd]])) {
#           dot[[nd]] <- .pmatch(dot[[nd]],unique(as.character(x@run.info[,2])))
#           w1 <- x@run.info[,2] %in% dot[[nd]]
#         } else if (nd == 'method' && !is.null(dot[[nd]])) {
#           dot[[nd]] <- .methodFix(dot[[nd]])
#           w2 <- x@run.info[,3] %in% dot[[nd]]
#         }
#         else if (nd == 'replication' && !is.null(dot[[nd]])) {
#           if (length(x@replicates) != 0) {
#             dot[[nd]] <- .replicate.methodFix(dot[[nd]])
#             w3 <- x@run.info[,4] %in% dot[[nd]]
#           }
#         } else if (nd == 'run' && !is.null(dot[[nd]])) {
#           if (!is.null(dot[[nd]])) {
#             if (length(x@replicates) != 0) {
#               r <- unlist(lapply(x@replicates[[1]],function(x) x$method))
#               ru <- unique(r)
#               names(ru) <- ru
#               rID <- lapply(ru,function(x) which(r == x))
#               w4 <- c()
#               for (i in 1:length(rID)) {
#                 w4 <- c(w4,rID[[i]][c(1:length(rID[[i]])) %in% dot[[nd]]])
#               }
#               w4 <- x@run.info[,5] %in% w4
#             }
#           }
#         }
#       }
#       x@run.info[w1 & w2 & w3 & w4,]
#     } else x@run.info
#   } else x@run.info[x@run.info[,1] %in% w,]
# }
#--------
.getModel.info <- function(x,w=NULL,species=NULL,method=NULL,replication=NULL,run=NULL) {
  if (missing(w) || is.null(w)) {
    
    w1 <- w2 <- w3 <- w4 <- TRUE
    if (!is.null(species)) {
      species <- .pmatch(species,unique(as.character(x@run.info[,2])))
      w1 <- x@run.info[,2] %in% species
    }
    
    if (!is.null(method)) {
      method <- .methodFix(method)
      w2 <- x@run.info[,3] %in% method
    }
    
    
    if (!is.null(replication)) {
      if (length(x@replicates) != 0) {
        replication <- .replicate.methodFix(replication)
        w3 <- x@run.info[,4] %in% replication
      }
    }
    
    if (!is.null(run)) {
      if (!is.null(run)) {
        if (length(x@replicates) != 0) {
          r <- unlist(lapply(x@replicates[[1]],function(x) x$method))
          ru <- unique(r)
          names(ru) <- ru
          rID <- lapply(ru,function(x) which(r == x))
          w4 <- c()
          for (i in 1:length(rID)) {
            w4 <- c(w4,rID[[i]][c(1:length(rID[[i]])) %in% run])
          }
          
          w4 <- x@run.info[,5] %in% w4
        }
      }
    }
    x@run.info[w1 & w2 & w3 & w4,]
  } else x@run.info[x@run.info[,1] %in% w,]
}

#--------


.getModel.info2 <- function(x,w=NULL,species=NULL,method=NULL,replication=NULL,run=NULL,wtest=NULL) {
  # comparing to .getModel.info: In this, only one species is allowed!
  # x: sdmModels
  if (!is.null(w)) {
    mi <- .getModel.info(x,w)
  } else {
    mi <- .getModel.info(x)
    u <- as.character(unique(mi[,2]))
    m <- as.character(unique(mi[,3]))
    r <- unique(mi[,4])
    if (!is.null(species)) {
      if (length(species) > 1) {
        species <- species[1]
        warning('only the first species is considered!')
      }
      if (is.numeric(species)) {
        if (length(u) <= species) species <- u[species]
        else stop('The specified species is not recognised!')
      } else {
        species <- .pmatch(species,u)
        if (is.na(species)) stop('The specified species is not recognised!')
      }
    } else {
      if (length(u) > 1) stop('This object contains models for more than one species; in species argument spcify the name of species!')
      else species <- u
    }
    
    if (!is.null(method)) {
      method <- .sdmMethods$fixNames(method)
      wm <- method %in% m
      if (any(!wm)) {
        if (all(!wm)) stop('the specified methods do not exist in the object!')
        warning(paste('Methods',paste(method[!wm],collapse=', '),'do not exsit in the object, and are excluded!'))
        method <- method[wm]
      }
    } else method <- m
    
    if (!is.null(replication)) replication <- .replicate.methodFix(replication)
    else replication <- r
    
    mi <- .getModel.info(x,species=species,method=method,replication=replication,run=run)
  }
  
  mi
}
#---------