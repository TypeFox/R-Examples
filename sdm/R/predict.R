# Author: Babak Naimi, naimi.b@gmail.com
# Date :  Aptil 2016
# Version 2.3
# Licence GPL v3
# 
# .raster2data.table <- function(r) {
#   if (inherits(r,'RasterBrick'))  {
#     o <- data.table(r@data@values)
#     o$cellnr <- 1:ncell(r)
#   } else if (inherits(r,'RasterLayer')) {
#     o <- data.table(r@data@values)
#     colnames(o) <- names(r)
#     o$cellnr <- 1:ncell(r)
#   } else {
#     o <- data.table(as.data.frame(r))
#     o$cellnr <- 1:ncell(r)
#   }
#   o
# }
#---------

.getlevels <- function(x) {
  o <- NULL
  if (inherits(x,'sdmdata')) {
    if (!is.null(x@factors)) {
      o <- x@factors
      names(o) <- o
      o <- lapply(o,function(i) levels(x@features[[i]]))
    }
  } else if (inherits(x,'sdmModels')) {
    
    if (!is.null(x@data@factors)) {
      o <- x@data@factors
      names(o) <- o
      o <- lapply(o,function(i) levels(x@data@features[[i]]))
    }
  } else if (inherits(x,'data.frame')) {
    f <- .where(is.factor,x)
    if (any(f)) {
      f <- names(f)[f]
      o <- vector('list',length(f))
      names(o) <- f
      for (i in seq_along(o)) o[[i]] <- levels(x[[f[i]]])
    }
  }
  o
}
#-------------

.getTotal.object.size <- function() {
  # only the size of objects in R_GlobalEnv
  paste(as.character(round(sum(unlist(lapply(ls(all.names = TRUE,envir=parent.frame()),function(x) object.size(get(x)))))/1024/1024/1024,4)),'Gb')
}

#-----------
.raster2df <- function(x,level) {
  d <- data.frame(cellnr=1:ncell(x),as.data.frame(x))
  bb <- rep(TRUE,nrow(d))
  for (i in 2:ncol(d)) bb <- bb & !is.na(d[,i])
  if (length(which(bb)) == 0) stop('raster object has no data...!')
  d <- d[bb,]
  rm(bb)
  
  if (!missing(level) && !is.null(level)) {
    n <- names(level)
    for (i in seq_along(n)) {
      l <- level[[i]]
      u <- sort(unique(d[,n[i]]))
      if (any(u %in% c(1:length(l)))) {
        u <- u[u %in% c(1:length(l))]
        l <- l[u]
        d[,n[i]] <- factor(d[,n[i]])
        levels(d[,n[i]]) <- l
      } else stop('the grid values in categorical rasters does not match with the factor levels in the model')
    }
  }
  d
}
#----------------

.generateName <- function(x) {
  paste(c(x,'_',sample(c(letters,1:9),6,replace=T)),collapse='')
}
#-----
.domClass <- function(v) {
  # decreasing sort of dominant classes in a character vector
  v <- as.character(v)
  u <- unique(v)
  names(u) <- u
  u <- unlist(lapply(u,function(x) length(which(v == x))))
  names(u)[order(u,decreasing = TRUE)]
}
#---------


if (!isGeneric("predict")) {
  setGeneric("predict", function(object, ...)
    standardGeneric("predict"))
}	
# 
# .predict=function(obj,pred,pred.par,dt=dt) {
#   pred.par[[1]] <- obj
#   pred.par[[2]] <- dt
#   m <- try(pred(pred.par),silent=TRUE)
#   options(warn=0)
#   m
# }




.generateWLP <- function(x,newdata,w=NULL,species=NULL,method=NULL,replication=NULL,run=NULL,ncore=NULL) {
  
  mi <- .getModel.info(x,w=w,species=species,method=method,replication=replication,run=run)
  
  s <- mi$success
  
  if (!all(s)) {
    if (!any(s)) stop('There is no model objects that were successfully fitted!')
    if (length(which(!s)) == 1) {
      warning(paste('1 model from the total',length(s),'models was NOT successfully fitted, so it is excluded in prediction!'))
    } else {
      warning(paste(length(which(!s)),'models from the total',length(s),'models were NOT successfully fitted, so they are excluded in prediction!'))
    }
    
    mi <- mi[s,]
    
  }
  mi <- mi[,1:5]
  #-----------
  nr <- nrow(mi)
  if (nr == 0) stop('the specified models do not exist!')
  species <- as.character(unique(mi[,2]))
  
  m <- unique(as.character(mi[,3]))
  #----
  w <- new('.workloadP',runTasks=mi)
  
  if (is.null(ncore)) w$ncore <- 1L
  else {
    w$ncore <- parallel::detectCores()
    if (ncore < w$ncore) w$ncore <- ncore
    # temporary until the HPC for windows is implemented:
    if (.is.windows()) w$ncore <- 1L
  }
  
  w$newdata$raster <- NULL
  nf <- nFact <- NULL
  
  if (inherits(newdata,'data.frame')) {
    n <- colnames(newdata)
    if (!all(x@setting@featuresFrame@vars %in% n)) stop('the data does not contain some or all of the variables that the model needs...')
    w$newdata$data.frame <- newdata
    w$modelFrame <- .getModelFrame(x@setting@featuresFrame,w$newdata$data.frame,response=species)
  } else if (inherits(newdata,'Raster')) {
    n <- names(newdata)
    if (!all(x@setting@featuresFrame@vars %in% n)) stop('the data does not contain some or all of the variables that the model needs...')
    w$newdata$raster <- newdata
    w$newdata$data.frame <- .raster2df(newdata,.getlevels(x))
    
    w$modelFrame <- .getModelFrame(x@setting@featuresFrame,w$newdata$data.frame,response=species)
    
    #b <- brick(raster(newdata))
    #mr <- rep(NA,ncell(b))
    
  } else stop('newdata should be a Raster* object or a data.frame...!')
  
  
  w$funs <- .sdmMethods$getPredictFunctions(m)
  w$arguments <- .sdmMethods$getPredictArguments(m)
  w$dataObject.names <- unique(unlist(lapply(x@setting@methods, .sdmMethods$getDataArgumentNames)))
  for (mo in m) {
    wc <- unlist(lapply(w$arguments[[mo]]$params,function(x) is.character(x)))
    
    if (any(!wc)) {
      if (!all(unlist(lapply(w$arguments[[mo]]$params[!wc],function(x) is.function(x))))) stop(paste('parameter definition for the model',m,'in the model container is not correctly defined!'))
      for (n in names(w$arguments[[mo]]$params[!wc])) {
        #if (!all(names(formals(w$arguments$predict[[mo]]$params[[n]])) %in% reserved.names)) stop(paste('the input argument for the function generates the parameter for model',m,'is unknown (not in the reseved objects)'))
        w$params[[n]] <- do.call(w$arguments[[mo]]$params[[n]],w$generateParams(names(formals(w$arguments[[mo]]$params[[n]])),sp))
        w$arguments$predict[[mo]]$params[[n]] <- n
      }
    }
  }
  
  #######
  nf <- .getFeatureNamesTypes(x@setting@featuresFrame)
  if ('factor' %in% nf) {
    nFact <- names(nf)[nf == 'factor']
    nf <- .excludeVector(names(nf),nFact)
    id <- mi[,1]
    dd <- as.data.frame(x@data)
    ddf <- .getModelFrame(x@setting@featuresFrame,dd,response=species)
    dd <- dd$rID
    for (sp in species) {
      mj <- mi[which(mi[,2] == sp),]
      r <- mj[,5]
      
      if (!is.null(ddf$specis_specific)) {
        d1 <- cbind(ddf$features,ddf$specis_specific[[sp]])
        d2 <- cbind(w$modelFrame$features,w$modelFrame$specis_specific[[sp]])
        nf1 <- colnames(w$modelFrame$features)
        nf2 <- colnames(w$modelFrame$specis_specific[[sp]])
      } else {
        d1 <- ddf$features
        d2 <- w$modelFrame$features
        nf1 <- colnames(w$modelFrame$features)
        nf2 <- NULL
      }
      
      if (!all(is.na(r))) {
        r <- unique(r)
        o <- data.frame(matrix(ncol=3,nrow=0))
        for (j in r) {
          ddd <- .factorFixW(d1[dd %in% .getRecordID(x@recordIDs,x@replicates[[sp]][[j]]$train,sp=sp,train=TRUE),],d2,nf = nf,nFact=nFact)
          
          if (length(ddd) > 0) {
            for (i in seq_along(ddd)) {
              o <- rbind(o,data.frame(f=ddd[[i]][[1]],old=ddd[[i]][[2]],new=ddd[[i]][[3]]))
            }
          }
        }
        
        if (nrow(o) > 0) {
          for (n in as.character(unique(o[,1]))) {
            wn <- which(o$f == n)
            oc <- o[wn,]
            u <- as.character(unique(oc[,2]))
            un <- as.character(unique(o[wn,3]))
            for (uu in u) {
              wc <- which(oc$old == uu)
              nc <- .domClass(as.character(oc$new[wc]))
              if (nc[1] %in% u && length(nc) > 1) nc <- nc[2]
              else nc <- nc[1]
              ww <- which(w$modelFrame$features[,n] == uu)
              w$modelFrame$features[ww,n] <- nc
            }
            w$modelFrame$features[,n] <- factor(w$modelFrame$features[,n])
          }
        }
      }
    }
  }
  #--------
  w$obj <- x@models
  w$runTasks$species <- as.character(w$runTasks$species)
  w$runTasks$method <- as.character(w$runTasks$method)
  sp <- as.character(mi[,2])
  nw <- unique(sp)
  w$runTasks$speciesID <- unlist(lapply(sp,function(x) {which(nw == x)}))
  
  m <- as.character(mi[,3])
  nw <- names(w$funs)
  w$runTasks$methodID <- unlist(lapply(m,function(x) {which(nw == x)}))
  w$runTasks$mIDChar <- as.character(mi[,1])
  w
}


setMethod('predict', signature(object='sdmModels'), 
          function(object, newdata, filename="",w=NULL,species=NULL,method=NULL,replication=NULL,run=NULL,mean=FALSE,control=NULL,overwrite=FALSE,nc=1L,obj.size=1,err=FALSE,...) {
            if (missing(newdata)) stop('mewdata is missing...')
            #---
            if (!.sdmOptions$getOption('sdmLoaded')) .addMethods()
            #---
            if (missing(method) || is.null(method)) method <- object@setting@methods
            pkgs <- .sdmMethods$getPackageNames(method)
            
            
            .sdm...temp <- NULL; rm(.sdm...temp)
            pos <- 1
            
            tmp <- sapply(pkgs, function(x) '.temp' %in% x)
            
            if (any(tmp)) {
              if (".sdm...temp" %in% ls(pattern='^.sdm..',pos=1,all.names = TRUE)) {
                ww <- ls(.sdmMethods$userFunctions)
                rm(list=ww,pos=1)
                rm(.sdm...temp,pos=1)
              }
              #----
              for (i in which(tmp)) pkgs[[i]] <- pkgs[[i]][which(pkgs[[i]] != '.temp')]
              #----
              tmp<- ls(.sdmMethods$userFunctions)
              
              if (length(tmp) > 0) {
                assign('.sdm...temp',c(),envir = as.environment(pos))
                for (ww in tmp) {
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
                stop(paste('There is no installed packages rquired by the selected methods. Package names:',paste(unlist(pkgs),collapse=', ')))
              } else {
                warning(paste('There is no installed packages rquired by the methods:',paste(object@setting@methods[!ww],collapse=', '),'; These methods are excluded! The packages need to be installed for these methods:',paste(unlist(pkgs[!ww]),collapse=', ')))
                method <- method[ww] 
              }
            }
            
            
            b <- NULL
            w <- .generateWLP(x = object,newdata=newdata,w=w,species=species,method=method,replication=replication,run=run,ncore=nc)
            #w <- sdm:::.generateWLP(x = object,newdata=newdata,w=NULL,species=NULL,method=NULL,replication=NULL,run=NULL,ncore=1)
            
            if (!is.null(w$newdata$raster)) {
              if (filename == '') filename <- .generateName('sdm_prediction')
              if (extension(filename) == '') filename <- paste(filename,'.grd',sep='')
              b <- brick(raster(w$newdata$raster))
              mr <- rep(NA,ncell(b))
            }
            
            mid <- w$runTasks$modelID
            #----
            nr <- nrow(w$runTasks)
            
            obj.size <- obj.size * 1073741824 # Gb to byte
            #--------------
            success <- rep(TRUE,nr)
            rnames <- fullnames <- c()
            errLog <- list()
            
            options(warn=-1)
            
            if (inherits(b,'Raster')) {
              if (nr > 1) {
                b <- brick(b,nl=nr)
                writeRaster(b,filename=filename,overwrite=overwrite)
                b <- brick(filename)
              } else {
                b <- raster(b)
                writeRaster(b,filename=filename,overwrite=overwrite)
                b <- raster(filename)
              }
              
            } else {
              mtx <- matrix(NA,nrow=nrow(w$newdata$data.frame),ncol=0)
            }
            #-------
            memreq <- (object.size(w$newdata$data.frame)[[1]] / ncol(w$newdata$data.frame))*nr
            co <- 1
            
            if (memreq <= obj.size) {
              o <- parallel::mclapply(w$runTasks$modelID,function(i,...) w$predictID(i),mc.cores = w$ncore,w=w)
              if (inherits(b,'Raster')) {
                if (length(o) > 1) {
                  for (j in seq_along(o)) {
                    if (!inherits(o[[j]],'try-error')) {
                      mr[w$newdata$data.frame$cellnr] <- o[[j]]
                      b <- update(b,mr,cell=1,band=co)
                      co <- co+1
                      rnames <- c(rnames,paste('id_',w$runTasks$modelID[j],'-sp_',w$runTasks$speciesID[j],'-m_',w$runTasks$method[j],if (!is.na(w$runTasks$replication[j])) paste('-re_',paste(strsplit(w$runTasks$replication[j],'')[[1]][1:4],collapse=''),sep=''),sep=''))
                      fullnames <- c(fullnames,paste('id_',w$runTasks$modelID[j],'-species_',w$runTasks$species[j],'-method_',w$runTasks$method[j],if (!is.na(w$runTasks$replication[j])) paste('-replication_',w$runTasks$replication[j],sep=''),sep=''))
                    } else {
                      success[j] <- FALSE
                      errLog <- c(errLog,o[[j]])
                    }
                  }
                } else {
                  if (!inherits(o[[1]],'try-error')) {
                    mr[w$newdata$data.frame$cellnr] <- o[[1]]
                    b <- update(b,mr,cell=1)
                    rnames <- paste('id_',w$runTasks$modelID[1],'-sp_',w$runTasks$speciesID[1],'-m_',w$runTasks$method[1],sep='')
                    fullnames <- paste('id_',w$runTasks$modelID[1],'-species_',w$runTasks$species[1],'-method_',w$runTasks$method[1],if (!is.na(w$runTasks$replication[1])) paste('-replication_',w$runTasks$replication[1],sep=''),sep='')
                  } else {
                    success[1] <- FALSE
                    errLog <- c(errLog,o[[1]])
                  }
                }
              } else {
                for (j in seq_along(o)) {
                  if (!inherits(o[[j]],'try-error')) {
                    mtx <- cbind(mtx,o[[j]])
                    co <- co+1
                    rnames <- c(rnames,paste('id_',w$runTasks$modelID[j],'-sp_',w$runTasks$speciesID[j],'-m_',w$runTasks$method[j],if (!is.na(w$runTasks$replication[j])) paste('-re_',paste(strsplit(w$runTasks$replication[j],'')[[1]][1:4],collapse=''),sep=''),sep=''))
                  } else {
                    success[j] <- FALSE
                    errLog <- c(errLog,o[[j]])
                  }
                } 
              }
            } else {
              memdiv <- ceiling(obj.size / (memreq /  nr))
              ii <- ceiling(nr/memdiv)
              for (i in 1:ii) {
                id <- (i-1) * memdiv + c(1:memdiv)
                id <- id[id %in% 1:nr]
                o <- parallel::mclapply(w$runTasks$modelID[id],function(i,...) w$predictID(i),mc.cores = w$ncore,w=w)
                
                if (inherits(b,'RasterBrick')) {
                  for (j in seq_along(o)) {
                    if (!inherits(o[[j]],'try-error')) {
                      mr[w$newdata$data.frame$cellnr] <- o[[j]]
                      b <- update(b,mr,cell=1,band=co)
                      co <- co+1
                      rnames <- c(rnames,paste('id_',w$runTasks$modelID[id[j]],'-sp_',w$runTasks$speciesID[id[j]],'-m_',w$runTasks$method[id[j]],if (!is.na(w$runTasks$replication[id[j]])) paste('-re_',paste(strsplit(w$runTasks$replication[id[j]],'')[[1]][1:4],collapse=''),sep=''),sep=''))
                      fullnames <- c(fullnames,paste('id_',w$runTasks$modelID[id[j]],'-species_',w$runTasks$species[id[j]],'-method_',w$runTasks$method[id[j]],if (!is.na(w$runTasks$replication[id[j]])) paste('-replication_',w$runTasks$replication[id[j]],sep=''),sep=''))
                    } else {
                      success[id[j]] <- FALSE
                      errLog <- c(errLog,o[[id[j]]])
                    }
                  }
                } else if (inherits(b,'RasterLayer')) {
                  if (!inherits(o[[1]],'try-error')) {
                    mr[w$newdata$data.frame$cellnr] <- o[[1]]
                    b <- update(b,mr,cell=1)
                    rnames <- paste('id_',w$runTasks$modelID[1],'-sp_',w$runTasks$speciesID[1],'-m_',w$runTasks$method[1],sep='')
                    fullnames <- paste('id_',w$runTasks$modelID[1],'-species_',w$runTasks$species[1],'-method_',w$runTasks$method[1],if (!is.na(w$runTasks$replication[1])) paste('-replication_',w$runTasks$replication[1],sep=''),sep='')
                  } else {
                    success[1] <- FALSE
                    errLog <- c(errLog,o[[1]])
                  }
                } else {
                  for (j in seq_along(o)) {
                    if (!inherits(o[[j]],'try-error')) {
                      mtx <- cbind(mtx,o[[j]])
                      co <- co+1
                      rnames <- c(rnames,paste('id_',w$runTasks$modelID[id[j]],'-sp_',w$runTasks$speciesID[id[j]],'-m_',w$runTasks$method[id[j]],if (!is.na(w$runTasks$replication[id[j]])) paste('-re_',paste(strsplit(w$runTasks$replication[id[j]],'')[[1]][1:4],collapse=''),sep=''),sep=''))
                    } else {
                      success[id[j]] <- FALSE
                      errLog <- c(errLog,o[[id[j]]])
                    }
                  }
                }
              }
            }
            #------------
            options(warn=0)
            if (!any(success)) {
              if (inherits(b,'Raster')) {
                rm(b)
                unlink(filename)
                stop('Error in prediction....!')
              } else {
                rm(mtx)
                stop('Error in prediction....!')
              }
            }
            
            if (!all(success)) {
              warning(paste(length(which(!success)),' models (out of ',length(success),') failed in the prediction!',sep=''))
              b <- b[[1:length(which(success))]]
              mid <- mid[success]
              w$runTasks <-  w$runTasks[success,]
            }
            
            #----------
            
            
            if (nr > 1 & mean & !is.na(w$runTasks$replication[1])) {
              species <- unique(as.character(w$runTasks$species))
              m <- unique(as.character(w$runTasks$method))
              
              if (inherits(b,'Raster')) {
                rnames <- c()
                newfullnames <- c()
                newr <- raster(b)
                
                for (sp in species) {
                  m1 <- w$runTasks[which(w$runTasks$species == sp),]
                  for (mo in m) {
                    m2 <- m1[which(m1$method == mo),]
                    re <- unique(as.character(m2$replication))
                    for (r in re) {
                      m3 <- m2[which(m2$replication == r),]
                      w <- which(mid %in% m3$modelID)
                      
                      if (length(w) > 1) {
                        temp <- calc(b[[w]],function(x) mean(x,na.rm=TRUE))
                        newr <- addLayer(newr,temp)
                      } else newr <- addLayer(newr,b[[w]])
                      rnames <- c(rnames,paste('sp_',m3$speciesID[1],'-m_',mo,paste('-re_',paste(strsplit(r,'')[[1]][1:4],collapse=''),sep=''),sep=''))
                      newfullnames <- c(newfullnames,paste('species_',sp,'-method_',mo,paste('-replication (Mean)_',m3$replication[1],sep=''),sep=''))
                    }
                  }
                }
                rm(b)
                #unlink(filename)
                b <- brick(newr,filename=filename,values=TRUE,overwrite=TRUE) 
                fullnames <- newfullnames
                rm(newr,newfullnames)
              } else {
                newmtx <- matrix(NA,nrow=nrow(w$newdata$data.frame),ncol=0)
                rnames <- c()
                fullnames <- c()
                for (sp in species) {
                  m1 <- w$runTasks[which(w$runTasks$species == sp),]
                  for (mo in m) {
                    m2 <- m1[which(m1$method == mo),]
                    re <- unique(as.character(m2$replication))
                    for (r in re) {
                      m3 <- m2[which(m2$replication == r),]
                      w <- which(mid %in% m3$modelID)
                      
                      if (length(w) > 1) {
                        temp <- apply(mtx[,w],1,function(x) mean(x,na.rm=TRUE))
                        newmtx <- cbind(newmtx,temp)
                      } else newmtx <- cbind(newmtx,mtx[,w])
                      rnames <- c(rnames,paste('sp_',m3$speciesID[1],'-m_',mo,paste('-re_',paste(strsplit(r,'')[[1]][1:4],collapse=''),sep=''),sep=''))
                      #fullnames <- c(fullnames,paste('species_',sp,'-method_',mo,paste('-replication (Mean)_',m3$replication[1],sep=''),sep=''))
                    }
                  }
                }
                mtx <- newmtx
              }
            }
            #------------
            if (err && length(errLog) > 0) {
              for (i in seq_along(errLog)) cat(errLog[[i]],'\n')
            }
            
            if (".sdm...temp" %in% ls(pattern='^.sdm..',pos=1,all.names = TRUE)) {
              w <- ls(.sdmMethods$userFunctions)
              rm(list=w,pos=1)
              rm(.sdm...temp,pos=1)
            }
            #if (".sdmMethods$userFunctions" %in% search()) detach('.sdmMethods$userFunctions')
            
            if (!inherits(b,'Raster')) {
              colnames(mtx) <- rnames
              return(mtx)
            } else {
              names(b) <- rnames
              b <- setZ(b,fullnames,name='fullname')
              return(b)
            }
          }
)

