# Author: Babak Naimi, naimi.b@gmail.com
# Date :  March 2016
# Version 2.2
# Licence GPL v3

#------
.newgroup <- function(name,values,index) {
  g <- new('.group',name=name)
  if (!is.null(values)) {
    if ((length(values) != length(index)) || !is.list(index)) stop('group values and index list are not match')
    g@values <- data.frame(indexID=1:length(values),values=values)
    g@indices <- index
  } else if (!is.null(index) && is.list(index) && !is.null(names(index))) {
    g@values <- data.frame(indexID=1:length(index),values=names(index))
    g@indices <- index
  }
  g
}
#-------
###################

# tbl <- c(122,'test','train') (example1 : tbl is vector)
# 
# tbl <- data.frame(matrix(NA,ncol=3,nrow=3)) (example 2: tbl is data.frame)
# tbl[,1] <- c(4,5,6)
# tbl[,2] <- rep('test',3)
# tbl[,3] <- rep('train',3)
# x: sdmdata
# the group for the specified ids will be changed from the first group to the second
.updateGroup <- function(x,tbl,sp=NULL) {
  if (is.vector(tbl)) {
    id <- as.numeric(tbl[1])
    g1 <- tbl[2]
    g2 <- tbl[3]
  } else {
    id <- as.numeric(tbl[,1])
    g1 <- tbl[,2]
    g2 <- tbl[,3]
  }
  gr <- .getGroupNames(x)
  for (g in gr) {
    for (i in seq_along(id)) {
      if (g1[i] %in% names(x@groups[[g]]@indices)) {
        if (id[i] %in% .getGroupIndex(x,g1[i])) {
          x@groups[[g]]@indices[[g1[i]]] <- x@groups[[g]]@indices[[g1[i]]][-which(x@groups[[g]]@indices[[g1[i]]] == id[i])]
          x@groups[[g]]@indices[[g2[i]]] <- c(x@groups[[g]]@indices[[g2[i]]] , id[i])
        }
      }
    }
  }
  x
}
#-----------
.getSpeciesNames <- function(d,n=NULL) {
  if (is.null(n)) d@species.names
  else d@species.names[d@species.names %in% n]
}
#-----
.getGroupLevels <- function(d,g=NULL) {
  if (!is.null(g)) {
    g <- g[1]
    if (!g %in% .getGroupNames(d)) stop(paste('group',g,'does not exist!'))
    as.character(d@groups[[g]]@values[,2])
  }
}
#-----
.getSpeciesDF <- function(d,sp,id) {
  if (missing(sp)) sp <- d@species.names[1]
  o <- list()
  for (s in sp) {
    o[[s]] <- data.frame(rID=id,value=NA)
    if (d@species[[s]]@type == 'Presence-Absence') {
      o[[s]][id %in% d@species[[s]]@presence,2] <- 1
      o[[s]][id %in% d@species[[s]]@absence,2] <- 0
    } else if (d@species[[s]]@type == 'Presence-Only') {
      o[[s]][id %in% d@species[[s]]@presence,2] <- 1
    } else if (d@species[[s]]@type == 'Abundance') {
      o[[s]][id %in% d@species[[s]]@abundance$rID,2] <- sapply(id[id %in% d@species[[s]]@abundance$rID],function(x) {d@species[[s]]@abundance[d@species[[s]]@abundance$rID == x,2]})
    } else if (d@species[[s]]@type == 'Absence-Only!') {
      o[[s]][id %in% d@species[[s]]@absence,2] <- 0
    } else if (d@species[[s]]@type == 'Abundance-constant!') {
      o[[s]][id %in% d@species[[s]]@abundance$rID,2] <- d@species[[s]]@abundance[1,2]
    }
  }
  o
}

#-----
.getSpeciesIndex <- function(d,sp=NULL) {
  id <- c()
  sp <- .getSpeciesNames(d,sp)
  for (i in seq_along(sp)) {
    id <- c(id,d@species[[sp[i]]]@presence,d@species[[sp[i]]]@absence,d@species[[sp[i]]]@abundance[,1],d@species[[sp[i]]]@background,d@species[[sp[i]]]@Multinomial[,1])
  }
  sort(unique(id))
}
#-----
.getGroupIndex=function(d,g=NULL) {
  if (!is.null(g)) {
    id <- c()
    for (gg in g) {
      gl <- unlist(lapply(unlist(strsplit(gg,':')),.trim))
      if (length(gl) > 1) {
        if (!gl[1] %in% .getGroupNames(d)) stop(paste('group',gl[1],'does not exist!'))
        if (!gl[2] %in% .getGroupNames(d,TRUE)) stop(paste('group level',gl[2],'does not exist!'))
        id <- c(id,d@groups[[gl[1]]]@indices[[gl[2]]])
      } else {
        if (!gl %in% c(.getGroupNames(d),.getGroupNames(d,TRUE))) stop(paste(gl,' is neither a group nor a group level!'))
        if (gl %in% .getGroupNames(d)) {
          for (gv in d@groups[[gl]]@values[,2]) id <- c(id,d@groups[[gl]]@indices[[gv]])
        } else {
          for (gn in .getGroupNames(d)) {
            if (gl %in% d@groups[[gn]]@values[,2]) id <- c(id,d@groups[[gn]]@indices[[gl]])
          }
        }
      }
    }
    unique(id)
  }
}
#-----
.getTimeIndex <- function(d,t=NULL) {
  if (!is.null(t)) {
    id <- c()
    #g <- d@groups
    #d@groups$training@values
    #sdm:::.getGroupNames(d)
    for (gg in t) {
      gl <- unlist(lapply(unlist(strsplit(gg,':')),.trim))
      if (length(gl) > 1) {
        if (!gl[1] %in% .getGroupNames(d)) stop(paste('group',gl[1],'does not exist!'))
        if (!gl[2] %in% .getGroupNames(d,TRUE)) stop(paste('group level',gl[2],'does not exist!'))
        id <- c(id,d@groups[[gl[1]]]@indices[[gl[2]]])
      } else {
        if (!gl %in% c(.getGroupNames(d),.getGroupNames(d,TRUE))) stop(paste(gl,' is neither a group nor a group level!'))
        if (gl %in% .getGroupNames(d)) {
          for (gv in d@groups[[gl]]@values[,2]) id <- c(id,d@groups[[gl]]@indices[[gv]])
        } else {
          for (gn in .getGroupNames(d)) {
            if (gl %in% d@groups[[gn]]@values[,2]) id <- c(id,d@groups[[gn]]@indices[[gl]])
          }
        }
      }
    }
    unique(id)
  }
}
#-----
.getIndex <- function(d,sp=NULL,groups=NULL,time=NULL) {
  id1 <- id2 <- NULL
  if (is.null(sp)) id1 <- d@features$rID
  else {
    id1 <- .getSpeciesIndex(d,sp)
  }
  
  if (!is.null(groups)) {
    id2 <- .getGroupIndex(d,groups)
  }
  
  if (is.null(id2)) id1
  else id1[id1 %in% id2]
}
#-----
.getGroupNames <- function(d,levels=FALSE) {
  if (levels) {
    if (!is.null(names(d@groups))) {
      nn <- c()
      for (n in names(d@groups)) nn <- c(nn,as.character(d@groups[[n]]@values[,2]))
      nn
    } else NULL
  } else names(d@groups)
}
#-----

.newGroup <- function(d,name,values=NULL,index=NULL) {
  d@groups[[name]] <- .newgroup(name,values,index)
  d
}
#-----
.getFeature <- function(d,n,type='l',id,...) {
  if (missing(id)) id <- .getIndex(d)
  rid <- which(d@features$rID %in% id)
  n <- n[1]
  type <- tolower(type)
  dot <- list(...)
  if (!n %in% d@features.name) stop('the variable does not exist!')
  if (type %in% c('l','linear')) {
    d@features[rid,n]
  } else if (type %in% c('q','quad','quadratic')) {
    .getFeature.quad(d@features[rid,n])
  } else if (type %in% c('c','cub','cubic')) {
    .getFeature.cubic(d$features[rid,n])
  } else if (type %in% c('poly')) {
    if ('degree' %in% names(dot)) degree <- dot[['degree']]
    else degree <- 3
    if ('raw' %in% names(dot)) raw <- dot[['raw']]
    else raw <- TRUE
    o <- .getFeature.poly(d@features[rid,n],degree=degree,raw=raw)
    colnames(o) <- paste(colnames(o),'.',n,sep='')
    o
  } else if (type %in% c('h','hing','hinge')) {
    if ('th' %in% names(dot)) th <- dot[['th']]
    else th <- NULL
    if (is.null(th)) {
      if ('species' %in% names(dot)) {
        if (is.numeric(dot[['species']]) && dot[['species']] <= length(d@species.names)) s <- dot[['species']]
        else if (is.character(dot[['species']]) && dot[['species']] %in% d@species.names) s <- which(d@species.names == dot[['species']])
        else stop('species is not identified!')
      } else {
        if (length(d@species.names) > 1) {
          s <- 1
          warning('to detect the threshold for the hinge feature, the species name is needed; since it is not specified the first species is used')
        } else s <- 1
      }
      s <- d@species.names[s]
      hP <- .getHingeParams(d@features[rid,n],y=.getSpeciesDF(d,sp=s,id=id)[[1]]$value)
      .getFeature.hinge(d@features[rid,n],increasing=hP$increasing,th=hP$threshold)
    } else .getFeature.hinge(d@features[rid,n],th=th)
    
  } else if (type %in% c('th','threshold')) {
    if ('th' %in% names(dot)) th <- dot[['th']]
    else th <- NULL
    if (is.null(th)) {
      if ('species' %in% names(dot)) {
        if (is.numeric(dot[['species']]) && dot[['species']] <= length(d@species.names)) s <- dot[['species']]
        else if (is.character(dot[['species']]) && dot[['species']] %in% d@species.names) s <- which(d@species.names == dot[['species']])
        else stop('species is not identified!')
      } else {
        if (length(d@species.names) > 1) {
          s <- 1
          warning('to detect the threshold for the threshold feature, the species name is needed; since it is not specified the first species is used')
        } else s <- 1
      }
      s <- d@species.names[s]
      thP <- .getThresholdParams(d@features[rid,n],y=.getSpeciesDF(d,sp=s,id=id)[[1]]$value)
      .getFeature.threshold(d@features[rid,n],th=thP$threshold,increasing=thP$increasing)
    } else .getFeature.threshold(d@features[rid,n],th=th)
    
    
    
  }
}
#-----------

if (!isGeneric('.addLog<-')) {
  setGeneric('.addLog<-', function(d,value)
    standardGeneric('.addLog<-'))
}
#---
setReplaceMethod('.addLog','sdmdata', 
                 function(d,value) {
                   d@errorLog <- c(d@errorLog,value)
                   d
                 }
)
#----------

.isBinomial <- function(x) {
  if (is.numeric(x)) {
    u <- unique(x)
    if (length(u) > 2) return(FALSE)
    else if (length(u) == 2) return(all(sort(u) == c(0,1)) | all(sort(u) == c(-1,1)))
    else return(u == 1)
  } else if (is.logical(x)) return(TRUE)
  else {
    x <- as.character(x)
    u <- unique(x)
    if (length(u) > 2) return(FALSE)
    else if (length(u) == 2) return(all(sort(u) == c('0','1')) | all(sort(u) == c('-1','1')))
    else return(u == '1')
  }
}
#----------

# detect the type of species data (i.e., pa, po, ab); whan all are 0 returns ao, and when variance is 0, returns ab_constant 
.speciesType <- function(x) {
  u <- unique(x)
  if (is.numeric(x)) {
    if (length(u) > 2) return('Abundance')
    else if (length(u) == 2) {
      if ((all(sort(u) == c(0,1)) || all(sort(u) == c(-1,1)))) return('Presence-Absence')
      else return('Abundance')
    } else {
      if (u == 1) return('Presence-Only')
      else if (u == 0 || u == -1) return('Absence-Only!') # ao is absence only! ONLY to detect and handle this kind of records
      else return('Abundance_constant!') # ab_constant to detect the sepcies data when the variance is 0
    } 
  } else if (is.logical(x)) {
    if (length(u) == 2) return('Presence-Absence')
    else {
      if (u) return('Presence-Only')
      else return('Absence-Only!')
    }
  } else {
    x <- as.character(x)
    u <- unique(x)
    if (length(u) > 2) return('Presence-Only')
    else if (length(u) == 2) {
      if (all(sort(u) == c('0','1')) || all(sort(u) == c('-1','1'))) return('Presence-Absence')
      else return('Presence-Only')
    } else return('Presence-Only')
  }
}

#----------
# check whether the names (vars) do exist in data (data.frame)
.varExist <-function(data,vars) {
  all(vars %in% names(data))
}
#----------

# get species dara.frame from the input data:
.getSpecies <- function(data,nsp,bg=FALSE,id.train=NULL,id.test=NULL) {
  species <- list()
  for (n in nsp) {
    if (!is.null(id.test)) {
      typ1 <- .speciesType(data[id.train,n])
      typ2 <- .speciesType(data[id.test,n])
      if (typ1 == typ2) typ <- typ1
      else stop('train and test data have different type (for example, one maybe presence-only while the other is presence-absence)!')
    } else typ <- .speciesType(data[,n])
    
    if (typ == 'Presence-Absence') {
      species[[n]] <- new('.species.data')
      species[[n]]@name <- n
      w <- as.numeric(data[,n])
      species[[n]]@presence <- data$rID[which(w == 1)]
      if (bg) {
        species[[n]]@background <- data$rID[which(w %in% c(0,-1))]
        species[[n]]@type <- 'Presence-Background'
      } else {
        species[[n]]@absence <- data$rID[which(w %in% c(0,-1))]
        species[[n]]@type <- typ
      }
    } else if (typ == 'Presence-Only') {
      if (is.numeric(data[,n]) || is.logical(data[,n])) {
        species[[n]] <- new('.species.data')
        species[[n]]@name <- n
        species[[n]]@presence <- data$rID
        species[[n]]@type <- typ
      } else {
        w <- as.character(data[,n])
        u <- unique(w)
        for (uu in u) {
          species[[uu]] <- new('.species.data')
          species[[uu]]@name <- uu
          species[[uu]]@presence <- data$rID[which(w == uu)]
          species[[uu]]@type <- typ
        }
      }
    } else if (typ == 'Abundance') {
      species[[n]] <- new('.species.data')
      species[[n]]@name <- n
      species[[n]]@abundance <- data.frame(rID=data$rID,abundance=data[,n])
      species[[n]]@type <- typ
    } else if (typ == 'Abundance_constant!') {
      species[[n]] <- new('.species.data')
      species[[n]]@name <- n
      species[[n]]@abundance <- data.frame(rID=data$rID,abundance=data[,n])
      species[[n]]@type <- typ
      warning(paste('for species',n,', the variance in abundance data is ZERO!'))
    } else if (typ == 'Multinomial') {
      species[[n]] <- new('.species.data')
      species[[n]]@name <- n
      species[[n]]@multinomial <- data.frame(rID=data$rID,name=data[,n])
      species[[n]]@type <- typ
    } else if (typ == 'Absence-Only!') {
      species[[n]] <- new('.species.data')
      species[[n]]@name <- n
      if (bg) {
        species[[n]]@background <- data$rID
        species[[n]]@type <- 'Background'
      } else {
        species[[n]]@absence <- data$rID
        species[[n]]@type <- typ
      }
      warning(paste('for species',n,', there is no presence record (all values are 0 or absence)!'))
    }
  }
  species
}

#----------
# remove duplicate records, and the rows that species columns contin NA OR all columns contain NA
.dataClean <- function(x,nsp) {
  rm.na <-0; rm.duplicate <- 0
  w <- nrow(x)
  x <- unique(x)
  if (nrow(x) < w) rm.duplicate <- w - nrow(x)
  
  w <- nrow(x)
  if (!missing(nsp)) {
    ww <- which(apply(x[which(colnames(x) %in% nsp),],1,function(x){all(is.na(x))}))
    if (length(ww) > 0) {
      x <- x[-ww,]
      rm.na <- w - nrow(x)
    }
  } else {
    ww <- which(apply(x,1,function(x){all(is.na(x))}))
    if (length(ww) > 0) {
      x <- x[-ww,]
      rm.na <- w - nrow(x)
    }
  }
  list(x,c(na=rm.na,duplicate=rm.duplicate))
}
#----------
# given a vector of colnames, their correponding col numbers are retrurned
.colNumber <- function(d,n) {
  unlist(lapply(n,function(x) which(colnames(d) == x)))
}
#----------
.char2time <- function(d,...) {
  if (length(list(...)) > 0) {
    tst <- try(as.POSIXct(d[1],...),silent=TRUE)
    if (!inherits(tst, "try-error") && !is.na(tst)) return(as.POSIXct(d,...))
    else {
      tst <- try(as.POSIXct(d[1]),silent=TRUE)
      if (!inherits(tst, "try-error") && !is.na(tst)) return(as.POSIXct(d))
      else {
        tst <- try(as.Date(d[1],...),silent=TRUE)
        if (!inherits(tst, "try-error") && !is.na(tst)) return(as.Date(d,...))
        else {
          tst <- try(as.Date(d[1]),silent=TRUE)
          if (!inherits(tst, "try-error") && !is.na(tst)) return(as.Date(d))
          else return(NA)
        }
      }
    }
  } else {
    tst <- try(as.POSIXct(d[1]),silent=TRUE)
    if (!inherits(tst, "try-error") && !is.na(tst)) return(as.POSIXct(d))
    else {
      tst <- try(as.Date(d[1]),silent=TRUE)
      if (!inherits(tst, "try-error") && !is.na(tst)) return(as.Date(d))
      else return(NA)
    }
  }
}
#-----------
#----
.where <- function(f, x) {
  vapply(x, f, logical(1))
}
#------
.int.to.numeric <- function(data) {
  w <- which(unlist(lapply(data,is.integer)))
  if (length(w) > 0) {
    for (i in w) data[,i] <- as.numeric(data[,i])
  }
  data
}
#----------
#--------------------

.which.is.coords <- function(n) {
  nxy <- NULL
  w <- tolower(n) %in% c('x','y','coords.x1','coords.x2','coords.x','coords.y','lon','long','longitude','lat','latitude')
  if (any(w)) {
    nxy <- n[w]
    if (length(nxy) == 2) {
      nxy <- c(nxy[tolower(nxy) %in% c('x','coords.x1','coords.x','lon','long','longitude')],nxy[tolower(nxy) %in% c('y','coords.x2','coords.y','lat','latitude')])
    } else nxy <- NULL
  }
  nxy
}
#---------------

.normalize <- function(x,except=NULL) {
  w <- !.where(is.factor,x)
  if (!is.null(except)) {
    w[except] <- FALSE
  }
  if (any(w)) {
    xx <- x
    for (i in seq_along(w)) {
      if (w[i]) {
        xx[,i] <- xx[,i] - mean(xx[,i],na.rm=TRUE)
        if (sd(x[,i],na.rm=TRUE) != 0) xx[,i] <- xx[,i] / sd(x[,i],na.rm=TRUE)
      }
    }
  }
  xx
}
#-----------

.speciesDetect <- function(data) {
  # to detect species columns with presence/absence data
  # also detect factors, and lon/lat coordinate columns
  nsp <-  nxy <- nFact <- nf <- nt <- NULL
  w <- which(unlist(lapply(data,.isBinomial)))
  if (length(w) == 0) {
    stop ('No species variable is detected, spcify the species variable in the formula or use an appropriate data structure...')
  } else {
    varNames <- colnames(data)
    nsp <- varNames[w]
    if (length(varNames) > length(nsp)) {
      nf <- varNames[-w]
      w <- which(unlist(lapply(data[,nf],function(x) class(x) %in% c("POSIXct","POSIXt","Date","yearmon","yearqtr"))))
      if (length(w) > 0) {
        nt <- nf[w]
        nf <- .excludeVector(nf,nt)
      }
      w <- which(unlist(lapply(data[,nf],function(x) class(x) %in% c('character','factor'))))
      if (length(w) > 0) {
        nFact <- nf[w]
        nf <- .excludeVector(nf,nFact)
      }
    }
    w <- tolower(varNames) %in% c('lon','long','longitude')
    if (any(w) && length(which(w) == 1)) {
      nx <- varNames[w]
      w <- tolower(varNames) %in% c('lat','latitude')
      if (any(w) && length(which(w) == 1)) nxy <- c(nx,varNames[w])
    }
  }
  list(nsp=nsp,nf=nf,nFact=nFact,nxy=nxy,nt=nt)
}
#--------

# create sdmdata object:
.createSdmdata <- function(train,formula=NULL,test,bg=NULL,crs=NULL,author=NULL,website=NULL,citation=NULL,help=NULL,description=NULL,date=NULL,license=NULL) {
  if (missing(test)) test <- NULL
  nFact <- nf <- nxy <- nsp <- ng <- nt <- ni <- NULL
  
  if (is.null(formula)) {
    nnFact <- nnf <- nnxy <- nnt <- NULL
    w <- .speciesDetect(train)
    if (!is.null(w$nFact)) {
      nFact <- w$nFact
      nnFact <- paste(paste('f(',w$nFact,')',sep=''),collapse='+')
    }
    if (!is.null(w$nf)) {
      nf <- w$nf
      nnf <- paste(w$nf,collapse='+')
    }
    if (!is.null(w$nxy)) {
      nxy <- w$nxy
      nnxy <- paste(paste('coords(',paste(w$nxy,collapse='+'),')',sep=''),collapse='+')
    }
    if (!is.null(w$nt)) {
      nt <- w$nt
      nnt <- paste(paste('time(',w$nt,')',sep=''),collapse='+')
    }
    formula <- as.formula(paste(paste(w$nsp,collapse="+"),'~',paste(c(nnf,nnFact,nnxy,nnt),collapse='+')),env = parent.frame())
  }
  
  exf <- .exFormula(formula,train)
  nall <- c(exf@vars,exf@species)
  
  
  if (!.varExist(train,nall)) stop('one or more specified variables in the formula do not exist in the train data!')
  
  nsp <- exf@species
  
  d <- new('sdmdata')
  
  d@sdmFormula <- exf
  
  w <- .dataClean(train,nsp)
  if (any(w[[2]] > 0)) {
    train <- w[[1]]
    ww <- c()
    if (w[[2]][1] > 0) .addLog(d) <- paste(w[[2]][1],'records with NA from the train data are removed')
    if (w[[2]][2] > 0) .addLog(d) <- paste(w[[2]][2],'duplicarted records from the train data are removed')
  }
  
  train$rID <- 1:nrow(train)
  
  train <- .int.to.numeric(train)
  
  
  if (!is.null(bg)) {
    w <- which(!nsp %in% colnames(bg))
    if (length(w) > 0) {
      nnsp <- nsp[w]
      nnsp <- matrix(0,nrow=nrow(bg),ncol=length(nnsp))
      colnames(nnsp) <- nsp[w]
      bg <- cbind(bg,nnsp)
    }
    bg <- as.data.frame(bg)
    
    if (!.varExist(data.frame(bg),nall)) stop('one or more predictor variables do not exist in the background data!')
    
    w <- .dataClean(bg,nsp)
    if (any(w[[2]] > 0)) {
      bg <- w[[1]]
      ww <- c()
      if (w[[2]][1] > 0) .addLog(d) <- paste(w[[2]][1],'records with NA from the background data are removed')
      if (w[[2]][2] > 0) .addLog(d) <- paste(w[[2]][2],'duplicarted records from the background data are removed')
    }
    
    bg$rID <- (nrow(train)+1):(nrow(bg)+nrow(train))
    train <- rbind(train[,c('rID',nall)],bg[,c('rID',nall)])
    bg <- bg$rID
  }
  
  if (!is.null(test)) {
    if (!.varExist(test,nall)) stop('one or more specified variables in the formula does not exist in the test data!')
    w <- .dataClean(test,nsp)
    if (any(w[[2]] > 0)) {
      test <- w[[1]]
      ww <- c()
      if (w[[2]][1] > 0) .addLog(d) <- paste(w[[2]][1],'records with NA from the test data are removed')
      if (w[[2]][2] > 0) .addLog(d) <- paste(w[[2]][2],'duplicarted records from the test data are removed')
    }
    
    test$rID <- (nrow(train)+1):(nrow(test)+nrow(train))
    
    if (!is.null(bg)) {
      w <- unlist(lapply(nsp,function(x) .speciesType(test[,x])))
      w <- unique(w)
      if (length(w) > 1) stop(paste('Independent test data has different types of records including',paste(w,collapse=', ')))
      else if (w == "Presence-Only") {
        d <- .newGroup(d,'training',index=list(train=train$rID,test=c(test$rID,bg)))
        cat('WARNING:\n Independent test dataset contains only presence records, so, background data (pseuso-absences) are used as Absence in the dataset!\n')
      } else d <- .newGroup(d,'training',index=list(train=train$rID,test=test$rID))
    } else d <- .newGroup(d,'training',index=list(train=train$rID,test=test$rID))
    test <- .int.to.numeric(test)
  } else {
    d <- .newGroup(d,'training',index=list(train=train$rID))
  }
  
  
  #-------
  
  if (is.null(nf) & is.null(nFact)) {
    nf <- exf@vars
    w <- .where(is.factor,train[,nf]) | .where(is.character,train[,nf])
    if (any(w)) nFact <- nf[w]
    
    if (!is.null(exf@model.terms)) {
      w <- unlist(lapply(exf@model.terms,class))
      if ('.factor' %in% w) {
        ww <- exf@model.terms[w == '.factor']
        if (length(ww) > 0) {
          for (i in seq_along(ww)) {
            w <- as.character(ww[[i]]@x)
            w <- .excludeVector(w,'+')
            nFact <- unique(c(nFact,w))
          }
        }
      }
    }
    
    nf <- .excludeVector(nf,nFact)
    
    if (!is.null(exf@data.terms)) {
      w <- unlist(lapply(exf@data.terms,class))
      if (".coord.vars" %in% w) nxy <- exf@data.terms[[which(w == ".coord.vars")]]@xy
      
      if (".grouping" %in% w) {
        ng <- c()
        ww <- exf@data.terms[which(w == ".grouping")]
        for (i in seq_along(ww)) ng <- c(ng,ww[[i]]@group.var)
        nf <- .excludeVector(nf,ng)
        nFact <- .excludeVector(nFact,ng)
      }
      
      if ('.time' %in% w) {
        nt <- c()
        ww <- exf@data.terms[which(w == ".time")]
        for (i in seq_along(ww)) nt <- c(nt,as.character(ww[[i]]@terms[[1]]))
        nf <- .excludeVector(nf,nt)
        nFact <- .excludeVector(nFact,nt)
      }
      
      if ('.Info' %in% w) {
        ni <- c()
        ww <- exf@data.terms[which(w == ".Info")]
        for (i in seq_along(ww)) ni <- c(ni,ww[[i]]@names)
        nf <- .excludeVector(nf,ni)
        nFact <- .excludeVector(nFact,ni)
      }
      
    } else {
      w <- !colnames(train) %in% c(nall,'rID')
      if (any(w)) {
        nxy <- .which.is.coords(colnames(train)[w])
        if (!is.null(test) && !.varExist(test,nxy)) nxy <- NULL
        ww <- unlist(lapply(which(w),function(x) class(train[,x]))) %in% c("POSIXct","POSIXt","Date","yearmon","yearqtr")
        if (any(ww)) {
          nt <- colnames(train)[w[which(ww)]]
          nf <- .excludeVector(nf,nt)
          nFact <- .excludeVector(nFact,nt)
        }
      }
    }
    nf <- .excludeVector(nf,nxy)
    nFact <- .excludeVector(nFact,nxy)
  }
  #-------
  nall <- c(nsp,nf,nFact,nxy,ng,ni,nt,'rID')
  if (!is.null(test)) {
    train <- rbind(train[,nall],test[,nall])
    rm(test)
    species <- .getSpecies(train,nsp,bg=!is.null(bg),id.train = d@groups$training@indices$train,id.test = d@groups$training@indices$test)
  } else {
    train <- train[,nall]
    species <- .getSpecies(train,nsp,bg=!is.null(bg))
  }
  #----
  
  if (!is.null(ng)) {
    for (n in ng) {
      ww <- as.character(train[,n])
      u <- unique(ww)
      if (length(u) == 1) warning(paste('the grouping variable',n,'is ignored; it is constant!'))
      else {
        w <- list()
        for (uu in u) {
          w[[uu]] <- train$rID[which(ww == uu)]
        }
        d <- .newGroup(d,n,index=w)
      }
    }
  }
  #------
  if (!is.null(c(nxy,ni,nt,website,help,description,date,license)) || !is.null(c(citation,author))) {
    d@info <- new('.info')
    if (!is.null(nxy)) {
      d@info@coords <- as.matrix(train[,c('rID',nxy)])
      if (!is.null(crs) && inherits(crs,'CRS')) d@info@crs <- crs
    }
    if (!is.null(ni)) d@info@info <- train[,c('rID',ni)]
    
    if (!is.null(nt)) {
      dt <- data.frame(rID=train$rID)
      for (n in nt) {
        if ((class(train[,n]) %in% c("POSIXct","POSIXt","Date","yearmon","yearqtr"))[1]) {
          dt <- cbind(dt,train[,n])
        } else {
          w <- unlist(lapply(exf@data.terms,class))
          w <- exf@data.terms[which(w == ".time")]
          w <- w[[which(unlist(lapply(w,function(x) x@terms[[1]] == n)))]]
          if (length(w@terms) > 1 && is.null(names(w@terms)) && length(names(w@terms)) == length(w@terms)) {
            w <- do.call(.char2time,c(list(d=train[,n]),w@terms[2:length(w@terms)]))
            if ((class(w) %in% c("POSIXct","POSIXt","Date","yearmon","yearqtr"))[1]) dt[,n] <- w
            else warning(paste('a time-based format is not detected for variable',n,", so it is IGNORED; it must have a detectable character format, or being one of time-based classes including: 'POSIXct', 'POSIXt', 'Date', 'yearmon','yearqtr'"))
          } else {
            w <- .char2time(train[,n])
            if ((class(w) %in% c("POSIXct","POSIXt","Date","yearmon","yearqtr"))[1]) dt[,n] <- w
            else warning(paste('a time-based format is not detected for variable',n,", so it is IGNORED; it must have a detectable character format, or being one of time-based classes including: 'POSIXct', 'POSIXt', 'Date', 'yearmon','yearqtr'"))
          }
        }
      }
      if (ncol(train) > 1) d@info@time <- dt
    }
    
    if (!is.null(c(website,help,description,date,license)) || !is.null(c(citation,author))) {
      d@info@metadata <- .newMetadata(authors=author,web=website,cit=citation,desc=description,date=date,license=license,help=help)
    }
  }
  #----------
  d@features.name <- c(nf,nFact)
  if (!is.null(nFact)) d@factors <- nFact
  d@species <- species
  d@species.names <- names(species)
  for (n in nFact) train[,n] <- factor(train[,n])
  if (!is.null(d@features.name)) d@features <- train[,c('rID',nf,nFact)]
  d
}
#-------
# 
# # create sdmdata object:
# .createSdmdata_base <- function(train,test=NULL,nsp=NULL,nf=NULL,nFact=NULL,nxy=NULL,ng=NULL,nt=NULL,ni=NULL,crs=NULL,author=NULL,website=NULL,citation=NULL,help=NULL,description=NULL,date=NULL,license=NULL) {
# #   nnFact <- nnf <- nnxy <- nnt <- NULL
# #   if (!is.null(nFact)) nnFact <- paste(paste('f(',nFact,')',sep=''),collapse='+')
# #   if (!is.null(nf)) nnf <- paste(nf,collapse='+')
# #   if (!is.null(nxy)) nnxy <- paste(paste('coords(',paste(nxy,collapse='+'),')',sep=''),collapse='+')
# #   if (!is.null(nt)) nnt <- paste(paste('time(',nt,')',sep=''),collapse='+')
# #   formula <- as.formula(paste(paste(w$nsp,collapse="+"),'~',paste(c(nnf,nnFact,nnxy,nnt),collapse='+')),env = parent.frame())
# #   
#   #exf <- .exFormula(formula,train)
#   #nall <- c(exf@vars,exf@species)
#   nall <- c(nsp,nf,nFact,nxy,ng,ni,nt)
#   if (!.varExist(train,nall)) stop('one or more specified variables in the formula does not exist in the train data!')
#   
#   d <- new('sdmdata')
#   
#   #d@sdmFormula <- exf
#   
#   #nsp <- exf@species
#   
#   w <- .dataClean(train,nsp)
#   if (any(w[[2]] > 0)) {
#     train <- w[[1]]
#     ww <- c()
#     if (w[[2]][1] > 0) .addLog(d) <- paste(w[[2]][1],'records with NA from the train data are removed')
#     if (w[[2]][2] > 0) .addLog(d) <- paste(w[[2]][2],'duplicarted records from the train data are removed')
#   }
#   
#   train$rID <- 1:nrow(train)
#   
#   train <- .int.to.numeric(train)
#   
#   if (!is.null(test)) {
#     if (!.varExist(test,nall)) stop('one or more specified variables in the formula does not exist in the test data!')
#     w <- .dataClean(test,nsp)
#     if (any(w[[2]] > 0)) {
#       test <- w[[1]]
#       ww <- c()
#       if (w[[2]][1] > 0) .addLog(d) <- paste(w[[2]][1],'records with NA from the test data are removed')
#       if (w[[2]][2] > 0) .addLog(d) <- paste(w[[2]][2],'duplicarted records from the test data are removed')
#     }
#     
#     test$rID <- (nrow(train)+1):(nrow(test)+nrow(train))
#     d <- .newGroup(d,'training',index=list(train=train$rID,test=test$rID))
#     test <- .int.to.numeric(test)
#   } else {
#     d <- .newGroup(d,'training',index=list(train=train$rID))
#   }
#   #-------
#   
#   w <- !colnames(train) %in% c(nall,'rID')
#   if (any(w)) {
#     if (is.null(nxy)) {
#       nxy <- .which.is.coords(colnames(train)[w])
#       if (!is.null(test) && !.varExist(test,nxy)) nxy <- NULL
#     }
#     
#     ww <- unlist(lapply(which(w),function(x) class(train[,x]))) %in% c("POSIXct","POSIXt","Date","yearmon","yearqtr")
#     if (any(ww)) {
#       nt <- colnames(train)[w[which(ww)]]
#     }
#   }
#   nf <- .excludeVector(nf,c(nt,nxy))
#   nFact <- .excludeVector(nFact,c(nt,nxy))
#   
#   #-------
#   nall <- c(nsp,nf,nFact,nxy,ng,ni,nt,'rID')
#   if (!is.null(test)) {
#     train <- rbind(train[,nall],test[,nall])
#     rm(test)
#     species <- .getSpecies(train,nsp,id.train = d@groups$training@indices$train,id.test = d@groups$training@indices$test)
#   } else {
#     train <- train[,nall]
#     species <- .getSpecies(train,nsp)
#   }
#   #----
#   
#   if (!is.null(ng)) {
#     for (n in ng) {
#       ww <- as.character(train[,n])
#       u <- unique(ww)
#       if (length(u) == 1) warning(paste('the grouping variable',n,'is ignored; it is constant!'))
#       else {
#         w <- list()
#         for (uu in u) {
#           w[[uu]] <- train$rID[which(ww == uu)]
#         }
#         d <- .newGroup(d,n,index=w)
#       }
#     }
#   }
#   #------
#   if (!is.null(c(nxy,ni,nt,website,help,description,date,license)) || !is.null(c(citation,author))) {
#     d@info <- new('.info')
#     if (!is.null(nxy)) {
#       d@info@coords <- as.matrix(train[,c('rID',nxy)])
#       if (!is.null(crs) && inherits(crs,'CRS')) d@info@crs <- crs
#     }
#     if (!is.null(ni)) d@info@info <- train[,c('rID',ni)]
#     
#     if (!is.null(nt)) {
#       dt <- data.frame(rID=train$rID)
#       for (n in nt) {
#         if ((class(train[,n]) %in% c("POSIXct","POSIXt","Date","yearmon","yearqtr"))[1]) {
#           dt <- cbind(dt,train[,n])
#         } else {
#           w <- .char2time(train[,n])
#           if ((class(w) %in% c("POSIXct","POSIXt","Date","yearmon","yearqtr"))[1]) dt[,n] <- w
#           else warning(paste('a time-based format is not detected for variable',n,", so it is IGNORED; it must have a detectable character format, or being one of time-based classes including: 'POSIXct', 'POSIXt', 'Date', 'yearmon','yearqtr'"))
#         }
#       }
#       if (ncol(train) > 1) d@info@time <- dt
#     }
#     
#     if (!is.null(c(website,help,description,date,license)) || !is.null(c(citation,author))) {
#       d@info@metadata <- .newMetadata(authors=author,web=website,cit=citation,desc=description,date=date,license=license,help=help)
#     }
#   }
#   #----------
#   d@species <- species
#   d@species.names <- names(species)
#   
#   if (!is.null(c(nf,nFact))) {
#     d@features.name <- c(nf,nFact)
#     if (!is.null(nFact)) d@factors <- nFact
#     for (n in nFact) train[,n] <- factor(train[,n])
#     d@features <- train[,c('rID',nf,nFact)]
#   } 
#   d
# }

#------

.Extract <- function(x,cells,factors) {
  n <- names(x)
  
  if (length(factors) == 1) {
    x2 <- values(x[[factors]])[cells]
  } else x2 <- values(x[[factors]])[cells,]
  
  if (length(n) > length(factors)) {
    x1 <- x[[-factors]][cells]
    d <- data.frame(x1,x2)
    colnames(d) <- c(n[-factors],n[factors])
  } else {
    d <- data.frame(x2)
    colnames(d) <- n[factors]
  }
  
  for (i in 1:length(factors)) d[,i] <- factor(d[,i])
  return(d )
}
#---------
#--- pseudo-absence based on random distribution in geographic space:
# !is.null(p) : removes the points that are located in the locations with presence record
.pseudo_gRandom <- function(preds,n=1000,p=NULL) {
  s <- sampleRandom(preds,n,cells=TRUE,xy=TRUE)
  if (!is.null(p) && ncol(p) == 2) {
    p.cells <- cellFromXY(preds,p)
    if (length(p.cells) > 0) {
      s <- s[which(!s[,'cell'] %in% p.cells),]
    }
  }
  s[,-1]
}
#----------
.pseudo_eRandom <- function(preds,n=1000,p=NULL) {
  #
}
#----------
.pseudo_gDist <- function(preds,n=1000,p=NULL) {
  #
}
#----------
.pseudo_eDist <- function(preds,n=1000,p=NULL) {
  #
}
#----------
.pseudo <- function(preds,n=1000,method='gRandom',p=NULL) {
  if (method == 'gRandom') .pseudo_gRandom(preds,n=n,p=p)
  else if (method == 'eRandom') .pseudo_eRandom(preds,n=n,p=p)
  else if (method == 'gDist') .pseudo_gDist(preds,n=n,p=p)
  else if (method == 'eDist') .pseudo_eDist(preds,n=n,p=p)
}



if (!isGeneric("sdmData")) {
  setGeneric("sdmData", function(formula, train, test, predictors,bg, filename, crs,...)
    standardGeneric("sdmData"))
}


setMethod('sdmData', signature(train='data.frame',predictors='missing'), 
          function(formula,train,test=NULL,predictors,bg=NULL,filename=NULL,crs=NULL,...) {
            if(missing(test)) test <- NULL
            if(missing(filename)) filename <- NULL
            if(missing(crs)) crs <- NULL
            if(missing(formula)) formula <- NULL
            if(missing(bg)) bg <- NULL
            #---
            if (!.sdmOptions$getOption('sdmLoaded')) .addMethods()
            #---
            dot <- list(...)
            n <- tolower(names(dot))
            for (i in seq_along(n)) {
              if (any(!is.na(pmatch(c("aut"),n[i])))) n[i] <- 'author'
              else if (any(!is.na(pmatch(c("web"),n[i])))) n[i] <- 'website'
              else if (any(!is.na(pmatch(c("cit"),n[i])))) n[i] <- 'citation'
              else if (any(!is.na(pmatch(c("hel"),n[i])))) n[i] <- 'help'
              else if (any(!is.na(pmatch(c("des"),n[i])))) n[i] <- 'description'
              else if (any(!is.na(pmatch(c("dat"),n[i])))) n[i] <- 'date'
              else if (any(!is.na(pmatch(c("lic"),n[i])))) n[i] <- 'license'
            }
            names(dot) <- n
            author <- dot[['author']]
            website <- dot[['website']]
            citation <- dot[['citation']]
            help <- dot[['help']]
            description <- dot[['description']]
            date <- dot[['date']]
            license <- dot[['license']]
            
            .createSdmdata(train = train, formula = formula, test = test,bg=bg,crs = crs,author=author,website=website,citation=citation,help=help,description=description,date=date,license=license)
          }
)

#-------  

setMethod('sdmData', signature(formula='data.frame',train='formula',predictors='missing'), 
          function(formula,train,test=NULL,predictors,bg=NULL,filename=NULL,crs=NULL,...) {
            # to make it user friendly
            if(missing(test)) test <- NULL
            if(missing(filename)) filename <- NULL
            if(missing(crs)) crs <- NULL
            if(missing(bg)) bg <- NULL
            sdmData(formula=train, train=formula,test=test,bg=bg,filename=filename,crs=crs,...)
          }
)

setMethod('sdmData', signature(formula='data.frame',train='missing',predictors='missing'), 
          function(formula,train,test=NULL,predictors,bg=NULL,filename=NULL,crs=NULL,...) {
            # to make it user friendly
            if(missing(test)) test <- NULL
            if(missing(filename)) filename <- NULL
            if(missing(crs)) crs <- NULL
            if(missing(bg)) bg <- NULL
            sdmData(train=formula,test=test,bg=bg,filename=filename,crs=crs,...)
          }
)

setMethod('sdmData', signature(train='SpatialPoints',predictors='missing'), 
          function(formula,train,test=NULL,predictors,bg=NULL,filename=NULL,crs=NULL,...) {
            
            if(missing(test)) test <- NULL
            if(missing(filename)) filename <- NULL
            if(missing(crs)) crs <- NULL
            if(missing(bg)) bg <- NULL
            
            #---
            if (!.sdmOptions$getOption('sdmLoaded')) .addMethods()
            #---
            
            nxy <- coordnames(train)
            
            if (!is.null(test)) {
              if (class(train) == class(test)) {
                nxyt <- coordnames(test)
                if (as.character(class(test) == 'SpatialPoints')) test <- data.frame(SPECIES=rep(1,length(test)),as(test,'data.frame'))
                else test <- as(test,'data.frame')
                if (all(nxy != nxyt)) {
                  colnames(test)[unlist(lapply(nxyt,function(x) which(colnames(test) ==x)))] <- nxy
                } 
              }
            }
            
            if (!is.na(proj4string(train))) crs <- CRS(proj4string(train))
            if (as.character(class(train) == 'SpatialPoints')) train <- data.frame(SPECIES=rep(1,length(train)),as(train,'data.frame'))
            else train <- as(train,'data.frame')
            
            if (!missing(formula)) {
              if (!all(nxy %in% all.vars(formula))) {
                if ('.' %in% all.vars(formula)) {
                  ww <- .exFormula(formula,train)
                  nw <- .excludeVector(colnames(train),c(ww@species,nxy))
                  if (length(nw) > 0) {
                    w <- as.character(.exFormula(formula,train,FALSE)@species)
                    if (length(w) == 0) formula[[2]] <- terms.formula(formula,data=train[,nw])[[2]]
                    else {
                      formula[[3]] <- terms.formula(formula,data=train[,nw])[[3]]
                      if (colnames(train)[1] == 'SPECIES') {
                        colnames(train)[1] <- w[1]
                        if (!is.null(test) && colnames(test)[1] == 'SPECIES') colnames(test)[1] <- w[1]
                      }
                    }
                     
                  }
                }
                
                formula <- update(formula,as.formula(paste('~ . + coods(',paste(nxy,collapse='+'),')',sep='')))
              }
            } else formula <- as.formula(paste('~ . + coods(',paste(nxy,collapse='+'),')',sep=''))
            
            dot <- list(...)
            n <- tolower(names(dot))
            for (i in seq_along(n)) {
              if (any(!is.na(pmatch(c("aut"),n[i])))) n[i] <- 'author'
              else if (any(!is.na(pmatch(c("web"),n[i])))) n[i] <- 'website'
              else if (any(!is.na(pmatch(c("cit"),n[i])))) n[i] <- 'citation'
              else if (any(!is.na(pmatch(c("hel"),n[i])))) n[i] <- 'help'
              else if (any(!is.na(pmatch(c("des"),n[i])))) n[i] <- 'description'
              else if (any(!is.na(pmatch(c("dat"),n[i])))) n[i] <- 'date'
              else if (any(!is.na(pmatch(c("lic"),n[i])))) n[i] <- 'license'
            }
            names(dot) <- n
            author <- dot[['author']]
            website <- dot[['website']]
            citation <- dot[['citation']]
            help <- dot[['help']]
            description <- dot[['description']]
            date <- dot[['date']]
            license <- dot[['license']]
            
            .createSdmdata(train = train, formula = formula, test = test,bg=bg,crs = crs,author=author,website=website,citation=citation,help=help,description=description,date=date,license=license) 
            
          }
)


setMethod('sdmData', signature(train='SpatialPoints',predictors='Raster'), 
          function(formula,train,test=NULL,predictors,bg=NULL,filename=NULL,crs=NULL,...) {
            if(missing(test)) test <- NULL
            if(missing(filename)) filename <- NULL
            if(missing(crs)) crs <- NULL
            if(missing(bg)) bg <- NULL
            #---
            if (!.sdmOptions$getOption('sdmLoaded')) .addMethods()
            #---
            testc <- as.character(class(test))
            trainc <- as.character(class(train))
            
            trainSP <-NULL
            errLog <- list()
            
            nxy <- coordnames(train)
            wF <- is.factor(predictors)
           
            if (!is.null(test)) {
              if (inherits(test,'SpatialPoints')) {
                nxyt <- coordnames(test)
                if (testc == 'SpatialPointsDataFrame') test <- as(test,'data.frame')
                else test <- data.frame(coordinates(test))
                if (all(nxy != nxyt)) {
                  colnames(test)[unlist(lapply(nxyt,function(x) which(colnames(test) ==x)))] <- nxy
                } 
              } else if (!is.data.frame(test)) stop('test data should be a data.frame or in the same class as the train data!')
              
              if (nxy[1] %in% colnames(test) & nxy[2] %in% colnames(test)) {
                cells <- cellFromXY(predictors,test[,nxy])
                cNA <- is.na(cells)
                if (any(cNA)) {
                  if (all(cNA)) stop('Test dataset has no overlap with the predictors...!')
                  wNA <- which(cNA)
                  test <- test[-wNA,]
                  cells <- cells[-wNA]
                  errLog <- c(errLog,paste(length(wNA),'records were removed from the test dataset because of no overlap with the predictors.'))
                  rm(cNA,wNA)
                }
                rm(cNA)
                if (!any(wF)) test.p <- data.frame(predictors[cells])
                else test.p <- .Extract(predictors,cells,which(wF))
                rm(cells)
                colnames(test.p) <- names(predictors)
                
                w <- colnames(test) %in% colnames(test.p)
                if (any(w)) {
                  test <- test[,colnames(test)[!w]]
                  errLog <- c(errLog,paste('WARNING: The variables',colnames(test)[w],'were removed from the test dataset as they exist in the predictors as well.'))
                }
              } else stop('the coordinates name are not match in the train and test datasets!')
            }
            
            if (!is.na(proj4string(train))) crs <- CRS(proj4string(train))
            
            if (trainc == 'SpatialPointsDataFrame') train <- as(train,'data.frame')
            else train <- coordinates(train)
            
            cells <- cellFromXY(predictors,train[,nxy])
            cNA <- is.na(cells)
            if (any(cNA)) {
              if (all(cNA)) stop('Train data has no overlap with the predictors...!')
              wNA <- which(cNA)
              train <- train[-wNA,]
              cells <- cells[-wNA]
              errLog <- c(errLog,paste(length(wNA),'records were removed from the train dataset because of no overlap with the predictors.'))
            }
            rm(cNA)
            
            if (!any(wF)) train.p <- data.frame(predictors[cells])
            else train.p <- .Extract(predictors,cells,which(wF))
            rm(cells)
            colnames(train.p) <- names(predictors)
            w <- colnames(train) %in% colnames(train.p)
            if (any(w)) {
              train <- train[,colnames(train)[!w]]
              errLog <- c(errLog,paste('WARNING: The variables',colnames(train)[w],'were removed from the train dataset as they exist in the predictors as well.'))
            }
            
            if (trainc == 'SpatialPointsDataFrame') {
              train <- data.frame(train,train.p)
              rm(train.p)
              if (!is.null(test)) {
                test <- data.frame(test,test.p)
                rm(test.p)
              }
              
              if (!missing(formula)) {
                if (!all(nxy %in% all.vars(formula))) {
                  if ('.' %in% all.vars(formula)) {
                    ww <- .exFormula(formula,train)
                    nw <- .excludeVector(colnames(train),c(ww@species,nxy))
                    w <- as.character(.exFormula(formula,train,FALSE)@species)
                    if (length(w) == 0) formula[[2]] <- terms.formula(formula,data=train[,nw])[[2]]
                    else formula[[3]] <- terms.formula(formula,data=train[,nw])[[3]]
                  }
                  formula <- update(formula,as.formula(paste('~ . + coods(',paste(nxy,collapse='+'),')',sep='')))
                }
              } else formula <- as.formula(paste('~ . + coods(',paste(nxy,collapse='+'),')',sep=''))
            } else {
              if (!missing(formula)) {
                ww <- as.character(.exFormula(formula,train.p,FALSE)@species)
                if (length(ww) > 1) {
                  warning('While SpatialPoints can be used for only 1 species, more names are defined in the formula! The first name is considered!')
                  errLog <- c(errLog,'WARNING: While SpatialPoints can be used for only 1 species, more names are defined in the formula! The first name is considered!')
                  w <- ww[1]
                } else if (length(ww) == 0) w <- 'SPECIES'
                
                if (!all(nxy %in% all.vars(formula))) {
                  
                  if ('.' %in% all.vars(formula)) {
                    if (length(ww) == 0) formula[[2]] <- terms.formula(formula,data=train.p)[[2]]
                    else formula[[3]] <- terms.formula(formula,data=train.p)[[3]]
                  }
                  formula <- update(formula,as.formula(paste('~ . + coods(',paste(nxy,collapse='+'),')',sep='')))
                }
              } else {
                formula <- as.formula(paste('SPECIES ~ . + coods(',paste(nxy,collapse='+'),')',sep=''))
                w <- 'SPECIES'
              }
              
              train <- data.frame(SPECIES=rep(1,nrow(train)),train)
              colnames(train)[1] <- w
              train <- data.frame(train,train.p)
              rm(train.p)
              
              if (!is.null(test)) {
                test <- data.frame(SPECIES=rep(1,nrow(test)),test,test.p)
                colnames(test)[1] <- w
                rm(test.p)
              }
            } 
            
            if (!is.null(bg)) {
              if (is.list(bg)) {
                nbg <- names(bg)
                nbg <- .pmatch(nbg,c('n','method','remove'))
                if ('n' %in% nbg) n <- bg[['n']]
                else n <- 1000
                if ('method' %in% nbg) {
                  if (.pmatch(bg[['method']],c('gRandom','random','rnd')) %in% c('gRandom','random','rnd')) {
                    m <- 'gRandom'
                  } else if (.pmatch(bg[['method']],c('eRandom','envrandom','ernd')) %in% c('eRandom','envrandom','ernd')) {
                    m <- 'eRandom'
                  } else if (.pmatch(bg[['method']],c('gDistance','geo')) %in% c('gDistance','geo')) {
                    m <- 'gDist'
                  } else if (.pmatch(bg[['method']],c('eDistance','environ','envDist')) %in% c('eDistance','environ','envDist')) {
                    m <- 'eDist'
                  }
                } else m <- 'gRandom'
                if ('remove' %in% nbg && is.logical(bg[['remove']])) r <- bg[['remove']]
                else r <- FALSE
                
                bg <- .pseudo(predictors,n=n,method = m,p = if (r) train[,nxy] else NULL)
                colnames(bg)[1:2] <- nxy
              } else if (is.numeric(bg)) {
                bg <- .pseudo(predictors,n=bg,method = 'gRandom',p = NULL)
                colnames(bg)[1:2] <- nxy
              } else if (!is.data.frame(bg)) bg <- NULL
            }
            
            dot <- list(...)
            n <- tolower(names(dot))
            for (i in seq_along(n)) {
              if (any(!is.na(pmatch(c("aut"),n[i])))) n[i] <- 'author'
              else if (any(!is.na(pmatch(c("web"),n[i])))) n[i] <- 'website'
              else if (any(!is.na(pmatch(c("cit"),n[i])))) n[i] <- 'citation'
              else if (any(!is.na(pmatch(c("hel"),n[i])))) n[i] <- 'help'
              else if (any(!is.na(pmatch(c("des"),n[i])))) n[i] <- 'description'
              else if (any(!is.na(pmatch(c("dat"),n[i])))) n[i] <- 'date'
              else if (any(!is.na(pmatch(c("lic"),n[i])))) n[i] <- 'license'
            }
            names(dot) <- n
            author <- dot[['author']]
            website <- dot[['website']]
            citation <- dot[['citation']]
            help <- dot[['help']]
            description <- dot[['description']]
            date <- dot[['date']]
            license <- dot[['license']]
            
            d <- .createSdmdata(train = train, formula = formula, test = test,bg=bg,crs = crs,author=author,website=website,citation=citation,help=help,description=description,date=date,license=license)
            
            if (length(errLog) > 0) {
              for (i in seq_along(errLog)) .addLog(d) <- errLog[[i]]
            }
            d
          }
)

