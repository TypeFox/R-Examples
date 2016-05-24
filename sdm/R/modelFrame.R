# Author: Babak Naimi, naimi.b@gmail.com
# Date :  September 2015
# Version 1.0
# Licence GPL v3

.getFeature.linear <- function(x) {
  x
}


.getFeature.quad <- function(x) {
  x * x
}

.getFeature.cubic <- function(x) {
  x * x * x
}

.getFeature.poly <- function(x,degree=3,raw=TRUE) {
  d <- as.data.frame(poly(x,degree=degree,raw=raw))
  colnames(d) <- paste('poly',1:degree,sep='')
  d
}



#-------
.getFeature.hinge <- function(x,th,increasing) {
  if (increasing) {
    .hinge(x,th)
  } else {
    .invhinge(x,th)
  }
}
#-------
.getFeature.threshold <- function(x,th,increasing) {
  if (increasing) {
    .thresh(x,th)
  } else {
    .invthresh(x,th)
  }
}

#-------
.hinge <- function(x,th) {
  ifelse(x < th,0,(x - th) / (max(x,na.rm=TRUE) - th))
}

.invhinge <- function(x,th) {
  ifelse(x > th,0,1 - ((x - min(x,na.rm=TRUE)) / (th - min(x,na.rm=TRUE))))
}
#---------
.invthresh <- function(x,th) {
  ifelse(x > th,0,1)
}

.thresh <- function(x,th) {
  ifelse(x <= th,0,1)
}
#------
.Jackard <- function(x,y) {
  sum(apply(data.frame(x,y),1,min)) / sum(apply(data.frame(x,y),1,max))
}
#----------

.getFeature.product <- function(data) {
  apply(data,1,function(x) {
    xx <- 1
    for (i in seq_along(x))  xx <- xx * x[i]
    xx
  })
}
#----------

.getHingeParams <- function(x,y,method='deviance') {
  n <- min(20,length(x))
  th <- as.vector(quantile(x,0:n/n))
  
  if (coefficients(lm(y~x))[2] < 0) {
    j <- rep(NA,length(th))
    if (method %in% c('deviance','dev','Deviance','Dev','devience')) {
      for (i in 1:length(th)) j[i] <- .deviance_binomial(y,.invhinge(x,th[i]))
      w <- which(j == min(j,na.rm=TRUE))
    } else if (method %in% c('jackard','J','jac','Jackard')) {
      for (i in 1:length(th)) j[i] <- .Jackard(y,.invhinge(x,th[i]))
      w <- which(j == max(j,na.rm=TRUE))
    }
    list(increasing=FALSE,threshold=th[w[1]])
  } else {
    j <- rep(NA,length(th))
    if (method %in% c('deviance','dev','Deviance','Dev','devience')) {
      for (i in 1:length(th)) j[i] <- .deviance_binomial(y,.hinge(x,th[i]))
      w <- which(j == min(j,na.rm=TRUE))
    } else if (method %in% c('jackard','J','jac','Jackard')) {
      for (i in 1:length(th)) j[i] <- .Jackard(y,.hinge(x,th[i]))
      w <- which(j == max(j,na.rm=TRUE))
    }
    list(increasing=TRUE,threshold=th[w[length(w)]])
  }
}
#-----------

.getThresholdParams <- function(x,y,method='deviance') {
  n <- min(20,length(x))
  th <- as.vector(quantile(x,0:n/n))
  
  if (coefficients(lm(y~x))[2] < 0) {
    j <- rep(NA,length(th))
    if (method %in% c('deviance','dev','Deviance','Dev','devience')) {
      for (i in 1:length(th)) j[i] <- .deviance_binomial(y,.invthresh(x,th[i]))
      w <- which(j == min(j,na.rm=TRUE))
    } else if (method %in% c('jackard','J','jac','Jackard')) {
      for (i in 1:length(th)) j[i] <- .Jackard(y,.invthresh(x,th[i]))
      w <- which(j == max(j,na.rm=TRUE))
    }
    list(increasing=FALSE,threshold=th[w[1]])
  } else {
    j <- rep(NA,length(th))
    if (method %in% c('deviance','dev','Deviance','Dev','devience')) {
      for (i in 1:length(th)) j[i] <- .deviance_binomial(y,.thresh(x,th[i]))
      w <- which(j == min(j,na.rm=TRUE))
    } else if (method %in% c('jackard','J','jac','Jackard')) {
      for (i in 1:length(th)) j[i] <- .Jackard(y,.thresh(x,th[i]))
      w <- which(j == max(j,na.rm=TRUE))
    }
    list(increasing=TRUE,threshold=th[w[length(w)]])
  }
}

#--------

# test the hinge feature against the null model; a vector equal to the number of simulation is returned,
# lower values close to 0 is better, and higher, i.e., closer to 1 implies random or worse than random

.hingeNull_binomial <- function(x,y,n=99,method='mont') {
  m1 <- .getHingeParams(x,y)
  o <- rep(NA,n)
  if (method %in% c('mont','Mont','M','m','Monte Carlo','monte','Monte')) yn <- replicate(n,sample(y))
  else if (method %in% c('boot','Boot','B','b','bootstrp','Bootstrap')) yn <- replicate(n,sample(y,length(y),replace=TRUE))
  
  if (m1$increasing) {
    h <- .hinge(x,m1$threshold)
    b <- .deviance_binomial(y,h)
    for (i in 1:n) o[i] <- .deviance_binomial(yn[,i],.hinge(x,m1$threshold))
  } else {
    h <- .invhinge(x,m1$threshold)
    b <- .deviance_binomial(y,h)
    for (i in 1:n) o[i] <- .deviance_binomial(yn[,i],.invhinge(x,m1$threshold))
  }
  o <- b / o
  w1 <- quantile(o,0.99)
  w2 <- quantile(o,0.01)
  w1 <- which(o > w1)
  w2 <- which(o < w2)
  if (length(w1) > 0) o <- o[-w1]
  if (length(w2) > 0) o <- o[-w2]
  o
}
#----------


.getSingleFeatureFrame <- function(n,cls,params=NULL,response=NULL) {
  o <- new('.featureFrame',var=n)
  if (cls == '.var') {
    o@feature.name <- n
    o@type <- 'linear'
  } else if (cls == '.quad') {
    o@feature.name <- paste0(n,'.quad')
    o@type <- 'quad'
  } else if (cls == '.cubic') {
    o@feature.name <- paste0(n,'.cubic')
    o@type <- 'cubic'
  } else if (cls == '.factor') {
    o@feature.name <- n
    o@type <- 'factor'
  } else if (cls == '.product') {
    o@feature.name <- paste(c('product_',n),collapse='.')
    o@type <- 'product'
  } else if (cls == '.log') {
    o@feature.name <- paste0(params[[1]],'.',n)
    o@type <- 'log'
    o@params <- params
  } else if (cls == '.poly') {
    d <- params$degree
    if (is.null(d)) d <- 3
    o@feature.name <- paste0(paste0(n,'.poly'),1:d)
    o@type <- 'poly'
    o@params <- params
  } else if (cls == '.hinge') {
    o@feature.name <- paste0(n,'.hinge')
    o@type <- 'hinge'
    o@params <- params
    o@response <- response
  } else if (cls == '.threshold') {
    o@feature.name <- paste0(n,'.threshold')
    o@type <- 'threshold'
    o@params <- params
    o@response <- response
  } else if (cls == '.func') {
    o@var <- character()
    o@feature.name <- deparse(params[[1]])
    o@type <- 'func'
    o@params <- params
  } else if (cls == '.simple.func') {
    o@var <- character()
    o@feature.name <- deparse(params[[1]])
    o@params <- params
    o@type <- 'simple.func'
  } else {
    o@feature.name <- paste0(cls,'.',n)
    o@type <- cls
    o@params <- params
    o@response <- response
  }
  o
}
#-----------

# need to be revised to support when there is multispecies with different indexes (records) for each species
.getFeaturetype <- function(d,f) {
  if (missing(f)) f <- d@sdmFormula
  if (inherits(f,'formula')) f <- .exFormula(f,as.data.frame(d))
  ff <- new('featuresFrame')
  o <- list()
  os <- om <- NULL
  for (m in f@model.terms) {
    fc <- class(m)
    if (fc %in% c('.var','.quad','.cubic','.factor')) o <- c(o,.getSingleFeatureFrame(slot(m,slotNames(m)[1]),fc))
    else if (fc == '.poly') o <- c(o,.getSingleFeatureFrame(m@x,fc,list(degree=m@degree,raw=m@raw)))
    else if (fc == '.product') o <- c(o,.getSingleFeatureFrame(m@x,fc))
    else if (fc == '.log') o <- c(o,.getSingleFeatureFrame(m@x,fc,list(as.character(m@term[[1]]))))
    else if (fc == '.func') {
      l <- as.list(m@x)
      if (as.character(l[[1]]) %in% c(':','*')) {
        o <- c(o,.getSingleFeatureFrame(as.character(.split.formula(m@x,l[[1]])),'.product'))
      } else if (as.character(l[[1]]) %in% c('>','<') && length(l) == 3) {
        o <- c(o,.getSingleFeatureFrame(as.character(l[[2]]),'.threshold',list(threshold=l[[3]],increasing=as.character(l[[1]]) == '>')))
      } else if (as.character(l[[1]]) %in% c('^') && length(l) == 3) {
        if (l[[3]] == 2) o <- c(o,.getSingleFeatureFrame(as.character(l[[2]]),'.quad'))
        else if (l[[3]] == 3) o <- c(o,.getSingleFeatureFrame(as.character(l[[2]]),'.cubic'))
        else o <- c(o,.getSingleFeatureFrame('xxx',fc,list(m@term)))
      } else {
        o <- c(o,.getSingleFeatureFrame('xxx',fc,list(m@term)))
      }
      
    } else if (fc == '.simple.func') {
      o <- c(o,.getSingleFeatureFrame('xxx',fc,list(m@term)))
    } else if (fc == '.hinge') {
      if (length(m@threshold) > 0 & length(m@increasing) > 0) o <- c(o,.getSingleFeatureFrame(m@x,fc,list(threshold=m@threshold,increasing=m@increasing)))
      else {
        if (length(f@species) == 0)  stop('The parameters for the hinge feature cannot be estimated without response variable!')
        else {
          if (!is.null(os)) os <- list()
          if (!all(f@species %in% d@species.names)) stop('the specified species in the formula has NO record in the data!')
          for (s in f@species) {
            ind <- .getIndex(d,sp=s)
            p <- .getHingeParams(x=.getFeatureRecords(d,ind=ind,n=m@x)[,2],y=.getSpeciesRecords(d,sp=s,ind=ind)[,2])
            os <- c(os,.getSingleFeatureFrame(m@x,fc,p,s))
          }
        }
      }
      
    } else if (fc == '.threshold') {
      if (length(m@threshold) > 0 & length(m@increasing) > 0) o <- c(o,.getSingleFeatureFrame(m@x,fc,list(threshold=m@threshold,increasing=m@increasing)))
      else {
        if (length(f@species) == 0) stop('The parameters for the threshold feature cannot be estimated without response variable!')
        else {
          if (!f@species %in% d@species.names) stop('the specified species in the formula has NO record in the data!')
          if (!is.null(os)) os <- list()
          for (s in f@species) {
            ind <- .getIndex(d,sp=s)
            p <- .getThresholdParams(x=.getFeatureRecords(d,ind=ind,n=m@x)[,2],y=.getSpeciesRecords(d,sp=s,ind=ind)[,2])
            os <- c(os,.getSingleFeatureFrame(m@x,fc,p,s))
          }
        }
      }
      
    } else if (fc == '.selectFrame') {
      # .getFeaturetype(d,f) # recursive
      # need to be checked and tested to select the final set!
    } else if (fc == '.nestedModel') {
      # need to be completed
    }
  }
  ff@vars <- unique(unlist(lapply(o,function(x) x@var)))
  ff@feature.types <- o
  if (!is.null(os)) {
    ff@vars <- unique(c(ff@vars,unlist(lapply(os,function(x) x@var))))
    ff@response.specific <- os
  }
  
  if (!is.null(om)) {
    ff@vars <- unique(c(ff@vars,unlist(lapply(om,function(x) x@var))))
    ff@model.specific <- om
  }
  
  ff
}




.getModelFrame <- function(x,data,response=NULL,dummy=FALSE) {
  o <- NULL
  if (all(x@vars %in% colnames(data))) {
    d <- data.frame(matrix(nrow=nrow(data),ncol=0))
    
    n <- lapply(x@feature.types,function(x) x@var)
    fn <- lapply(x@feature.types,function(x) x@feature.name)
    ft <- lapply(x@feature.types,function(x) x@type)
    
    un <- unique(n)
    u <- list()
    for (i in un) u[[i]] <- list()
    for (i in 1:length(n)) {
      u[[n[[i]]]][[length(u[[n[[i]]]])+1]] <-list(ft=ft[[i]],fn=fn[[i]],nr=i)
    }
    nrm <- c()
    for (i in 1:length(u)) {
      ft2 <- lapply(u[[i]],function(x) x$ft)
      fn2 <-lapply(u[[i]],function(x) x$fn)
      nr <- sapply(u[[i]],function(x) x$nr)
      
      if ('poly' %in% ft2) {
        w <- which(ft2 == 'poly')
        dg <- x@feature.types[[nr[w]]]@params$degree
        if (dg == 1) {
          if ('linear' %in% ft2) {
            nrm <- c(nrm,nr[w])
          } else {
            fn[[nr[w]]] <- n[[nr[w]]]
            ft[[nr[w]]] <- 'linear'
          }
        } else if (dg > 1) {
          if ('quad' %in% ft2) {
            w <- which(ft2 == 'quad')
            nrm <- c(nrm,nr[w])
          }
          
          if ('linear' %in% ft2) {
            w <- which(ft2 == 'linear')
            nrm <- c(nrm,nr[w])
          }
          if (dg > 2) {
            if ('cubic' %in% ft2) {
              w <- which(ft2 == 'cubic')
              nrm <- c(nrm,nr[w])
            }
          }
        }
      }
    }
    if (length(nrm) > 0) {
      ft <- ft[-nrm]
      fn <- fn[-nrm]
      n <- n[-nrm]
      x@feature.types <- x@feature.types[-nrm]
    }
    #---
    
    for (i in 1:length(n)) {
      if (ft[[i]] == 'linear') d[[fn[[i]]]] <- data[,n[[i]]]
      else if (ft[[i]] == 'factor') d[[fn[[i]]]] <- factor(data[,n[[i]]])
      else if (ft[[i]] == 'quad') d[[fn[[i]]]] <- .getFeature.quad(data[,n[[i]]])
      else if (ft[[i]] == 'cubic') d[[fn[[i]]]] <- .getFeature.cubic(data[,n[[i]]])
      else if (ft[[i]] == 'poly') {
        temp <- do.call(".getFeature.poly",c(list(x=data[,n[[i]]]),x@feature.types[[i]]@params))
        colnames(temp) <- x@feature.types[[i]]@feature.name
        d <- cbind(d,temp)
      }
      else if (ft[[i]] == 'product') d[[fn[[i]]]] <- .getFeature.product(data[,x@feature.types[[i]]@var])
      else if (ft[[i]] == 'log') d[[fn[[i]]]] <- do.call(x@feature.types[[i]]@params[[1]],list(x=data[,n[[i]]]))
      else if (ft[[i]] %in% c('func','simple.func')) {
        temp <- model.frame(as.formula(paste('~',deparse(x@feature.types[[i]]@params[[1]]))),data=data)
        n <- colnames(temp)
        d[[n]] <- as.vector(temp[,1])
      } else if (ft[[i]] == 'threshold') d[[fn[[i]]]] <- .getFeature.threshold(data[,n[[i]]],th = x@feature.types[[i]]@params$threshold,increasing=x@feature.types[[i]]@params$increasing)
      else if (ft[[i]] == 'hinge') d[[fn[[i]]]] <- .getFeature.hinge(data[,n[[i]]],th = x@feature.types[[i]]@params$threshold,increasing=x@feature.types[[i]]@params$increasing)
    }
    if (!is.null(x@response.specific)) {
      sp <- lapply(x@response.specific, function(x) x@response)
      if (!is.null(response)) {
        w <- response %in% sp
        if (any(w)) sp <- response[w]
        else warning('the specified response variable(s) do not exist in the data, all the existing responses are considered!')
      }
      o <- list()
      for (s in sp) {
        dd <- data.frame(matrix(nrow=nrow(data),ncol=0))
        n <- lapply(x@response.specific,function(x) x@var)
        fn <- lapply(x@response.specific,function(x) x@feature.name)
        ft <- lapply(x@response.specific,function(x) x@type)
        for (i in 1:length(n)) {
          if (ft[[i]] == 'threshold') dd[[fn[[i]]]] <- .getFeature.threshold(data[,n[[i]]],th = x@response.specific[[i]]@params$threshold,increasing=x@response.specific[[i]]@params$increasing)
          else if (ft[[i]] == 'hinge') dd[[fn[[i]]]] <- .getFeature.hinge(data[,n[[i]]],th = x@response.specific[[i]]@params$threshold,increasing=x@response.specific[[i]]@params$increasing)
        }
        o[[s]] <- dd
      }
    }
  } else stop('some of the specified variables in the formula do not exist in the data!')
  list(features=d,specis_specific=o) 
}
#-----------


.getFeatureNamesTypes <- function(x,merged=TRUE) {
  # get the name of features from a featureFrame onject (x)
  # If merged is FALSE, separately reported in a list
  n1 <- unlist(lapply(x@feature.types,function(x) x@feature.name))
  n2 <- unlist(lapply(x@response.specific ,function(x) x@feature.name))
  
  f1 <- unlist(lapply(x@feature.types,function(x) x@type))
  f2 <- unlist(lapply(x@response.specific,function(x) x@type))
  names(f1) <- n1
  names(f2) <- n2
  if (merged) c(f1,f2)
  else list(features=f1,response.specific=f2)
}
