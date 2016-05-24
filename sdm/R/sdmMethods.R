# Author: Babak Naimi, naimi.b@gmail.com
# Date :  March 2016
# Version 2.0
# Licence GPL v3


if (!isGeneric("add")) {
  setGeneric("add", function(x,w,echo,...)
    standardGeneric("add"))
}


setMethod('add', signature(x='list',w='character'), 
          function(x,w='sdm',echo=TRUE,...) {
            #dot <- list(...)
            if (missing(w)) w <- 'sdm'
            if (missing(echo) || !is.logical(echo)) echo <- TRUE
            
            if (w %in% c('sdm','sdmMethod','model','sdmCorrelativeMethod')) {
              m <- do.call('.create.sdmCorrelativeMethod',x)
              if (inherits(m,"sdmCorrelativeMethod")) {
                .sdmMethods$addMethod(m,echo)
              } else stop('The method is not added to the sdmMethods contrainer!')
            }
          }
)



if (!isGeneric("getmethod")) {
  setGeneric("getmethod", function(x,w,...)
    standardGeneric("getmethod"))
}


setMethod('getmethod', signature(x='character'), 
          function(x,w,...) {
            #dot <- list(...)
            if (missing(w)) w <- 'sdm'
            if (w == 'sdm') {
              x <- .methodFix(x)
              if (!is.na(x)) .sdmMethods$Methods[[x]]
              else stop('the specified method does not exist!')
            }
            #n <- getmethodNames('sdm',FALSE)
          }
)

if (!isGeneric("getmethodNames")) {
  setGeneric("getmethodNames", function(w,...)
    standardGeneric("getmethodNames"))
}


setMethod('getmethodNames', signature(w='ANY'), 
          function(w,alt,...) {
            if (missing(w)) w <- 'sdm'
            if (missing(alt)) alt <- TRUE
            if (w == 'sdm') .sdmMethods$getMethodNames(alt=alt)
          }
)


if (!isGeneric("getModelInfo")) {
  setGeneric("getModelInfo", function(x,...)
    standardGeneric("getModelInfo"))
}


setMethod('getModelInfo', signature(x='sdmModels'), 
          function(x,w,...) {
            if (missing(w)) w <- NULL
            .getModel.info(x,w,...)
            
          }
)



.addMethods <- function() {
  methodInfo <- NULL
  n <- getmethodNames('sdm',alt=FALSE)
  lst <- list.files(system.file("methods/sdm", package="sdm"),pattern='R$',full.names = TRUE)
  for (l in lst) {
    source(l,local=TRUE)
    pkg <- methodInfo$packages
    pkg <- pkg[!pkg == '.tmp']
    if (!methodInfo$name[1] %in% n && all(.is.installed(pkg))) {
      add(x=methodInfo,'sdm',echo=FALSE)
    }
  }
  .sdmOptions$addOption('sdmLoaded',TRUE)
}

.is.installed <- function(n) {
  inst <- utils::installed.packages()[,1]
  nn <- n %in% inst
  names(nn) <- n
  nn
}


.create.sdmCorrelativeMethod <- function(name,packages=NULL,modelTypes=NULL,fitParams,fitSettings=NULL,settingRules=NULL,fitFunction,predictParams=NULL,predictSettings=NULL,predictFunction=NULL,tuneParams=NULL,metadata=NULL,...) {
  m <- new('sdmCorrelativeMethod',name=name[1])
  if (length(name) > 1) m@aliases <- name[2:length(name)]
  
  Installed <- TRUE
  
  if (!is.null(packages) && is.character(packages)) {
    m@packages <- packages
    #w <- rep(FALSE,length(m@packages))
    
    w <- .is.installed(m@packages)
    
    #if (!all(w)) print(paste('warning: packages (',paste(m@packages[!w],collapse=', '),') need to be installed to get this method working!',sep=''))
    if (!all(w)) Installed <- FALSE
    else {
      for (i in seq_along(m@packages)) w[i] <- require(m@packages[i],character.only=TRUE)
    }
  } else m@packages <- NULL
  
  if (!is.null(modelTypes)) {
    modelTypes <- tolower(modelTypes)
    for (i in 1:length(modelTypes)) {
      if (modelTypes[i] %in% c('po','presenceonly','presence-only','presence')) {
        modelTypes[i] <- 'po'
      } else if (modelTypes[i] %in% c('pa','presenceabsence','presence-absence')) {
        modelTypes[i] <- 'pa'
      } else if (modelTypes[i] %in% c('pb','presenceb','presence-background','presence-pseudo','presence-pseudoabsence','ppa','psa')) {
        modelTypes[i] <- 'pb'
      } else if (modelTypes[i] %in% c('ab','abundance')) {
        modelTypes[i] <- 'ab'
      } else if (modelTypes[i] %in% c('n','nominal','multinominal')) {
        modelTypes[i] <- 'n'
      } else {
        warning(paste('modelType',modelTypes[i],'is unknown, it is ignored!'))
        modelTypes[i] <- NA
      }
    }
    m@modelTypes <- modelTypes
  }
  
  #-------
  if (is.list(fitParams)) {
    n <- names(fitParams)
    if (is.null(n)) stop('fitParams is not appropriately defined; example: list(formula="standard.formula",data="sdmDataFrame")')
    m@fitParams <- fitParams
  } else stop('fitParams should be a list')
  #------
  if (!is.null(fitSettings)) {
    if (!is.list(fitSettings)) stop('fitSettings should be a list!')
    n <- names(fitSettings)
    if (is.null(n)) stop('fitSettings is not appropriately defined; example: list(family=link(binomial),ntrees=1000)')
    if ('' %in% n) {
      w <- which(n == '')
      for (ww in w) {
        if (is.character(n[ww])) names(fitSettings)[ww] <- n[w]
        else stop('fitSettings is not appropriately defined; example: list(family=link(binomial),ntrees=1000)')
      }
    }
    m@fitSettings <- fitSettings
  }
  #------
  if (!is.null(settingRules)) {
    if (!is.function(settingRules)) stop('settingRules should be a function!')
    m@settingRules <- settingRules
  }
  #-----
  
  if (class(fitFunction) == "character") {
    if (length(strsplit(fitFunction,'::')[[1]]) == 2) fitFunction <- strsplit(fitFunction,'::')[[1]][2]
    
    if (exists(fitFunction,mode='function')) {
      
      if (environmentName(environment(get(fitFunction))) == "R_GlobalEnv") {
        # assign to the environment in the container of methods!
        if (is.null(m@.temp.env)) m@.temp.env <- new.env()
        assign(fitFunction,get(fitFunction),envir = m@.temp.env) ####
        m@packages <- unique(c(m@packages,'.temp'))
        # when environment if attached, the conflict with the existing object in .GlobalEnv should be resolved!
      }
      m@fitFunction <- eval(parse(text=paste("function(params) {do.call(",fitFunction,",params)}")))
    } else if (length(strsplit(fitFunction,':::')[[1]]) == 2 && class(eval(parse(text=fitFunction))) == 'function') {
      if (!exists(strsplit(fitFunction,':::')[[1]][2],mode='function')) {
        m@fitFunction <- eval(parse(text=paste("function(params) {do.call(",strsplit(fitFunction,':::')[[1]][2],",params)}")))
      } else {
        if (is.null(m@.temp.env)) m@.temp.env <- new.env()
        assign(strsplit(fitFunction,':::')[[1]][2],eval(parse(text=fitFunction)),envir = m@.temp.env) ####
        m@packages <- unique(c(m@packages,'.temp'))
        m@fitFunction <- eval(parse(text=paste("function(params) {do.call(",strsplit(fitFunction,':::')[[1]][2],",params)}")))
      }
      
    } else if (!Installed) {
      m@fitFunction <- eval(parse(text=paste("function(params) {do.call(",fitFunction,",params)}")))
    } else stop('fitFunction cannot be identified!')
    
  } else if (class(fitFunction) == 'function') {
    if (is.null(m@.temp.env)) m@.temp.env <- new.env()
    if (!paste(m@name,'.fit',sep='') %in% ls(envir=m@.temp.env)) fn <-paste(m@name,'.fit',sep='') ####
    else stop('the user defined function in fitFunction cannot be registered because an object with a similar name exists in the container!')
    assign(fn,fitFunction,envir = m@.temp.env) ####
    m@packages <- unique(c(m@packages,'.temp'))
    m@fitFunction <- eval(parse(text=paste("function(params) {do.call(",fn,",params)}")))
  } else stop('fitFunction cannot be identified!')
  #---------
  
  if (!is.null(predictFunction)) {
    if (class(predictFunction) == "character") {
      if (length(strsplit(predictFunction,'::')[[1]]) == 2) predictFunction <- strsplit(predictFunction,'::')[[1]][2]
      if (exists(predictFunction,mode='function')) {
        if (environmentName(environment(get(predictFunction))) == "R_GlobalEnv") {
          # assign to the environment in the container of methods!
          if (is.null(m@.temp.env)) m@.temp.env <- new.env()
          assign(predictFunction,get(predictFunction),envir = m@.temp.env) ####
          m@packages <- unique(c(m@packages,'.temp'))
          # when environment if attached, the conflict with the existing object in .GlobalEnv should be resolved!
        }
        m@predictFunction <- eval(parse(text=paste("function(params) {do.call(",predictFunction,",params)}")))
      } else if (length(strsplit(predictFunction,':::')[[1]]) == 2 && class(eval(parse(text=predictFunction))) == 'function') {
        if (!exists(strsplit(predictFunction,':::')[[1]][2],mode='function')) {
          m@predictFunction <- eval(parse(text=paste("function(params) {do.call(",strsplit(predictFunction,':::')[[1]][2],",params)}")))
        } else {
          if (is.null(m@.temp.env)) m@.temp.env <- new.env()
          assign(strsplit(predictFunction,':::')[[1]][2],eval(parse(text=predictFunction)),envir = m@.temp.env) ####
          m@packages <- unique(c(m@packages,'.temp'))
          m@predictFunction <- eval(parse(text=paste("function(params) {do.call(",strsplit(predictFunction,':::')[[1]][2],",params)}")))
        }
        
      } else if (!Installed) {
        m@predictFunction <- eval(parse(text=paste("function(params) {do.call(",predictFunction,",params)}")))
      } else stop('predictFunction cannot be identified!')
    } else if (class(predictFunction) == 'function') {
      if (is.null(m@.temp.env)) m@.temp.env <- new.env()
      if (!paste(m@name,'.predict',sep='') %in% ls(envir=m@.temp.env)) fn <-paste(m@name,'.predict',sep='') ####
      else stop('the user defined function in predictFunction cannot be registered because an object with a similar name exists in the container!')
      assign(fn,predictFunction,envir = m@.temp.env) ####
      m@packages <- unique(c(m@packages,'.temp'))
      m@predictFunction <- eval(parse(text=paste("function(params) {do.call(",fn,",params)}")))
    } else stop('predictFunction cannot be identified!')
    #---------
    if (!is.null(predictParams)) {
      if (is.list(predictParams)) {
        n <- names(predictParams)
        if (is.null(n)) stop('predictParams is not appropriately defined; example: list(newdata="sdmDataFrame")')
        m@predictParams <- predictParams
      } else stop('predictParams should be a list')
    }
    if (!is.null(predictSettings)) {
      if (is.list(predictSettings)) m@predictSettings <- predictSettings
      else stop('predictSettings should be a list!')
    }
  }
  #-------- 
  if (!is.null(tuneParams)) {
    if (!is.list(tuneParams)) stop('tuneParams should be a list; example: list(ntrees=seq(500,3000,by=200))')
    n <- names(tuneParams)
    if (is.null(n)) stop('tuneParams is not appropriately defined; example: list(ntrees=seq(500,3000,by=200))')
    m@tuneParams <- tuneParams
  } 
  #------------  
  if (inherits(metadata,'.Metadata')) m@metadata <- metadata
  else m@metadata <- .newMetadata(...)
  
  m
}
# 
# 
# .create.sdmCorrelativeMethod <- function(name,packages=NULL,modelTypes=NULL,fitParams,fitSettings=NULL,settingRules=NULL,fitFunction,predictParams=NULL,predictSettings=NULL,predictFunction=NULL,tuneParams=NULL,metadata=NULL,...) {
#   m <- new('sdmCorrelativeMethod',name=name[1])
#   if (length(name) > 1) m@aliases <- name[2:length(name)]
#   if (!is.null(packages) && is.character(packages)) {
#     m@packages <- packages
#     w <- rep(FALSE,length(m@packages))
#     for (i in seq_along(m@packages)) w[i] <- require(m@packages[i],character.only=TRUE)
#     if (!all(w)) print(paste('warning: packages (',paste(m@packages[!w],collapse=', '),') need to be installed to get this method working!',sep=''))
#   } else m@packages <- NULL
#   
#   if (!is.null(modelTypes)) {
#     modelTypes <- tolower(modelTypes)
#     for (i in 1:length(modelTypes)) {
#       if (modelTypes[i] %in% c('po','presenceonly','presence-only','presence')) {
#         modelTypes[i] <- 'po'
#       } else if (modelTypes[i] %in% c('pa','presenceabsence','presence-absence')) {
#         modelTypes[i] <- 'pa'
#       } else if (modelTypes[i] %in% c('pb','presenceb','presence-background','presence-pseudo','presence-pseudoabsence','ppa','psa')) {
#         modelTypes[i] <- 'pb'
#       } else if (modelTypes[i] %in% c('ab','abundance')) {
#         modelTypes[i] <- 'ab'
#       } else if (modelTypes[i] %in% c('n','nominal','multinominal')) {
#         modelTypes[i] <- 'n'
#       } else {
#         warning(paste('modelType',modelTypes[i],'is unknown, it is ignored!'))
#         modelTypes[i] <- NA
#       }
#     }
#     m@modelTypes <- modelTypes
#   }
#     
#   #-------
#   if (is.list(fitParams)) {
#     n <- names(fitParams)
#     if (is.null(n)) stop('fitParams is not appropriately defined; example: list(formula="standard.formula",data="sdmDataFrame")')
#     m@fitParams <- fitParams
#   } else stop('fitParams should be a list')
#   #------
#   if (!is.null(fitSettings)) {
#     if (!is.list(fitSettings)) stop('fitSettings should be a list!')
#     n <- names(fitSettings)
#     if (is.null(n)) stop('fitSettings is not appropriately defined; example: list(family=link(binomial),ntrees=1000)')
#     if ('' %in% n) {
#       w <- which(n == '')
#       for (ww in w) {
#         if (is.character(n[ww])) names(fitSettings)[ww] <- n[w]
#         else stop('fitSettings is not appropriately defined; example: list(family=link(binomial),ntrees=1000)')
#       }
#     }
#     m@fitSettings <- fitSettings
#   }
#   #------
#   if (!is.null(settingRules)) {
#     if (!is.function(settingRules)) stop('settingRules should be a function!')
#     m@settingRules <- settingRules
#   }
#   #-----
#   if (exists(as.character(substitute(fitFunction)),mode='function') || exists(as.character(substitute(fitFunction)), envir = asNamespace('sdm'), inherits = FALSE)) {
#     if (class(substitute(fitFunction)) == 'name') {
#       fn <- as.character(substitute(fitFunction))
#       if (environmentName(environment(fitFunction)) == "R_GlobalEnv") {
#         # assign to the environment in the container of methods!
#         if (is.null(m@.temp.env)) m@.temp.env <- new.env()
#         assign(fn,fitFunction,envir = m@.temp.env) ####
#         m@packages <- unique(c(m@packages,'.temp'))
#         # when environment if attached, the conflict with the existing object in .GlobalEnv should be resolved!
#       } 
#       m@fitFunction <- eval(parse(text=paste("function(params) {do.call(",fn,",params)}")))
#       
#     } else if (class(substitute(fitFunction)) == 'call') {
#       if (environmentName(environment(fitFunction)) == "R_GlobalEnv" || environmentName(environment(fitFunction)) == "sdm") {
#         if (is.null(m@.temp.env)) m@.temp.env <- new.env()
#         if (!paste(m@name,'.fit',sep='') %in% ls(envir=m@.temp.env)) fn <-paste(m@name,'.fit',sep='') ####
#         else stop('the user defined function in fitFunction cannot be registered because an object with a similar name exists in the container!')
#         
#         assign(fn,fitFunction,envir = m@.temp.env) ####
#         m@packages <- unique(c(m@packages,'.temp'))
#       } else {
#         if (as.character(substitute(fitFunction))[1] == "::") fn <- as.character(substitute(fitFunction))[3]
#         else stop('fitFunction cannot be identified!')
#       } 
#       m@fitFunction <- eval(parse(text=paste("function(params) {do.call(",fn,",params)}")))
#     } else stop('fitFunction cannot be identified!')
#     
#   } else if (exists(as.character(substitute(fitFunction)), envir = asNamespace('sdm'), inherits = FALSE)) {
#     fn <- as.character(substitute(fitFunction))
#     m@fitFunction <- eval(parse(text=paste("function(params) {do.call(",fn,",params)}")))
#   } else stop(paste("Function",as.character(substitute(fitFunction)),"is not found!"))
#   #---------
#   
#   if (!is.null(predictFunction)) {
#     if (exists(as.character(substitute(predictFunction)),mode='function') || exists(as.character(substitute(predictFunction)), envir = asNamespace('sdm'), inherits = FALSE)) {
#       
#       if (class(substitute(predictFunction)) == 'name') {
#         fn <- as.character(substitute(predictFunction))
#         if (environmentName(environment(predictFunction)) == "R_GlobalEnv" ) {
#           # assign to the environment in the container of methods!
#           if (is.null(m@.temp.env)) m@.temp.env <- new.env()
#           assign(fn,predictFunction,envir = m@.temp.env) ####
#           m@packages <- unique(c(m@packages,'.temp'))
#           # when environment if attached, the conflict with the existing object in .GlobalEnv should be resolved!
#         } 
#         m@predictFunction <- eval(parse(text=paste("function(params) {do.call(",fn,",params)}")))
#         
#       } else if (class(substitute(predictFunction)) == 'call') {
#         if (environmentName(environment(predictFunction)) == "R_GlobalEnv" || environmentName(environment(predictFunction)) == "sdm") {
#           if (is.null(m@.temp.env)) m@.temp.env <- new.env()
#           if (!paste(m@name,'.predict',sep='') %in% ls(envir=m@.temp.env)) fn <-paste(m@name,'.predict',sep='') ####
#           else stop('the user defined function in predictFunction cannot be registered because an object with a similar name exists in the container!')
#           assign(fn,predictFunction,envir = m@.temp.env) ####
#           m@packages <- unique(c(m@packages,'.temp'))
#         } else {
#           if (as.character(substitute(predictFunction))[1] == "::") fn <- as.character(substitute(predictFunction))[3]
#           else stop('predictFunction cannot be identified!')
#         } 
#         
#         m@predictFunction <- eval(parse(text=paste("function(params) {do.call(",fn,",params)}")))
#       } else if (exists(as.character(substitute(predictFunction)), envir = asNamespace('sdm'), inherits = FALSE)) {
#         fn <- as.character(substitute(predictFunction))
#         m@predictFunction <- eval(parse(text=paste("function(params) {do.call(",fn,",params)}")))
#       } else stop('predictFunction cannot be identified!')
#       
#     } else stop(paste("Function",as.character(substitute(predictFunction)),"is not found!"))
#   
#     if (!is.null(predictParams)) {
#       if (is.list(predictParams)) {
#         n <- names(predictParams)
#         if (is.null(n)) stop('predictParams is not appropriately defined; example: list(newdata="sdmDataFrame")')
#         m@predictParams <- predictParams
#       } else stop('predictParams should be a list')
#     }
#     if (!is.null(predictSettings)) {
#       if (is.list(predictSettings)) m@predictSettings <- predictSettings
#       else stop('predictSettings should be a list!')
#     }
#   }
#  #-------- 
#   if (!is.null(tuneParams)) {
#     if (!is.list(tuneParams)) stop('tuneParams should be a list; example: list(ntrees=seq(500,3000,by=200))')
#     n <- names(tuneParams)
#     if (is.null(n)) stop('tuneParams is not appropriately defined; example: list(ntrees=seq(500,3000,by=200))')
#     m@tuneParams <- tuneParams
#   } 
#   #------------  
#   if (inherits(metadata,'.Metadata')) m@metadata <- metadata
#   else m@metadata <- .newMetadata(...)
#   
#   m
# }

#-----

.update.sdmCorrelativeMethod <- function(m,...) {
  name=NULL;packages=NULL;modelTypes=NULL;fitParams=NULL;fitSettings=NULL;settingRules=NULL;fitFunction=NULL;predictParams=NULL;predictSettings=NULL;predictFunction=NULL;tuneParams=NULL;metadata=NULL
  
  dot <- list(...)
  n <- tolower(names(dot))
  for (i in seq_along(n)) {
    if (any(!is.na(pmatch(c("nam"),n[i])))) name <- dot[[i]]
    else if (any(!is.na(pmatch(c("pac"),n[i])))) packages <- dot[[i]]
    else if (any(!is.na(pmatch(c("mod"),n[i])))) modelTypes <- dot[[i]]
    else if (any(!is.na(pmatch(c("fits"),n[i])))) fitSettings <- dot[[i]]
    else if (any(!is.na(pmatch(c("fitp"),n[i])))) fitParams <- dot[[i]]
    else if (any(!is.na(pmatch(c("set"),n[i])))) settingRules <- dot[[i]]
    else if (any(!is.na(pmatch(c("fitf"),n[i])))) fitFunction <- dot[[i]]
    else if (any(!is.na(pmatch(c("predicts"),n[i])))) predictSettings <- dot[[i]]
    else if (any(!is.na(pmatch(c("predictp"),n[i])))) predictParams <- dot[[i]]
    else if (any(!is.na(pmatch(c("predictf"),n[i])))) predictFunction <- dot[[i]]
    else if (any(!is.na(pmatch(c("tun"),n[i])))) tuneParams <- dot[[i]]
    else if (any(!is.na(pmatch(c("met"),n[i])))) metadata <- dot[[i]]
  }
  #--------
  if (length(name) > 1) m@aliases <- name[2:length(name)]
  
  Installed <- TRUE
  
  if (!is.null(packages) && is.character(packages)) {
    m@packages <- packages
    #w <- rep(FALSE,length(m@packages))
    #for (i in seq_along(m@packages)) w[i] <- require(m@packages[i],character.only=TRUE)
    w <- .is.installed(m@packages)
    #if (!all(w)) print(paste('warning: packages (',paste(m@packages[!w],collapse=', '),') need to be installed to get this method working!',sep=''))
    if (!all(w)) Installed <- FALSE
  } else m@packages <- NULL
  
  if (!is.null(modelTypes)) {
    modelTypes <- tolower(modelTypes)
    for (i in 1:length(modelTypes)) {
      if (modelTypes[i] %in% c('po','presenceonly','presence-only','presence')) {
        modelTypes[i] <- 'po'
      } else if (modelTypes[i] %in% c('pa','presenceabsence','presence-absence')) {
        modelTypes[i] <- 'pa'
      } else if (modelTypes[i] %in% c('pb','presenceb','presence-background','presence-pseudo','presence-pseudoabsence','ppa','psa')) {
        modelTypes[i] <- 'pb'
      } else if (modelTypes[i] %in% c('ab','abundance')) {
        modelTypes[i] <- 'ab'
      } else if (modelTypes[i] %in% c('n','nominal','multinominal')) {
        modelTypes[i] <- 'n'
      } else {
        warning(paste('modelType',modelTypes[i],'is unknown, it is ignored!'))
        modelTypes[i] <- NA
      }
    }
    m@modelTypes <- modelTypes
  }
  
  #-------
  if (is.list(fitParams)) {
    n <- names(fitParams)
    if (is.null(n)) stop('fitParams is not appropriately defined; example: list(formula="standard.formula",data="sdmDataFrame")')
    m@fitParams <- fitParams
  } else stop('fitParams should be a list')
  #------
  if (!is.null(fitSettings)) {
    if (!is.list(fitSettings)) stop('fitSettings should be a list!')
    n <- names(fitSettings)
    if (is.null(n)) stop('fitSettings is not appropriately defined; example: list(family=link(binomial),ntrees=1000)')
    if ('' %in% n) {
      w <- which(n == '')
      for (ww in w) {
        if (is.character(n[ww])) names(fitSettings)[ww] <- n[w]
        else stop('fitSettings is not appropriately defined; example: list(family=link(binomial),ntrees=1000)')
      }
    }
    m@fitSettings <- fitSettings
  }
  #------
  if (!is.null(settingRules)) {
    if (!is.function(settingRules)) stop('settingRules should be a function!')
    m@settingRules <- settingRules
  }
  #-----
  
  if (class(fitFunction) == "character") {
    if (length(strsplit(fitFunction,'::')[[1]]) == 2) fitFunction <- strsplit(fitFunction,'::')[[1]][2]
    
    if (exists(fitFunction,mode='function')) {
      
      if (environmentName(environment(get(fitFunction))) == "R_GlobalEnv") {
        # assign to the environment in the container of methods!
        if (is.null(m@.temp.env)) m@.temp.env <- new.env()
        assign(fitFunction,get(fitFunction),envir = m@.temp.env) ####
        m@packages <- unique(c(m@packages,'.temp'))
        # when environment if attached, the conflict with the existing object in .GlobalEnv should be resolved!
      }
      m@fitFunction <- eval(parse(text=paste("function(params) {do.call(",fitFunction,",params)}")))
    } else if (length(strsplit(fitFunction,':::')[[1]]) == 2 && class(eval(parse(text=fitFunction))) == 'function') {
      if (!exists(strsplit(fitFunction,':::')[[1]][2],mode='function')) {
        m@fitFunction <- eval(parse(text=paste("function(params) {do.call(",strsplit(fitFunction,':::')[[1]][2],",params)}")))
      } else {
        if (is.null(m@.temp.env)) m@.temp.env <- new.env()
        assign(strsplit(fitFunction,':::')[[1]][2],eval(parse(text=fitFunction)),envir = m@.temp.env) ####
        m@packages <- unique(c(m@packages,'.temp'))
        m@fitFunction <- eval(parse(text=paste("function(params) {do.call(",strsplit(fitFunction,':::')[[1]][2],",params)}")))
      }
      
    } else if (!Installed) {
      m@fitFunction <- eval(parse(text=paste("function(params) {do.call(",fitFunction,",params)}")))
    } else stop('fitFunction cannot be identified!')
    
  } else if (class(fitFunction) == 'function') {
    if (is.null(m@.temp.env)) m@.temp.env <- new.env()
    if (!paste(m@name,'.fit',sep='') %in% ls(envir=m@.temp.env)) fn <-paste(m@name,'.fit',sep='') ####
    else stop('the user defined function in fitFunction cannot be registered because an object with a similar name exists in the container!')
    assign(fn,fitFunction,envir = m@.temp.env) ####
    m@packages <- unique(c(m@packages,'.temp'))
    m@fitFunction <- eval(parse(text=paste("function(params) {do.call(",fn,",params)}")))
  } else stop('fitFunction cannot be identified!')
  #---------
  
  if (!is.null(predictFunction)) {
    if (class(predictFunction) == "character") {
      if (length(strsplit(predictFunction,'::')[[1]]) == 2) predictFunction <- strsplit(predictFunction,'::')[[1]][2]
      if (exists(predictFunction,mode='function')) {
        if (environmentName(environment(get(predictFunction))) == "R_GlobalEnv") {
          # assign to the environment in the container of methods!
          if (is.null(m@.temp.env)) m@.temp.env <- new.env()
          assign(predictFunction,get(predictFunction),envir = m@.temp.env) ####
          m@packages <- unique(c(m@packages,'.temp'))
          # when environment if attached, the conflict with the existing object in .GlobalEnv should be resolved!
        }
        m@predictFunction <- eval(parse(text=paste("function(params) {do.call(",predictFunction,",params)}")))
      } else if (length(strsplit(predictFunction,':::')[[1]]) == 2 && class(eval(parse(text=predictFunction))) == 'function') {
        if (!exists(strsplit(predictFunction,':::')[[1]][2],mode='function')) {
          m@predictFunction <- eval(parse(text=paste("function(params) {do.call(",strsplit(predictFunction,':::')[[1]][2],",params)}")))
        } else {
          if (is.null(m@.temp.env)) m@.temp.env <- new.env()
          assign(strsplit(predictFunction,':::')[[1]][2],eval(parse(text=predictFunction)),envir = m@.temp.env) ####
          m@packages <- unique(c(m@packages,'.temp'))
          m@predictFunction <- eval(parse(text=paste("function(params) {do.call(",strsplit(predictFunction,':::')[[1]][2],",params)}")))
        }
        
      } else if (!Installed) {
        m@predictFunction <- eval(parse(text=paste("function(params) {do.call(",predictFunction,",params)}")))
      } else stop('predictFunction cannot be identified!')
    } else if (class(predictFunction) == 'function') {
      if (is.null(m@.temp.env)) m@.temp.env <- new.env()
      if (!paste(m@name,'.predict',sep='') %in% ls(envir=m@.temp.env)) fn <-paste(m@name,'.predict',sep='') ####
      else stop('the user defined function in predictFunction cannot be registered because an object with a similar name exists in the container!')
      assign(fn,predictFunction,envir = m@.temp.env) ####
      m@packages <- unique(c(m@packages,'.temp'))
      m@predictFunction <- eval(parse(text=paste("function(params) {do.call(",fn,",params)}")))
    } else stop('predictFunction cannot be identified!')
    #---------
    if (!is.null(predictParams)) {
      if (is.list(predictParams)) {
        n <- names(predictParams)
        if (is.null(n)) stop('predictParams is not appropriately defined; example: list(newdata="sdmDataFrame")')
        m@predictParams <- predictParams
      } else stop('predictParams should be a list')
    }
    if (!is.null(predictSettings)) {
      if (is.list(predictSettings)) m@predictSettings <- predictSettings
      else stop('predictSettings should be a list!')
    }
  }
  #-------- 
  if (!is.null(tuneParams)) {
    if (!is.list(tuneParams)) stop('tuneParams should be a list; example: list(ntrees=seq(500,3000,by=200))')
    n <- names(tuneParams)
    if (is.null(n)) stop('tuneParams is not appropriately defined; example: list(ntrees=seq(500,3000,by=200))')
    m@tuneParams <- tuneParams
  } 
  #------------  
  if (inherits(metadata,'.Metadata')) m@metadata <- metadata
  else m@metadata <- .newMetadata(...)
  
  m
}
#########
# .update.sdmCorrelativeMethod <- function(m,...) {
#   name=NULL;packages=NULL;modelTypes=NULL;fitParams=NULL;fitSettings=NULL;settingRules=NULL;fitFunction=NULL;predictParams=NULL;predictSettings=NULL;predictFunction=NULL;tuneParams=NULL;metadata=NULL
#   
#   dot <- list(...)
#   n <- tolower(names(dot))
#   for (i in seq_along(n)) {
#     if (any(!is.na(pmatch(c("nam"),n[i])))) name <- dot[[i]]
#     else if (any(!is.na(pmatch(c("pac"),n[i])))) packages <- dot[[i]]
#     else if (any(!is.na(pmatch(c("mod"),n[i])))) modelTypes <- dot[[i]]
#     else if (any(!is.na(pmatch(c("fits"),n[i])))) fitSettings <- dot[[i]]
#     else if (any(!is.na(pmatch(c("fitp"),n[i])))) fitParams <- dot[[i]]
#     else if (any(!is.na(pmatch(c("set"),n[i])))) settingRules <- dot[[i]]
#     else if (any(!is.na(pmatch(c("fitf"),n[i])))) fitFunction <- dot[[i]]
#     else if (any(!is.na(pmatch(c("predicts"),n[i])))) predictSettings <- dot[[i]]
#     else if (any(!is.na(pmatch(c("predictp"),n[i])))) predictParams <- dot[[i]]
#     else if (any(!is.na(pmatch(c("predictf"),n[i])))) predictFunction <- dot[[i]]
#     else if (any(!is.na(pmatch(c("tun"),n[i])))) tuneParams <- dot[[i]]
#     else if (any(!is.na(pmatch(c("met"),n[i])))) metadata <- dot[[i]]
#   }
#   #--------
#   if (!is.null(name)) {
#     if (length(name) > 1) m@aliases <- name[2:length(name)]
#   }
#   
#   if (!is.null(packages) && is.character(packages)) {
#     m@packages <- packages
#     w <- rep(FALSE,length(m@packages))
#     for (i in seq_along(m@packages)) w[i] <- require(m@packages[i],character.only=TRUE)
#     if (!all(w)) print(paste('warning: packages (',paste(m@packages[!w],collapse=', '),') need to be installed to get this method working!',sep=''))
#   }
#   
#   if (!is.null(modelTypes)) {
#     modelTypes <- tolower(modelTypes)
#     for (i in seq_along(modelTypes)) {
#       if (modelTypes[i] %in% c('po','presenceonly','presence-only','presence')) modelTypes[i] <- 'po'
#       else if (modelTypes[i] %in% c('pa','presenceabsence','presence-absence')) modelTypes[i] <- 'pa'
#       else if (modelTypes[i] %in% c('pb','presenceb','presence-background','presence-pseudo','presence-pseudoabsence','ppa','psa')) modelTypes[i] <- 'pb'
#       else if (modelTypes[i] %in% c('ab','abundance')) modelTypes[i] <- 'ab'
#       else if (modelTypes[i] %in% c('n','nominal','multinominal')) modelTypes[i] <- 'n'
#       else {
#         warning(paste('modelType',modelTypes[i],'is unknown, it is ignored!'))
#         modelTypes[i] <- NULL
#       }
#     }
#     m@modelTypes <- modelTypes
#   }
#   
#   #-------
#   if (!is.null(fitParams) && is.list(fitParams)) {
#     n <- names(fitParams)
#     if (is.null(n)) stop('fitParams is not appropriately defined; example: list(formula="standard.formula",data="sdmDataFrame")')
#     m@fitParams <- fitParams
#   }
#   #------
#   if (!is.null(fitSettings) && is.list(fitSettings)) {
#     n <- names(fitSettings)
#     if (is.null(n)) stop('fitSettings is not appropriately defined; example: list(family=link(binomial),ntrees=1000)')
#     if ('' %in% n) {
#       w <- which(n == '')
#       for (ww in w) {
#         if (is.character(n[ww])) names(fitSettings)[ww] <- n[w]
#         else stop('fitSettings is not appropriately defined; example: list(family=link(binomial),ntrees=1000)')
#       }
#     }
#     m@fitSettings <- fitSettings
#   }
#   #------
#   if (!is.null(settingRules)) {
#     if (!is.function(settingRules)) stop('settingRules should be a function!')
#     m@settingRules <- settingRules
#   }
#   #-----
#   if (!is.null(fitFunction) && exists(as.character(substitute(fitFunction)),mode='function')) {
#     if (class(substitute(fitFunction)) == 'name') {
#       fn <- as.character(substitute(fitFunction))
#       if (environmentName(environment(fitFunction)) == "R_GlobalEnv") {
#         # assign to the environment in the container of methods!
#         if (is.null(m@.temp.env)) m@.temp.env <- new.env()
#         assign(fn,fitFunction,envir = m@.temp.env) ####
#         m@packages <- unique(c(m@packages,'.temp'))
#         # when environment if attached, the conflict with the existing object in .GlobalEnv should be resolved!
#       } 
#       m@fitFunction <- eval(parse(text=paste("function(params) {do.call(",fn,",params)}")))
#       
#     } else if (class(substitute(fitFunction)) == 'call') {
#       
#       if (environmentName(environment(fitFunction)) == "R_GlobalEnv" || environmentName(environment(fitFunction)) == "sdm") {
#         if (is.null(m@.temp.env)) m@.temp.env <- new.env()
#         if (!paste(m@name,'.fit',sep='') %in% ls(envir=m@.temp.env)) fn <- paste(m@name,'.fit',sep='') ####
#         else stop('the user defined function in fitFunction cannot be registered because an object with a similar name exists in the container!')
#         
#         assign(fn,fitFunction,envir = m@.temp.env) ####
#         m@packages <- unique(c(m@packages,'.temp'))
#       } else {
#         if (as.character(substitute(fitFunction))[1] == "::") fn <- as.character(substitute(fitFunction))[3]
#         else stop('fitFunction cannot be identified!')
#       }
#       
#       m@fitFunction <- eval(parse(text=paste("function(params) {do.call(",fn,",params)}")))
#     } else stop('fitFunction cannot be identified!')
#     
#   } 
#   #---------
#   
#   if (!is.null(predictFunction)) {
#     if (exists(as.character(substitute(predictFunction)),mode='function')) {
#       
#       if (class(substitute(predictFunction)) == 'name') {
#         fn <- as.character(substitute(predictFunction))
#         if (environmentName(environment(predictFunction)) == "R_GlobalEnv" || environmentName(environment(predictFunction)) == "sdm") {
#           # assign to the environment in the container of methods!
#           if (is.null(m@.temp.env)) m@.temp.env <- new.env()
#           assign(fn,predictFunction,envir = m@.temp.env) ####
#           m@packages <- unique(c(m@packages,'.temp'))
#           # when environment if attached, the conflict with the existing object in .GlobalEnv should be resolved!
#         } 
#         m@predictFunction <- eval(parse(text=paste("function(params) {do.call(",fn,",params)}")))
#         
#       } else if (class(substitute(predictFunction)) == 'call') {
#         
#         if (environmentName(environment(predictFunction)) == "R_GlobalEnv") {
#           if (is.null(m@.temp.env)) m@.temp.env <- new.env()
#           if (!paste(m@name,'.predict',sep='') %in% ls(envir=m@.temp.env)) fn <-paste(m@name,'.predict',sep='') ####
#           else stop('the user defined function in predictFunction cannot be registered because an object with a similar name exists in the container!')
#           assign(fn,predictFunction,envir = m@.temp.env) ####
#           m@packages <- unique(c(m@packages,'.temp'))
#         } else {
#           if (as.character(substitute(predictFunction))[1] == "::") fn <- as.character(substitute(predictFunction))[3]
#           else stop('predictFunction cannot be identified!')
#         } 
#         
#         m@predictFunction <- eval(parse(text=paste("function(params) {do.call(",fn,",params)}")))
#       } else stop('predictFunction cannot be identified!')
#       
#     } else stop(paste("Function",as.character(substitute(predictFunction)),"is not found!"))
#     
#     if (!is.null(predictParams)) {
#       if (is.list(predictParams)) {
#         n <- names(predictParams)
#         if (is.null(n)) stop('predictParams is not appropriately defined; example: list(newdata="sdmDataFrame")')
#         m@predictParams <- predictParams
#       } else stop('predictParams should be a list')
#     }
#     
#     if (!is.null(predictSettings)) {
#       if (is.list(predictSettings)) m@predictSettings <- predictSettings
#       else stop('predictSettings should be a list!')
#     }
#   }
#   #-------- 
#   if (!is.null(tuneParams)) {
#     if (!is.list(tuneParams)) stop('tuneParams should be a list; example: list(ntrees=seq(500,3000,by=200))')
#     n <- names(tuneParams)
#     if (is.null(n)) stop('tuneParams is not appropriately defined; example: list(ntrees=seq(500,3000,by=200))')
#     m@tuneParams <- tuneParams
#   } 
#   #------------  
#   if (inherits(metadata,'.Metadata')) m@metadata <- metadata
#   else  me <- .newMetadata(...)
#   if(!is.null(me)) m@metadata <- me
#   m
# }
#-----------
.movEnv <- function(e1,e2) {
  n1 <- ls(envir = e1)
  for (n in n1) assign(n,e1[[n]],envir = e2)
  rm(list=n1,envir = e1)
  e2
}
#----------
.movEnv2sdm <- function(e) {
  n1 <- ls(envir = e)
  for (n in n1) {
    if (!exists(n,envir = as.environment("package:sdm"))) assign(n,e[[n]],envir = as.environment("package:sdm"))
  } 
  rm(list=n1,envir = e)
}
#---------

.newMetadata <- function(...) {
  dot <- list(...)
  if (length(dot) > 0) {
    w <- unlist(lapply(dot, function(x) inherits(x,'.Metadata')))
    if (any(w)) return(dot[[which(w)]])
    else {
      ndot <- unlist(lapply(tolower(names(dot)),function(x) paste(strsplit(x,'')[[1]][1:3],collapse='')))
      w <- ndot %in% c('ful','cre','aut','web','cit','hel','des','dat','lic','url')
      if (any(w)) {
        m <- new(".Metadata")
        
        w <- which(ndot == 'aut')
        if (length(w) > 0) {
          if (is.list(dot[[w]])) m@authors <- dot[[w]]
          else if (is.character(dot[[w]])) m@authors <- list(dot[[w]])
        }
        
        w <- which(ndot == 'cre')
        if (length(w) > 0) {
          if (is.list(dot[[w]])) m@creators <- dot[[w]]
          else if (is.character(dot[[w]])) m@creators <- list(dot[[w]])
        }
        
        w <- which(ndot == 'tit')
        if (length(w) > 0) m@title <- dot[[w]]
        
        w <- which(ndot == 'url' | ndot == 'web')
        if (length(w) > 0) m@url <- dot[[w]]
        
        w <- which(ndot == 'cit')
        if (length(w) > 0) m@citations <- dot[[w]]
        
        w <- which(ndot == 'hel')
        if (length(w) > 0) m@Help <- dot[[w]]
        
        w <- which(ndot == 'des')
        if (length(w) > 0) m@description <- dot[[w]]
        
        w <- which(ndot == 'dat')
        if (length(w) > 0) m@date <- dot[[w]]
        
        w <- which(ndot == 'lic')
        if (length(w) > 0) m@license <- dot[[w]]
        return(m)
      }
    }
  }
}


#################################################

###############
.sdmMethods <- new('.sdmMethodsContainer')


# 
# #################
# .sdmMethods$getFitFunctions(c('glm','brt'))
# 
# 

# 
# 
# 
# 
# 
# 
# 
# m@fitFunction
# m@predictFunction
# m@.temp.env
# m@metadata
# m@settingRules
# m <- eval(m)
# saveRDS(m,'glm_model.rds')
# saveRDS(.create.sdmCorrelativeMethod(name='glm',modelTypes = 'pa',fitParams = list(formula='sdmFormula',data='sdmDataFrame'),
#                                      fitSettings = list(family=binomial(link='logit'),weights=NULL,model=FALSE),fitFunction = glm,
#                                      predictParams=list(object='sdmModel',newdata='sdmDataFrame'),
#                                      predictSettings=list(type='response'),predictFunction=predict.glm),file='glm_model.rds')
# 
# readRDS('glm_model.rds')
# save(m,file='glm_model.RData')
# eval(get(load('glm_model.RData')))
# rm(m)
# sa
# m
# 
# #############
# m <- .create.sdmCorrelativeMethod(name='glm',modelTypes = 'pa',fitParams = list(formula='sdmFormula',data='sdmDataFrame'),
#                                   fitSettings = list(family=binomial(link='logit'),weights=NULL,model=FALSE),fitFunction = glm,
#                                   predictParams=list(object='sdmModel',newdata='sdmDataFrame'),
#                                   predictSettings=list(type='response'),predictFunction=predict.glm)
# 
# 
# a <- .sdmMethods$new()
# 
# a$addMethod(m)
# a$getPredictFunctions()


