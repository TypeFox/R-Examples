# Author: Babak Naimi, naimi.b@gmail.com
# Date :  March. 2016
# Version 1.0
# Licence GPL v3

#-------------
.detachPackage <- function(n,unload=TRUE,force=TRUE) {
  s <- search()
  ss <- unlist(lapply(strsplit(s,":"),function(x) x[2]))
  for (nn in n) {
    if (nn %in% ss) {
      detach(s[which(ss == nn)],force=force,character.only=TRUE,unload=unload)
    }
  }
}

.maxent <- function(formula,data,beta,prevalence,feat,args,...) {
  if (!is.null(.sdmOptions$getOption('maxJar')) && !.sdmOptions$getOption('maxJar')) {
    jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
    if (file.exists('maxent.jar')) {
      if (!file.copy('maxent.jar',jar)) stop('maxent.jar file should be copied in the right place in the package dismo; check the help page of maxent in dismo package!')
      else {
        .sdmOptions$addOption('maxJar',file.exists(jar))
        .detachPackage(c('dismo','rJava'))
        if (!.require('dismo') & !.require('rJava')) stop('package dismo OR rJava does not exist or could not be loaded!')
      }
    } else stop(paste('maxent.jar does not exist in the working directory or in ', paste(system.file(package="dismo"), "/java/", sep=''),'...!',sep=''))
  }
  options(java.parameters = "-Xmx1g" )
  #--------
  if (!missing(feat) && is.vector(feat)) {
    w <- which(substr(feat,1,1) == '-')
    if (length(w) > 0) feat[w] <- substring(feat[w],2)
    
    feat <- unique(.pmatch(feat,c('linear','quadratic','product','hinge','autofeature','threshold')))
    feat <- feat[!is.na(feat)]
    if (length(feat) == 0) feat <- NULL
    else if ('autofeature' %in% feat) feat <- '-A'
  } else feat <- NULL
  
  if (!missing(beta) && is.numeric(beta)) {
    beta <- paste0('betamultiplier=',beta)
  } else beta <- NULL
  if (!missing(prevalence) && is.numeric(prevalence)) {
    prevalence <- paste0('defaultprevalence=',prevalence)
  } else prevalence <- NULL
  if(missing(args) || !is.character(args)) args <- NULL
  
  args <- c(args,feat,beta,prevalence)
  #--------
  sp <- deparse(formula[[2]])
  w <- colnames(data) != sp
  dismo::maxent(x=data[,w,drop=FALSE],p=data[,!w,drop=TRUE],args=args,...)
}