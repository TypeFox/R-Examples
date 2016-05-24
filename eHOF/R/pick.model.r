pick.model <-  function(...) UseMethod("pick.model")

pick.model.HOF.list <- function(object, ...) {
   out <- sapply(object, function(x) pick.model.HOF(object=x,  ...))
  return(out)
}

pick.model.HOF <- function (object, level = 0.95, test = c('AICc', 'BIC', 'AIC','Dev'), modeltypes, penal = 'df', gam = FALSE, selectMethod = c('bootselect.lower', 'bootselect.always', 'IC.weight', 'pick.model'), silent = FALSE, ...) {
  selectMethod <- match.arg(selectMethod)
  if(is.null(object$bootstrapmodels) & selectMethod != 'pick.model') {
    if(getOption("eHOF.bootselectmessage")) {
       message('Bootselect or IC.weight method only possible after bootstrapping.')
       options(eHOF.bootselectmessage = FALSE)}
    selectMethod <- 'pick.model'
  }
  test <- match.arg(test)
  if(missing(modeltypes)) modeltypes <- eHOF.modelnames
  if(length(penal)==1) {
    penal <- switch(penal,
    df= sapply(object$models, function(x) length(x$par)),
    h = c(I=1, II=2, III=3, IV=4, V=5, VI=6, VII=7),
    m = c(I=1, II=2, III=2, IV=3, V=4, VI=4, VII=5)
    )
  }
  rejectedmodels <- sapply(object$models, function(x) is.na(x$deviance))
  if(gam) {
    gamfun <- function(x, bs = 'cr', k = -1, ...) gam(x$y ~ s(x$x, bs=bs, k=k),family = get(x$family), scale = 0, ...)
    pg <- gamfun(object, ...)
    object$models$GAM <- pg
    modeltypes <- c(modeltypes, 'GAM')
    if('GAM' %in% names(penal)) penal['GAM'] <- sum(pg$edf) else 
      penal <- c(penal, GAM=sum(pg$edf))
  }
  
#   if(gam) {
#     x <- object$x; y <- object$y; fam <- get(object$family)
#   	gamfun <- function(..., y ~ s(x, bs=bs, k=k), family = fam, scale = 0, bs = 'cr', k = -1, ...) gam(...)
#     pg <- gamfun(object = object, bs = bs, k = k, ...)
#   	object$models$GAM <- pg
#   	modeltypes <- c(modeltypes, 'GAM')
#   	if('GAM' %in% names(penal)) penal['GAM'] <- sum(pg$edf) else 
#   	penal <- c(penal, GAM=sum(pg$edf))
#   }
#    
  if (test == "AIC") {
      k <-  2
      p <- penal[match(names(penal), names(object$models))]
      ic <- -2 * logLik(object) + k * p
      ic <- ic[names(ic) %in% modeltypes]
      model <- names(ic[!rejectedmodels])[which.min(ic[!rejectedmodels])]
  }
  if (test == "AICc") {
      k <-  2
      p <- penal[match(names(penal), names(object$models))]
      ic <- -2 * logLik(object) + k * p + (2*k*(k + 1))/(object$nobs - k - 1)
      ic <- ic[names(ic) %in% modeltypes]
      model <- names(ic[!rejectedmodels])[which.min(ic[!rejectedmodels])]
  }
  if( test == 'BIC') {
      k <- log(object$nobs)
      p <- penal[match(names(penal), names(object$models))]
      ic <- -2 * logLik(object) + k * p
      ic <- ic[names(ic) %in% modeltypes]
      model <- names(ic[!rejectedmodels])[which.min(ic[!rejectedmodels])]
  }
  if (test == "Dev") {
      ic <- deviance(object)
      ic <- ic[names(ic) %in% modeltypes]
      model <- names(ic[!rejectedmodels])[which.min(ic[!rejectedmodels])]
  }

  if(selectMethod == 'bootselect.lower') {
    # modboot <- names(ic[!rejectedmodels])[which.min(ic[!rejectedmodels])]    
    modboot <- names(which.max(table(object$bootstrapmodels)))
    if(eHOF.modelnames[match(modboot, eHOF.modelnames)] < eHOF.modelnames[match(model, eHOF.modelnames)] )
      if(!silent) message(object$y.name, ': Most frequent bootstrap model (',modboot, ') lower than original choice (', model, ',',  object$bootstraptest,') by ', test, ' test.', sep='')
      if(!modboot %in% names(rejectedmodels)[rejectedmodels]) model <- modboot
  }
  if(selectMethod == 'bootselect.always') {
    modboot <- names(which.max(table(object$bootstrapmodels)))
    if(modboot != model )
      if(!silent) message(object$y.name, ': Most frequent bootstrap model (',modboot, ') not equal to original choice (', model, ',',  object$bootstraptest,') by ', test, ' test.', sep='')
    if(!modboot %in% names(rejectedmodels)[rejectedmodels]) model <- modboot
  }
if(selectMethod == 'IC.weights') {
#     cat('Stabilisation test for models V and VII: ')
#     print(table(object$bootstrapmodels)[modboot])
#     if(table(object$bootstrapmodels)[modboot] <= length(object$bootstrapmodels)/2) {
#     cat(': Chosen model (', modboot, ') does not exceed 50%, ', sep='')
      dev <- deviance(object); ll <- logLik(object)
      AICc <- -2 * ll + 2 * penal + 2 * penal *(penal + 1)/(object$nobs - penal - 1) 
      d.AICc <- AICc - min(AICc, na.rm=TRUE)
      AICc.W <- round(exp(-0.5*AICc)/ sum(exp(-0.5*AICc), na.rm=TRUE),4)
      model.weight <- names(which.max(AICc.W[!rejectedmodels]))
      if(eHOF.modelnames[match(model.weight, eHOF.modelnames)] < eHOF.modelnames[match(model, eHOF.modelnames)]) {
        if(!silent) message('Selection method: Highest ', test, ' weights.\n', object$y.name, ': Original model choice (', model, ') lower than the model with highest IC weight (', model.weight, ') is chosen instead.', sep='')
        if(!model.weight %in% names(rejectedmodels)[rejectedmodels]) model <- model.weight
	  }
  }
  return(model)
}
