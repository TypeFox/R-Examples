
#' Predict Conditional Effects
#' 
#' Predicts conditional treatment effects based on a fitted EffectLiteR model.
#' 
#' @param obj Object of class \code{effectlite}.
#' @param newdata An optional data.frame, containing the same continuous and 
#' categorical covariates as used when fitting the EffectLiteR model in
#' obj. Only covariates (and neither the dependent variable nor indicators for 
#' latent variables) should be included.
#' @return Object of class \code{"data.frame"}.
#' @examples
#' m1 <- effectLite(y="dv", z=c("z1"), k=c("k1","kateg2"), x="x", 
#' control="control", data=example01)
#' newdata <- data.frame(k1="male", kateg2="1", z1=2)
#' elrPredict(m1, newdata)
#' @export
elrPredict <- function(obj, newdata=NULL){
  
  stopifnot(inherits(obj, "effectlite"))
  
  ##TODO merge function with computeConditionalEffects()
  if(is.null(newdata)){return(obj@results@condeffects)}
  
  z <- obj@input@vnames$z
  k <- obj@input@vnames$k
  stopifnot(all(names(newdata) %in% c(z,k)))
  
  ##TODO error check very important!
  ##TODO add documentation
  ##TODO do the above things before making it public!
  
  ## required things
  lavresults <- obj@results@lavresults
  nk <- obj@input@nk
  nz <- obj@input@nz
  ng <- obj@input@ng
  
  
  #compute Kstar values first
  if(nk > 1){
    tmp <- obj@input@vlevels$levels.k.original
    tmp <- tmp[length(tmp):1]
    tmp <- expand.grid(tmp)
    tmp$kstar <- factor(obj@input@vlevels$kstar)
    
    ## add Kstar values to newdata
    newdata <- merge(newdata, tmp)
  }
  
  ## estimates and vcov
  est <- parameterEstimates(lavresults)$est ## parameter estimates
  names(est) <- parameterEstimates(lavresults)$label 
  vcov <- lavInspect(lavresults, "vcov.def", add.class = FALSE)
  
  ## compute formula and model.matrix  
  if(nz==0 & nk==1){
    formula <- as.formula(" ~ 1")
    modmat <- model.matrix(formula, data=newdata)
    kz <- "00"
    
  }else if(nz>0 & nk==1){
    formula <- as.formula(paste0(" ~ ", paste(z, collapse=" + ")))
    modmat <- model.matrix(formula, data=newdata)
    kz <- paste0("0",0:nz)
    
  }else if(nz==0 & nk>1){      
    formula <- as.formula(" ~ kstar")
    modmat <- model.matrix(formula, data=newdata)      
    kz <- paste0(1:nk-1,"0")
    
  }else if(nz>0 & nk>1){ 
    formula <- as.formula(paste0(" ~ ", 
                                 paste("kstar", z, sep="*", collapse=" + ")))
    modmat <- model.matrix(formula, data=newdata)            
    kz <- c(paste0(1:nk-1,"0"), paste0("0",1:nz))
    kz <- c(kz, paste0(rep(1:(nk-1),nz), rep(1:nz, each=nk-1)))
    
  }
  
  estimates <- est[paste0("g1",kz)]
  vcov_est <- vcov[paste0("g1",kz),paste0("g1",kz)]
  individualeffects <- cbind(modmat %*% estimates)
  individualeffects <- cbind(individualeffects,
                             apply(modmat,1,function(x){sqrt(t(x) %*% vcov_est %*% x)}))
  
  
  if(ng > 2){
    for(i in 3:ng){
      estimates <- est[paste0("g",i-1,kz)]
      vcov_est <- vcov[paste0("g",i-1,kz),paste0("g",i-1,kz)]
      individualeffects <- cbind(individualeffects, modmat %*% estimates)
      individualeffects <- cbind(individualeffects,
                                 apply(modmat,1,function(x){sqrt(t(x) %*% vcov_est %*% x)}))
    }      
  }  
  
  individualeffects <- as.data.frame(individualeffects)
  names(individualeffects) <- paste0(rep(c("","se_"), times=ng-1),
                                     "g",
                                     rep(2:ng-1, each=2))
  
  ## add true-outcomes
  estimates <- est[paste0("b0",kz)]
  trueoutcomes <- cbind(modmat %*% estimates)
  for(i in 1:(ng-1)){
    estimates <- est[paste0("b",i,kz)]
    trueoutcomes <- cbind(trueoutcomes, modmat %*% estimates)
  }
  trueoutcomes <- as.data.frame(trueoutcomes)
  names(trueoutcomes) <- paste0("Eygx", 0:(ng-1), "kz")
  individualeffects <- cbind(individualeffects,trueoutcomes)
  
  return(individualeffects)
}

## maybe add an option that true outcomes and orginal data is included...
computeConditionalEffects <- function(obj, est, vcov, m1){
  
  current.na.action <- options('na.action')
  on.exit(options(current.na.action))
  
  options(na.action='na.pass')
  
  ## required things
  z <- obj@input@vnames$z
  k <- obj@input@vnames$k
  x <- obj@input@vnames$x
  data <- obj@input@data  
  mm <- obj@input@measurement 
  nz <- obj@input@nz
  nk <- obj@input@nk
  ng <- obj@input@ng
  
  data$id <- 1:nrow(data)
  latentz <- z[which(!z %in% names(data))]
  
  ## add factor scores
  if(length(latentz) > 0){
    fscores <- data.frame(do.call("rbind", lavPredict(m1)))
    fscores <- subset(fscores, select=latentz)
    fscores$id <- unlist(lavInspect(m1, "case.idx"))
    data <- merge(data,fscores)
  }
  
  ## compute formula and model.matrix  
  if(nz==0 & nk==1){
    formula <- as.formula(" ~ 1")
    modmat <- model.matrix(formula, data=data)
    kz <- "00"
    dsub <- data.frame(matrix(vector(),nrow=nrow(data),ncol=0))
    
  }else if(nz>0 & nk==1){
    formula <- as.formula(paste0(" ~ ", paste(z, collapse=" + ")))
    modmat <- model.matrix(formula, data=data)
    kz <- paste0("0",0:nz)
    dsub <- data[,c(x,z)]
    
  }else if(nz==0 & nk>1){      
    formula <- as.formula(" ~ kstar")
    modmat <- model.matrix(formula, data=data)      
    kz <- paste0(1:nk-1,"0")
    dsub <- data[,c(x,"kstar",k)]
    names(dsub)[2] <- "K"
    
  }else if(nz>0 & nk>1){ 
    formula <- as.formula(paste0(" ~ ", 
                                 paste("kstar", z, sep="*", collapse=" + ")))
    modmat <- model.matrix(formula, data=data)            
    kz <- c(paste0(1:nk-1,"0"), paste0("0",1:nz))
    kz <- c(kz, paste0(rep(1:(nk-1),nz), rep(1:nz, each=nk-1)))
    dsub <- data[,c(x,"kstar",k,z)]
    names(dsub)[2] <- "K"
    
  }
  
  estimates <- est[paste0("g1",kz)]
  vcov_est <- vcov[paste0("g1",kz),paste0("g1",kz)]
  condeffects <- cbind(modmat %*% estimates)
  condeffects <- cbind(condeffects,
                       apply(modmat,1,function(x){sqrt(t(x) %*% vcov_est %*% x)}))
  
  
  if(ng > 2){
    for(i in 3:ng){
      estimates <- est[paste0("g",i-1,kz)]
      vcov_est <- vcov[paste0("g",i-1,kz),paste0("g",i-1,kz)]
      condeffects <- cbind(condeffects, modmat %*% estimates)
      condeffects <- cbind(condeffects,
                           apply(modmat,1,function(x){sqrt(t(x) %*% vcov_est %*% x)}))
    }      
  }  
  
  condeffects <- as.data.frame(condeffects)
  names(condeffects) <- paste0(rep(c("","se_"), times=ng-1),
                               "g",
                               rep(2:ng-1, each=2))
  condeffects <- cbind(dsub,condeffects)
  
  ## add true-outcomes
  estimates <- est[paste0("b0",kz)]
  trueoutcomes <- cbind(modmat %*% estimates)
  for(i in 1:(ng-1)){
    estimates <- est[paste0("b",i,kz)]
    trueoutcomes <- cbind(trueoutcomes, modmat %*% estimates)
  }
  trueoutcomes <- as.data.frame(trueoutcomes)
  names(trueoutcomes) <- paste0("Eygx", 0:(ng-1), "kz")
  condeffects <- cbind(condeffects,trueoutcomes)
  
  
  ## add variables used in the propscore model
  propscore <- obj@input@vnames$propscore
  if(!is.null(propscore)){
    
    d <- obj@input@data
    
    if(is(propscore, "formula")){      
      form <- propscore
    }else{
      form <- as.formula(paste0(x, " ~ ", paste0(propscore, collapse=" + ")))
    }
    
    dsub <- model.frame(form,data=d)
    condeffects <- condeffects[,-1]
    condeffects <- cbind(dsub, condeffects)
    
  }
  
  return(condeffects)
  
}

