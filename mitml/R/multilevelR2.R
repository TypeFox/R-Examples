multilevelR2 <- function(model, print=c("RB1","RB2","SB","MVP")){

  # print argument case insensitive
  print <- toupper(print)
  print <- match.arg(print, several.ok=TRUE)
  method <- NULL

  # select method
  cls <- ifelse(is.list(model), class(model[[1]]), class(model))
  if(any(grepl("^l?merMod$",cls))) method <- "lmer"
  if(any(grepl("^lme$",cls))) method <- "nlme"
  if(is.null(method)) stop("Calculation of R-squared statistics not supported for models of class")

  # calculate R-squared
  if(is.list(model)){

    out <- sapply(model, .getRsquared, print=print, method=method)
    if(is.null(dim(out))) out <- matrix(out,nrow=1)
    out <- rowMeans(out)
  
  }else{

    out <- .getRsquared(model, print, method)

  }
  out

}

.getRsquared <- function(model, print, method){
# R squared for single model fit (lme4)

  # check if refit is necessary
  refit <- any(c("RB1","RB2","SB")%in%print)
  
  if(method=="lmer"){

    # model terms
    trm <- terms(model)
    if(!as.logical(attr(trm,"intercept"))) stop("Model must contain intercept.")
    yvr <- as.character(attr(trm,"variables")[-1])[attr(trm,"response")]
    cvr <- names(lme4::getME(model,"flist"))
    if(length(cvr)>1) stop("Calculation of R-squared only support for models with a single cluster variable.")
    cvr <- cvr[1]
  
    if(refit){
      
      # fit null model
      fml0 <- formula(paste0(yvr,"~1+(1|",cvr,")"))
      model0 <- update(model, fml0)
      
      # variance components under null
      vc0 <- lme4::VarCorr(model0)
      s0 <- attr(vc0,"sc")^2 
      t0.0 <- vc0[[cvr]][1,1]
      
    }
  
    # alternative model components
    beta <- lme4::fixef(model)[-1]
    X <- lme4::getME(model,"X")[,-1,drop=F]
    Z <- lme4::getME(model,"mmList")[[1]][,-1,drop=F]
    muX <- colMeans(X)
    muZ <- colMeans(Z)
    vZ <- cov(Z)
      
    # predicted and total variance
    vc1 <- lme4::VarCorr(model)
    t0.1 <- vc1[[cvr]][1,1]
    t10.1 <- vc1[[cvr]][1,-1]
    t11.1 <- vc1[[cvr]][-1,-1,drop=F]
    s1 <- attr(vc1,"sc")^2

  }

  if(method=="nlme"){
  
    # model terms
    trm <- terms(model)
    if(!as.logical(attr(trm,"intercept"))) stop("Model must contain intercept.")
    yvr <- as.character(attr(trm,"variables")[-1])[attr(trm,"response")]
    cvr <- attr(nlme::getGroups(model),"label")
    if(length(nlme::getGroupsFormula(model,asList=T))>1) stop("Calculation of R-squared only support for models with a single cluster variable.")
  
    if(refit){
      
      # fit null model
      ffml0 <- formula(paste0(yvr,"~1"))
      rfml0 <- formula(paste0("~1|",cvr,""))
      if(is.null(nlme::getData(model))) stop("No data sets found in 'lme' fit. See '?testModels' for an example.")
      model0 <- update(model, fixed=ffml0, random=rfml0, data=model$data)
      
      # variance components under null
      vc0 <- nlme::getVarCov(model0)
      s0 <- model0$sigma^2
      t0.0 <- vc0[1,1]
      
    }
   
    # alternative model components
    beta <- nlme::fixef(model)[-1]
    fe <- model$terms
    X <- model.matrix(fe,nlme::getData(model))[,-1,drop=F]
    re <- attr(model$modelStruct$reStruct[[1]],"formula")
    Z <- model.matrix(re,nlme::getData(model))[,-1,drop=F]
    muX <- colMeans(X)
    muZ <- colMeans(Z)
    vZ <- cov(Z)
      
    # predicted and total variance
    vc1 <- nlme::getVarCov(model)
    t0.1 <- vc1[1,1]
    t10.1 <- vc1[1,-1]
    t11.1 <- vc1[-1,-1,drop=F]
    s1 <- model$sigma^2
  
  }

  # calculate R2
  vyhat <- var( X %*% beta ) 
  vy <- vyhat + t0.1 + 2*(muZ %*% t10.1) + muZ%*%t11.1%*%muZ + sum(diag(t11.1%*%vZ)) + s1

  if(refit){
    rb1 <- 1 - s1/s0
    rb2 <- 1 - t0.1/t0.0
    sb <-  1 - (s1+t0.1)/(s0+t0.0)
  }else{
    rb1 <- rb2 <- sb <- NA
  }
  mvp <- as.vector(vyhat/vy)

  c(RB1=rb1, RB2=rb2, SB=sb, MVP=mvp)[print]

}
