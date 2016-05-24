anovageePrim2 <- function(m1, m2,...){

  mm1 <- model.matrix(m1)
  mm2 <- model.matrix(m2)

  P1 <- mm1 %*% solve(t(mm1)%*%mm1) %*% t(mm1) 
  P2 <- mm2 %*% solve(t(mm2)%*%mm2) %*% t(mm2)
  e2 <- mm2 - P1 %*% mm2
  e1 <- mm1 - P2 %*% mm1

  m2inm1 <- all(apply(e2,2,var) < 1e-10)
  m1inm2 <- all(apply(e1,2,var) < 1e-10)
  
  if (!any(c(m2inm1,m1inm2)))
    cat("Models not nested\n")
  else 
    if (all(c(m2inm1,m1inm2)))
      cat("Models are identical\n")
    else {
      if (m1inm2){
        tmp <- m1
        m1 <- m2
        m2 <- tmp
      }
      ## Now mm2 < mm1
      mm1 <- model.matrix(m1)
      mm2 <- model.matrix(m2)

      ## What is this? I wonder
      mf1 <- paste(paste(formula(m1))[c(2,1,3)],collapse=" ")
      mf2 <- paste(paste(formula(m2))[c(2,1,3)],collapse=" ")

      ## Reparametrize the model
      mm <- cbind(mm2,mm1)
      qmm <- qr(mm)
      qmmq <- qr.Q(qmm)
      nymm1 <- as.data.frame(qmmq[,1:qmm$rank])
      colnames(nymm1) <- paste("parm",1:ncol(nymm1),sep=".")
      nymm2 <- nymm1[,1:ncol(mm2),drop=FALSE]

      formula1 <- formula(paste(formula(m1)[[2]],formula(m1)[[1]],
                                paste(c("-1",colnames(nymm1)),collapse="+"),collapse=""))
            
      m1call <- m1$call

      ## BUGFIX provided by Stefan Boehringer
      ##nymm1[,paste(formula(m1call)[[2]])] <- m1$y
      nymm1[, paste(formula(m1)[[2]])] <- m1$y

      nymm1[,paste(m1call$id)] <- m1$id
      m1call$offset <- m1$offset
      m1call$weights <- m1$weights
      m1call$formula <- formula1
      m1call$data <- nymm1
      m1ny <- eval(m1call)

      ## Calculate wald statistic
      beta <- coef(m1ny)
      vbeta <- summary(m1ny)$cov.unscaled
      df<- dim(mm1)[2]-dim(mm2)[2]
      rbeta<-rep(1,length(beta))
      rbeta[1:df]<-0
      beta0<-rev(rbeta)
      zeroidx <- beta0==0
      X2<-t(beta[zeroidx])%*% solve(vbeta[zeroidx,zeroidx,drop=FALSE])%*%beta[zeroidx]
      

      ## Make table with results
      topnote <- paste("Model 1", mf1,"\nModel 2", mf2)
      title <- "Analysis of 'Wald statistic' Table\n"      
      table <- data.frame(Df=df, X2=X2,p=1-pchisq(X2,df))
      dimnames(table) <- list("1", c("Df", "X2", "P(>|Chi|)"))      
      val <- structure(table, heading = c(title,topnote), class = c("anova", 
                                                      "data.frame"))
      return(val)
    }
}

anova.geeglmlist <- 
  function (object, ..., dispersion = NULL, test = NULL) 
{
  responses <- as.character(lapply(object, function(x) {
    deparse(formula(x)[[2]])
  }))
  sameresp <- responses == responses[1]
  if (!all(sameresp)) {
    object <- object[sameresp]
    warning("Models with response ", deparse(responses[!sameresp]), 
            " removed because response differs from ", "model 1")
  }
  
  ns <- sapply(object, function(x) length(x$residuals))
  if (any(ns != ns[1])) 
    stop("models were not all fitted to the same size of dataset")
  
  objects <- list(object,...)    
  m1 <- objects[[1]][[1]]
  if (length(objects[[1]])>1)
    m2 <- objects[[1]][[2]]
  else 
    m2 <- NULL
  
  value <- anovageePrim2(m1,m2)
  return(value)
}


anova.geeglm<-function (object, ..., dispersion = NULL, test = NULL) 
{
  dotargs <- list(...)
  named <- if (is.null(names(dotargs))) 
    rep(FALSE, length(dotargs))
  else (names(dotargs) != "")
  if (any(named)) 
    warning("The following arguments to anova.glm(..) are invalid and dropped: ", 
            paste(deparse(dotargs[named]), collapse = ", "))
  dotargs <- dotargs[!named]
  is.glm <- unlist(lapply(dotargs, function(x) inherits(x, "glm")))
  dotargs <- dotargs[is.glm]
  if (length(dotargs) > 0) 
    return(anova.geeglmlist(c(list(object), dotargs), dispersion = dispersion, 
                            test = test))
  
  varlist <- attr(object$terms, "variables")
  ##print(varlist)
  x <- if (n <- match("x", names(object), 0)) 
    object[[n]]
  else model.matrix(object)
  
  varseq <- attr(x, "assign")
  
  nvars <- max(0, varseq)
  betaList <- vbetaList <- NULL
  
  if (nvars > 1) {
    method <- object$method
    if (!is.function(method)) 
      method <- get(method, mode = "function", envir = parent.frame())
    for (i in 1:(nvars - 1)) {
      eprint("calling fit....")
      ##print(length(object$y))
      fit <- method(x = x[, varseq <= i, drop = FALSE], 
                    y = object$y, weights = object$prior.weights, 
                    corstr = object$corstr,
                    start = object$start, offset = object$offset,   id=object$id,
                    family = object$family, control = object$control)
      
      betaList <- c(betaList,list(fit$beta))
      vbetaList <- c(vbetaList,list(fit$vbeta))
    }
  }
  
  betaList <- c(betaList, list( object$geese$beta ))
  vbetaList <- c(vbetaList, list( object$geese$vbeta ))
  
  hasIntercept <- (length(grep("(Intercept)",names(betaList[[1]])))!=0)
  dimVec <- unlist(lapply(betaList,length))

  if (hasIntercept){
    dfVec <- dimVec[1]-1
  } else {
    dfVec <- dimVec[1]
  }

  if (length(dimVec)>1){
    for (i in 2:length(dimVec))
      dfVec <- c(dfVec,dimVec[i]-dimVec[i-1])
  }
  
  ##print(dfVec)
  X2Vec <- NULL
  ## Calculate Wald statistics
  for (i in 1:length(dfVec)){
    beta <- betaList[[i]]
    vbeta <- vbetaList[[i]]
    beta0 <- rep(1,length(beta))
    beta0[1:dfVec[i]] <- 0
    beta0 <- rev(beta0)
    zeroidx <- beta0==0
    X2 <- t(beta[zeroidx])%*%solve(vbeta[zeroidx,zeroidx,drop=FALSE])%*%beta[zeroidx]
    X2Vec <- c(X2Vec,X2)
  }
  
  resdf <- dfVec
  resdev <- X2Vec
  table <- data.frame(resdf, resdev, 1-pchisq(resdev,resdf))
  
  tl <- attr(object$terms, "term.labels")

  #print(table)
  if (length(tl) == 0) 
    table <- table[1, , drop = FALSE]
  
  dimnames(table) <- list(c(tl), c("Df", "X2", "P(>|Chi|)"))
  
  title <- paste("Analysis of 'Wald statistic' Table", "\nModel: ", 
                 object$family$family, ", link: ", object$family$link, 
                 "\nResponse: ", as.character(varlist[-1])[1],
                 "\nTerms added sequentially (first to last)\n", 
                 sep = "")
  
  structure(table, heading = title, class = c("anova", "data.frame"))
}





# anova.geeglm <- function(object, ...){
#   anovaPgee (object, ...)
# }


# anovaPgee <- function(object, ...){
#   #cat("anova.gee\n")
#   m1 <- object
#   objects <- list(object,...)
#   if (length(objects)>1)
#     m2 <- objects[[2]]
#   else 
#     m2 <- NULL
  
#   if (is.null(m2)){
#     term <- attr(object$terms,"term.labels")
#     resp <- paste(formula(object))[2]
#     rhs  <- lapply(1:length(term), function(i) paste(term[1:i],collapse=" + "))
#     print(rhs)
    
#     model.list <- c(paste(resp,"~ 1"), paste(resp,"~", rhs))

    
#     value <- NULL
#     for (i in 2:length(model.list)){
#       if (i==2){
#         mf1 <- model.list[i-1]
#         mf2 <- model.list[i]
#         ##print(mf1); print(mf2)

#         #print(mf1)
#         #print(object)
        
#         m1 <- update(object,formula=as.formula(mf1))
#         m2 <- update(object,formula=as.formula(mf2))
#       } else {
#         m1 <- m2
        
#         m2 <- update(object,formula=as.formula(model.list[i]))
#         ##print(formula(m1)[1:3]); print(formula(m2)[1:3])
#       }
#       value <- rbind(value,anovageePrim(m1,m2))
#     }
#     rownames(value) <- term
#     attr(value,"model1") <- NULL
#     attr(value,"model2") <- NULL
#   } else {
#     value <- anovageePrim(object,m2)
#   }
#   value[,3] <- round(value[,3],5)
#   return(value)
# }


# anovageePrim <- function(m1, m2,...){
#   mm1 <- model.matrix(m1)
#   mm2 <- model.matrix(m2)
#   P1 <- mm1 %*% solve(t(mm1)%*%mm1) %*% t(mm1) 
#   P2 <- mm2 %*% solve(t(mm2)%*%mm2) %*% t(mm2)
#   e2 <- mm2 - P1 %*% mm2
#   e1 <- mm1 - P2 %*% mm1

#   #print(mm1[c(1:5,100:105),]); print(mm2[c(1:5,100:105),])
#   m2inm1 <- all(apply(e2,2,var) < 1e-10)
#   m1inm2 <- all(apply(e1,2,var) < 1e-10)

#   #print(apply(e2,2,var))
#   #print(apply(e1,2,var))
#   #print(m2inm1)
#   #print(m1inm2)
#   if (!any(c(m2inm1,m1inm2)))
#     cat("Models not nested\n")
#   else 
#     if (all(c(m2inm1,m1inm2)))
#       cat("Models are identical\n")
#     else {
#       if (m1inm2){
#         tmp <- m1
#         m1 <- m2
#         m2 <- tmp
#       }
      
#       mm1 <- model.matrix(m1)
#       mm2 <- model.matrix(m2)
      
#       mf1 <- paste(paste(formula(m1))[c(2,1,3)],collapse=" ")
#       mf2 <- paste(paste(formula(m2))[c(2,1,3)],collapse=" ")

#       mm <- cbind(mm2,mm1)
#       qmm <- qr(mm)
#       qmmq <- qr.Q(qmm)
#       nymm1 <- as.data.frame(qmmq[,1:qmm$rank])
#       colnames(nymm1) <- paste("parm",1:ncol(nymm1),sep=".")
#       nymm2 <- nymm1[,1:ncol(mm2),drop=FALSE]

#       dimDiff <- ncol(nymm1)-ncol(nymm2)

#       D <- diag(dimDiff)
#       L <- cbind(matrix(0,ncol=ncol(nymm2),nrow=nrow(D)),D)

#       formula1 <- formula(paste(formula(m1)[[2]],formula(m1)[[1]],
#                                 paste(c("-1",colnames(nymm1)),collapse="+"),collapse=""))
      

      
#       m1call <- m1$call
#       #print(formula(m1call)[[2]])
#       #print(nymm1[1:10,])
#       #print(paste(m1call$formula[[2]]))
#       #nymm1[,paste(m1call$formula[[2]])] <- m1$y
#       nymm1[,paste(formula(m1call)[[2]])] <- m1$y
#       nymm1[,paste(m1call$id)] <- m1$id
#       m1call$offset <- m1$offset
#       m1call$weights <- m1$weights
#       m1call$formula <- formula1

#       m1call$data <- nymm1
#       m1ny <- eval(m1call)
#       #print(class(m1ny))

#       val <- esticon(m1ny,L,joint.test=TRUE)

#       rownames(val)<-""
#       class(val) <- c("anova.gee","data.frame")
#       attr(val,"model1") <- mf1
#       attr(val,"model2") <- mf2
#       return(val)
#     }
# }

# print.anova.geeglm <- function(x,...){
#   cat("Analysis table for GEE models\n\n")
#   if (!is.null(attr(x,"model1"))){
#     cat("Model 1: "); cat(attr(x,"model1"), "\n")
#     cat("Model 2: "); cat(attr(x,"model2"), "\n\n")
#   }
#   print.data.frame(x)
  
# }


