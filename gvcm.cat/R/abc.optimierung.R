abcfit <-
function( 
x,
y,
weights,
family,
tuning = c("AIC", "BIC"), 
control,
indices,
start = rep(0, nvars),
offset = rep(0, n),
...
)

{

# definitions zu design + names
  x <- as.matrix(x)
  xnames <- dimnames(x)[[2L]]
  ynames <- if (is.matrix(y))
      rownames(y) else names(y)    
  conv <- FALSE
  n <- NROW(y)
  p <- nvars <- ncol(x)

# check gewichte/offset
    if (is.null(weights))
        weights <- rep.int(1, n)
    if (is.null(offset))
        offset <- rep.int(0, n)
    if (length(offset)==1)
        offset <- rep.int(offset, n)
    
# check start 
  initials <- if (!is.null(start)) { # start value fuer ceofs da
                  if (length(start) != nvars)   # falsche Laenge?
                      {stop(gettextf("length of 'start' should equal %d and correspond to initial coefs for %s",
                        nvars, paste(deparse(xnames), collapse = ", ")),
                        domain = NA)}
                  else { # keine falsche Laenge!
                      isnogood <- is.na(start)
                      if (any(isnogood)) { start[isnogood] <- 0 }
                      start }
            } else # kein start fuer ceofs da => default value 0!
            { start <- rep(0, nvars) 
              start
            } 

# definitions
  index1 <- indices["index1",]
  index2 <- colSums(as.matrix(indices[c("index2", "index3", "index4", 
                                            "index5", "index7", "index8"),]))
  # index2 - indicator for special v
  # index3 - indicator for special p
  # index4 - indicator for grouped 
  # index5 - indicator for grouped.fused
  # index6 - indicator for pspline
  # index7 - indicator for SCAD
  # index8 - indicator for elastic
  # index9 - indicator for mspline
                    

# definitions tuning
  if (tuning == "AIC") {criterion <- function(A){0} } else
                       {criterion <- function(A){(log(n)-2)*dim(A)[2]} }    
  assured.intercept <- control$assured.intercept

                     
# matrix for variavle selection 
   
    A.model <- matrix(nrow=ncol(x),ncol=0)
    rownames(A.model) <- xnames
    
    # building basic A.model: non-varying elements
      if (length(which(index2 == 0))!=0) {
      if ( which(index2 == 0)[1]==1) {
           A.model <- cbind(A.model, diag(p)[,1])
           b <- 2 } else { b <- 1 }
      if ( b <= length(which(index2 == 0)) ) {
           for (i in b:length(which(index2 == 0))) {
                ind <- (cumsum(index1)[which(index2 == 0)[i]-1]+1):(cumsum(index1)[which(index2 == 0)[i]])
                A.model <- cbind(A.model, diag(p)[,ind] ) }
          }
      }
    
    # building basic A.model: intercept varies, but has to be in the model
      lom <- rep(index2,times=index1) # level of measurement
      if (((index2[1]!=0) && (assured.intercept==TRUE))==TRUE) {
      A.model <- cbind(A.model, rep(c(1,0),c(index1[1],p-index1[1])) )  
      index2[1] <- 0
      }

# initial values    
  suppressWarnings(try(option <- glm.fit(x=x%*%A.model, y=y, weights = weights, 
          offset = offset, start = as.vector(t(A.model)%*%initials), 
          family = family, intercept = FALSE))) 
  if(exists("option")==TRUE) { abc.option <- option$aic + criterion(A.model) } else
                         { option <- NA
                           abc.option <- Inf 
                         } 
  abc.model <- abc.option
  A.option <- A.model
  A <- as.matrix(abc.a.coefs(list(index1, index2)))    
  A.models <- list()

# selecting and splitting
  w <- 0
  while(abc.option <= abc.model) {
      
      model <- option
      abc.model <- abc.option
      A.model <- A.option
      
      A.aspirants <- list()
      # 1. adding 1df via new coefficients
      b <- ncol(A)
      if (length(b)>0) { if (b!=0) {
      for (i in 1:b) {
           A.aspirant <- cbind(A.model, matrix(A[,i], ncol=1))
           A.aspirants[[i]] <- A.aspirant
          } } }
    
      # 2. adding 1 df via splitting blocks
      splitting.aspirants <- which(colSums(A.model)>1)
      if (length(splitting.aspirants)>0) {
      if(length(b)>0) { if (b!=0) {j <- b+1} else {j <- 1} } else {j <- 1} 
      for (i in 1:length(splitting.aspirants)) {
           dismissed <- as.matrix(A.model[,splitting.aspirants[i]])
           ind1 <- sum(dismissed)
           ind2 <- lom[which(dismissed==1)[1]]
           splittings.1 <- abc.a.coefs(list(ind1,ind2),TRUE) 
           splittings.2 <- 1-splittings.1
           for (k in 1:ncol(splittings.1)){
               split.1 <- dismissed
               split.1[which(dismissed==1)] <- splittings.1[,k]
               split.2 <- dismissed
               split.2[which(dismissed==1)] <- splittings.2[,k]
               A.aspirant <- cbind(A.model[,-splitting.aspirants[i]], split.1, split.2)
               A.aspirants[[j]] <- A.aspirant
               j <- j+1 }
          } }    
    
      # possible models
      abc <- c()
      if (length(A.aspirants)>0) {
          for (i in 1:length(A.aspirants)) {
               rm(option)             
               suppressWarnings(try(option <- glm.fit(x=x%*%A.aspirants[[i]], y=y, weights = weights, 
                         offset = offset, start = as.vector(t(A.aspirants[[i]])%*%initials), 
                         family = family, intercept = FALSE))) 
               if(exists("option")==FALSE) { stop ("There are problems computing the model with a subset of coefficients only. \n") }
               abc <- c(abc, option$aic + criterion(A.aspirants[[i]]))
              } 
      
          # best model/ updates
          best <- which(abc==min(abc))[1]
          abc.option <- abc[best]
          A.option <- A.aspirants[[best]]
          option <- glm.fit(x=x%*%A.option, y=y, weights = weights, 
                         offset = offset, start = as.vector(t(A.option)%*%initials), 
                         family = family, intercept = FALSE)
          if(option$converged==FALSE){       #
            abc.option <- abc.model + 1
            A.option <- A.model
            option <- model
          }
          best.coefs <- which(rowSums(A.option)==1)
          if (length(dim(A)[2]) > 0){ if(dim(A)[2]!=0) {
              excluded <- which(colSums(matrix(A[best.coefs,],ncol=dim(A)[2],byrow=FALSE))>0) # unit == columns of A
              if(length(excluded)>0) {A <- as.matrix(A[,-excluded]) }}}
      } else { abc.option <- abc.option + 1 }    
      
      w <- w+1
      A.models[[w]] <- A.model
      
  } # while

# names A.model, X.model
namen <- c()
for (i in 1:dim(A.model)[2]) {
     f <- A.model[,i]
     f.w <- which(A.model[,i]==1)
     if (sum(f)==1) {namen <- c(namen, colnames(x)[f.w]) } else {
          namen.temp <- colnames(x)[f.w[1]]
          for (j in 2:length(f.w)) {namen.temp <- paste(namen.temp,colnames(x)[f.w[j]],sep="/")}
          namen <- c(namen, namen.temp)
        }
    }      
colnames(A.model) <- namen          

X.model <- x%*%A.model

beta.reduced <- round(model$coefficients, digits=control$accuracy)
names(beta.reduced) <- colnames(A.model)
coefficients <- round(A.model %*% beta.reduced, digits=control$accuracy)
names(coefficients) <- colnames(x)
 
# output
list(coefficients = coefficients, residuals = model$residuals, fitted.values = model$fitted.values,
    effects = model$effects, R = model$R,
    rank = round(model$rank, 2), 
    qr = model$QR,     
    family = model$family,
    linear.predictors = model$linear.predictors, deviance = model$deviance, aic = model$aic,
    null.deviance = model$null.deviance, iter = w, weights = model$weights, prior.weights = model$prior.weights,
    df.residual = model$df.residual, df.null = model$df.null, y = y, converged = model$converged,
    boundary = model$boundary,
    A.model=A.model, X.model=X.model, abc.model=abc.model, A.models=A.models
    )
# return(list(model = model, A.model=A.model, X.model=X.model, abc.model=abc.model, iter=w, A.models=A.models))

}

