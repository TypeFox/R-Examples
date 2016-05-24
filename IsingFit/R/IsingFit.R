IsingFit <-
  function(x, family='binomial', AND = TRUE, gamma = 0.25, plot = TRUE, progressbar = TRUE, lowerbound.lambda = NA,...){
    t0 <- Sys.time()
    xx <- x
    if (family!='binomial') 
      stop ("This procedure is currently only supported for binary (family='binomial') data")
    
    NodesToAnalyze <- apply(x,2,sd, na.rm=TRUE) != 0
    names(NodesToAnalyze) <- colnames(x)
    if (!any(NodesToAnalyze)) stop("No variance in dataset")
    if (any(!NodesToAnalyze))
    {
      warning(paste("Nodes without variance:",paste(colnames(x)[!NodesToAnalyze],collapse = ", ")))
    }
    
    x <- as.matrix(x)
    allthemeans <- colMeans(x)
    x <- x[,NodesToAnalyze,drop=FALSE]
    nvar <- ncol(x)
    p <- nvar - 1
    intercepts <- betas <- lambdas <- list(vector,nvar)
    nlambdas <- rep(0,nvar)
    for (i in 1: nvar){
      a <- glmnet(x[,-i], x[,i], family = family)
      intercepts[[i]] <- a$a0
      betas[[i]] <- a$beta
      lambdas[[i]] <- a$lambda
      nlambdas[i] <- length(lambdas[[i]])
    }
    
    if (progressbar==TRUE) pb <- txtProgressBar(max=nvar, style = 3)
    P <- logl <- sumlogl <- J <- matrix(0, max(nlambdas), nvar)
    for (i in 1:nvar)
    {
      J[1:ncol(betas[[i]]),i] <- colSums(betas[[i]]!=0)
    }
    logl_M <- P_M <- array(0, dim=c(nrow(x),max(nlambdas), nvar) )
    N <- nrow(x)
    for (i in 1:nvar){  # i <- 1
      betas.ii <- as.matrix( betas[[i]] )
      int.ii <- intercepts[[i]]
      y <- matrix( 0 , nrow=N , ncol= ncol(betas.ii) ) 
      xi <- x[,-i]
      NB <- nrow( betas.ii) # number of rows in beta
      for (bb in 1:NB){   # bb <- 1
        y <- y + betas.ii[rep(bb,N),] * xi[,bb]
      }
      y <- matrix( int.ii , nrow=N , ncol=ncol(y) , byrow=TRUE ) + y
      # number of NAs
      n_NA <- max(nlambdas)-ncol(y)
      if (n_NA > 0 ){ 
        for ( vv in 1:n_NA){ 
          y <- cbind( y , NA ) 
        } 
      }
      # calculate P matrix
      P_M[,,i] <- exp(y*x[,i])/(1+exp(y))
      logl_M[,,i] <- log(P_M[,,i])  
      if (progressbar==TRUE) setTxtProgressBar(pb, i)
    }
    
    logl_Msum <- colSums( logl_M , 1, na.rm=FALSE )
    if (progressbar==TRUE) close(pb)
    sumlogl <- logl_Msum 
    sumlogl[sumlogl==0]=NA
    penalty <- J * log(nrow(x)) + 2 * gamma * J * log(p)
    EBIC <- -2 * sumlogl + penalty
    
    lambda.mat <- matrix(NA,nrow(EBIC),ncol(EBIC))
    for (i in 1:nvar){
      lambda.mat[,i] <- c(lambdas[[i]],rep(NA,nrow(EBIC)-length(lambdas[[i]])))
    }
    
    if(!is.na(lowerbound.lambda)){
      EBIC <- EBIC/(lambda.mat>=lowerbound.lambda)*1
    }
    
    lambda.opt <- apply(EBIC,2,which.min)
    lambda.val <- rep(NA,nvar)
    thresholds <- 0
    for(i in 1:length(lambda.opt)){
      lambda.val[i] <- lambda.mat[lambda.opt[i],i]
      thresholds[i] <- intercepts[[i]][lambda.opt[i]]
    }
    weights.opt <- matrix(,nvar,nvar)
    for (i in 1:nvar){
      weights.opt[i,-i] <- betas[[i]][,lambda.opt[i]]
    }
    asymm.weights <- weights.opt
    diag(asymm.weights)=0
    if (AND==TRUE) {
      adj <- weights.opt
      adj <- (adj!=0)*1
      EN.weights <- adj * t(adj)
      EN.weights <- EN.weights * weights.opt
      meanweights.opt <- (EN.weights+t(EN.weights))/2
      diag(meanweights.opt) <- 0 
    } else {
      meanweights.opt <- (weights.opt+t(weights.opt))/2
      diag(meanweights.opt) <- 0
    }
    graphNew <- matrix(0,length(NodesToAnalyze),length(NodesToAnalyze))
    graphNew[NodesToAnalyze,NodesToAnalyze] <- meanweights.opt
    colnames(graphNew) <- rownames(graphNew) <- colnames(xx)
    threshNew <- ifelse(allthemeans > 0.5, -Inf, Inf)
    threshNew[NodesToAnalyze] <- thresholds
    if (plot==TRUE) notplot=FALSE else notplot=TRUE
    q <- qgraph(graphNew,layout='spring',labels=names(NodesToAnalyze),DoNotPlot=notplot,...)
    Res <- list(weiadj = graphNew, thresholds = threshNew, q = q, gamma = gamma, 
                AND = AND, time = Sys.time() - t0, asymm.weights = asymm.weights,
                lambda.values = lambda.val)
    class(Res) <- "IsingFit"
    return(Res)
  }

## Methods:
plot.IsingFit <- function(object,...) qgraph(object$q,DoNotPlot = FALSE, ...)

print.IsingFit <- function(x)
{
  cat("Estimated network:\n")
  
  print(round(x$weiadj,2))
  
  cat("\n\nEstimated Thresholds:\n")
  
  print(x$thresholds)  
}

summary.IsingFit <- function(object)
{
  cat("\tNetwork Density:\t\t", round(mean(object$weiadj[upper.tri(object$weiadj)]!=0),2),"\n",
      "Gamma:\t\t\t",round(object$gamma,2),"\n",
      "Rule used:\t\t",ifelse(object$AND,"And-rule","Or-rule"),"\n",
      "Analysis took:\t\t",format(object$time,format="%s"),"\n"
  )
}

exportNetLogo <- function(object,objectname,....)
{
  if (is.character(object))
  {
    object <- paste0('"',object,'"')
  }
  if (is.vector(object))
  {
    res=paste(paste0("[",paste(object, collapse = " "),"]"),collapse="\n")
  } else if (is.matrix(object))
  {
    res=paste("[\n",paste(apply(object,1,function(s)paste0("[",paste(s, collapse = " "),"]")),collapse="\n"),"\n]")
  } else stop("Object not supported")
  write.table(res,file=paste0(objectname,".txt"),row.names=FALSE, col.names=FALSE,quote=FALSE)
  return(res)
}


# 
# print(fit)
# plot(fit)
# summary(fit)
# export(fit)

