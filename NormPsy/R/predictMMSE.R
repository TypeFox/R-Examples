

predictMMSE <- function(model,VarTime,Timelim,nTime,Xprofile,methInteg="GH",nsim=20,draws=FALSE,ndraws=2000)
{
 #dyn.load("backtransformation.so")
 # require(lcmm)

 if(missing(model)) stop("The argument 'model' should be specified")
 #if(missing(newdata)) stop("The argument newdata should be specified")
 if(!inherits(model,"hlme")) stop("Use with hlme object only")
 #if (!inherits(newdata, "data.frame")) stop("newdata should be a data.frame object")
 #if(length(prm.transfo) != (length(nodes)+2)) stop(paste("With",length(nodes),"nodes,",length(nodes)+2,"parameters are required"))
 if(!(methInteg %in% c("GH","MC",0,1))) stop("Only Gauss-Hermite (0 or 'GH') and Monte-Carlo (1 or 'MC') integration methods are available")
 if(methInteg %in% c("GH",0) & !(nsim %in% c(5,7,9,15,20,30,40,50))) stop("For Gauss-Hermite integration method, 'nsim' should be either 5,7,9,15,20,30,40 or 50")
 if(isTRUE(draws) & model$conv!=1) stop("No confidence interval can be computed since the program did not converge properly ")

  
  nodes <- c(0,10,20,23,26,28,30)
  prm.transfo <- c(0,10.33921,8.06659,10.82235,17.69599,11.59447,14.93948,13.16887,13.37304)
  
  nbclasses <- model$ng
  
  #creer newdata
  #if(missing(newdata))
  #{
   time <- seq(Timelim[1],Timelim[2],length.out=nTime) #si 2 varialbes de temps ???
   newdata <- as.data.frame(sapply(Xprofile,rep,length.out=nTime))
   newdata <- cbind(time,newdata)
   colnames(newdata) <- c(VarTime,names(Xprofile))
  #}
  
   mu <- lcmm::predictY(model,newdata,var.time=VarTime)[[1]]
   nobs <- length(mu[,1])
   ch <- matrix(0,sum(model$idea0),sum(model$idea0))
   ch[lower.tri(ch,diag=TRUE)] <- model$cholesky

   if(is.null(model$call$random)) {call_random <- ~-1} else {call_random <- model$call$random }
   z <- all.names(call_random)
   ind_factor <- which(z=="factor")
   olddata <- eval(model$call$data)
   if(length(ind_factor))
   {
    nom.factor <- z[ind_factor+1]
    for (v in nom.factor)
    {
     mod <- levels(as.factor(olddata[,v]))
     if (!all(levels(as.factor(newdata1[,v])) %in% mod)) stop(paste("invalid level in factor", v))
     newdata1[,v] <- factor(newdata1[,v], levels=mod)
    }
   }
   call_random <- gsub("factor","",call_random)

   Z <- model.matrix(formula(paste("~",call_random,sep="")),data=newdata)

   VC <- Z %*% ch %*% t(ch) %*% t(Z) + model$best["stderr"]^2*diag(nobs)


  if(length(model$N)>4)
  {
   if(model$N[5]>0)
   {
     corChar <- as.character(model$call$cor)
     corTime <- newdata[,corChar[2]]

     if(corChar[1]=="BM")
     {
      Corr <- model$best["cor1"]^2 * outer(corTime,corTime,Vectorize(function(x,y) min(abs(x),abs(y))))
     }
     if(corChar[1]=="AR")
     {
      Corr <- model$best["cor2"]^2 * (1/outer(corTime,corTime,Vectorize(function(x,y) exp(abs(x-y)))))^model$best["cor1"]
     }

     VC <- VC + Corr
   }
  }



 nbzitr <- length(nodes)
 methInteg <- ifelse(methInteg %in% c("GH",0),0,1)
 Ymarg <- rep(0,nobs)

 if(methInteg==0)
 {
  VC0 <- diag(VC)
  VC1 <- rep(0,nobs*(nobs+1)/2)
 }
 else
 {
  VC0 <- rep(0,nobs)
  VC1 <- VC[upper.tri(VC,diag=TRUE)]
 }



  if(nbclasses==1)
  {
   if(!isTRUE(draws))
   {
    res <- .Fortran("backtransformation",as.double(mu),as.double(VC0),as.double(VC1),as.integer(nobs),as.double(prm.transfo),
         as.integer(nbzitr),as.double(nodes),as.integer(nsim),as.integer(methInteg),y=as.double(Ymarg),package="NormPsy")

    res$y[which(res$y==9999)] <- NA
    res.pred <-cbind(time, matrix(res$y,ncol=1))
    colnames(res.pred) <- c("time","MMSEpred") 
    class(res.pred) <- "predMMSE" 
    return(res.pred)
   }
   else
   {
    ndraws <- as.integer(ndraws)
    ydraws <- NULL

    Mat <- matrix(0,ncol=length(model$best),nrow=length(model$best))
    # que la partie sup utilis?e donc OK si rien en bas
    Mat[upper.tri(Mat,diag=TRUE)]<- model$V
    Chol <- chol(Mat)
    Chol <- t(Chol)


     doOneDraw <- function()
     {
      bdraw <- rnorm(length(model$best))
      bdraw <- model$best + Chol %*% bdraw
      modeldraw <- model
      modeldraw$best <- bdraw
      mudraw <- lcmm::predictY(modeldraw,newdata,var.time=VarTime)[[1]]

      Ymarg <- rep(0,nobs)
      out <- .Fortran("backtransformation",as.double(mudraw),as.double(VC0),as.double(VC1),as.integer(nobs),as.double(prm.transfo),
         as.integer(nbzitr),as.double(nodes),as.integer(nsim),as.integer(methInteg),y=as.double(Ymarg),package="NormPsy")

      out$y[out$y==9999] <- NA
     
      return(out$y)
     }

    ydraws <- replicate(ndraws,doOneDraw())
    ydistr <- apply(ydraws,1,FUN=function(x) quantile(x[!is.na(x)],probs=c(0.025,0.5,0.975)))
    Ypred_50 <- matrix(ydistr[2,],ncol=1,byrow=F)
    Ypred_2.5 <- matrix(ydistr[1,],ncol=1,byrow=F)
    Ypred_97.5 <- matrix(ydistr[3,],ncol=1,byrow=F)

    Ypred <- cbind(time,Ypred_50,Ypred_2.5,Ypred_97.5)
    colnames(Ypred) <- c("time","MMSEdistr_50%","MMSEdistr_2.5%","MMSEdistr_97.5%")
    class(Ypred) <- "predMMSE"
    return(Ypred)
   }
  }
  else # ie plusieurs classes
  {
   if(!isTRUE(draws))
   {
    res <- NULL
    for(g in 1:nbclasses)
    {
     gmu <- mu[,g]
     if(model$N[4]==0 | g==nbclasses) wg <- 1
     else wg <- (model$best[paste("varprop class",g,sep=" ")])^2

     gVC0 <- wg*VC0
     gVC1 <- wg*VC1

     gres <- .Fortran("backtransformation",as.double(gmu),as.double(gVC0),as.double(gVC1),as.integer(nobs),as.double(prm.transfo),
         as.integer(nbzitr),as.double(nodes),as.integer(nsim),as.integer(methInteg),y=as.double(Ymarg),package="NormPsy")

     gres$y[which(gres$y==9999)] <- NA
     res <- cbind(res,gres$y)
    }
     res <- cbind(time,res) 
    colnames(res) <- c("time",paste("MMSEclass",1:nbclasses,sep=""))
    class(res) <- "predMMSE" 
    return(res)
   }
   else
   {
    res <- NULL
    ndraws <- as.integer(ndraws)
    ydraws <- NULL

    Mat <- matrix(0,ncol=length(model$best),nrow=length(model$best))
    # que la partie sup utilis?e donc OK si rien en bas
    Mat[upper.tri(Mat,diag=TRUE)]<- model$V
    Chol <- chol(Mat)
    Chol <- t(Chol)

    for (j in 1:ndraws)
    {
     bdraw <- rnorm(length(model$best))
     bdraw <- model$best + Chol %*% bdraw
     modeldraw <- model
     modeldraw$best <- bdraw
     mudraw <- lcmm::predictY(modeldraw,newdata,var.time=VarTime)[[1]]

     for(g in 1:nbclasses)
     {
      gmudraw <- mudraw[,g]

      if(model$N[4]==0 | g==nbclasses) wg <- 1
      else wg <- (model$best[paste("varprop class",g,sep=" ")])^2

      gVC0 <- wg*VC0
      gVC1 <- wg*VC1

      Ymarg <- rep(0,nobs)

      gres <- .Fortran("backtransformation",as.double(gmudraw),as.double(gVC0),as.double(gVC1),as.integer(nobs),as.double(prm.transfo),
         as.integer(nbzitr),as.double(nodes),as.integer(nsim),as.integer(methInteg),y=as.double(Ymarg),package="NormPsy")

      gres$y[which(gres$y==9999)] <- NA
      res <- cbind(res,gres$y)
     }
    }
     
    classe <- factor(paste("class",rep(1:nbclasses,ndraws),sep=""))
    ydistr <- apply(res,1,function(v) unlist(tapply(v,classe,FUN=function(x) quantile(x[!is.na(x)],probs=c(0.025,0.5,0.975)))))
    ydistr <- t(ydistr)
    ydistr <- ydistr[,order(rep(c(2,1,3),nbclasses))]
    ydistr <- cbind(time,ydistr) 
    colnames(ydistr) <- c("time",paste("MMSEdistr_50_class",1:model$ng,sep=""),paste("MMSEdistr_2.5_class",1:model$ng,sep=""),paste("MMSEdistr_97.5_class",1:model$ng,sep=""))
    class(ydistr) <- "predMMSE"
    return(ydistr)
   }
  }
}



