Classify <-
function(x,y,xte=NULL,rho=0,beta=1,rhos=NULL,type=c("mle","deseq","quantile"), prior=NULL, transform=TRUE, alpha=NULL){
  if(is.null(xte)){
    xte <- x
    warning("Since no xte was provided, testing was performed on training data set.")
  }
  if(!is.null(rho) && length(rho)>1) stop("Can only enter 1 value of rho. If you would like to enter multiple values, use rhos argument.")
  type <- match.arg(type)
  if(!transform && !is.null(alpha)) stop("You have asked for NO transformation but have entered alpha.")
  if(transform && is.null(alpha)) alpha <- FindBestTransform(x)
  if(transform){
    if(alpha<=0 || alpha>1) stop("alpha must be between 0 and 1")
    x <- x^alpha
    xte <- xte^alpha
  }  
  if(is.null(prior)) prior <- rep(1/length(unique(y)), length(unique(y)))
  if(is.null(rho)&&is.null(rhos)) stop("Must enter rho or rhos.")
  null.out <- NullModel(x, type=type)
  ns <- null.out$n
  nste <- NullModelTest(null.out,x,xte,type=type)$nste
  uniq <- sort(unique(y))
  if(is.null(rhos)){
    ds <- GetD(ns,x,y,rho,beta)
    discriminant <- matrix(NA, nrow=nrow(xte), ncol=length(uniq))
    for(k in 1:length(uniq)){
      discriminant[,k] <- rowSums(scale(xte,center=FALSE,scale=(1/log(ds[k,])))) - rowSums(scale(nste,center=FALSE,scale=(1/ds[k,]))) + log(prior[k])
    }
    save <- list(ns=ns,nste=nste,ds=ds,discriminant=discriminant,ytehat=uniq[apply(discriminant,1,which.max)], alpha=alpha, rho=rho,x=x,y=y,xte=xte,alpha=alpha,type=type)
    class(save) <- "poicla"
    return(save)
  } else {
    save <- list()
    ds.list <- GetD(ns,x,y,rho=NULL, rhos=rhos,beta)
    for(rho in rhos){
      ds <- ds.list[[which(rhos==rho)]]
      discriminant <- matrix(NA, nrow=nrow(xte), ncol=length(uniq))
      for(k in 1:length(uniq)){
        discriminant[,k] <- rowSums(scale(xte,center=FALSE,scale=(1/log(ds[k,])))) - rowSums(scale(nste,center=FALSE,scale=(1/ds[k,]))) + log(prior[k])
      }
      save[[which(rhos==rho)]] <- (list(ns=ns,nste=nste,ds=ds,discriminant=discriminant,ytehat=uniq[apply(discriminant,1,which.max)], alpha=alpha, rho=rho,x=x,y=y,xte=xte,alpha=alpha,type=type))
    }
    class(save) <- "poicla"
    return(save)
  }
}

