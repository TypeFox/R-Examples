panImpute <- function(data, type, formula, n.burn=5000, n.iter=100, m=10,
                      group=NULL, prior=NULL, seed=NULL, save.pred=FALSE,
                      silent=FALSE){

# wrapper function for the Gibbs sampler in the pan package

  # .model.byFormula <- mitml:::.model.byFormula
  # .model.byType <- mitml:::.model.byType

  # empty objects to assign to
  clname <- yvrs <- y <- ycat <- zcol <- xcol <- pred <- clus <- psave <- 
    pvrs <- qvrs <- pnames <- qnames <- NULL

  # *** checks
  if(!missing(type) & !missing(formula)) stop("Only one of 'type' or 'formula' may be specified.")
  if(save.pred & !missing(type)){
    warning("Option 'save.pred' is ignored if 'type' is specified")
    save.pred=FALSE 
  }

  # preserve original order
  if(!is.data.frame(data)) as.data.frame(data)
  data <- cbind(data, original.order=1:nrow(data))
  if(!missing(type)) type <- c(type,0)

  # address additional grouping
  grname <- group
  if(!missing(formula)){
    if(is.null(group)){
      group <- rep(1,nrow(data))
    }else{
      group <- data[,group]
      if(length(group)!=nrow(data)) stop("Argument 'group' is not correctly specified.")
    }
  }
  if(!missing(type)){
    if(sum(type==-1)>1) stop("Argument 'group' is not correctly specified,")
    if(sum(type==-1)==0){
      group <- rep(1,nrow(data))
    }else{
      group <- data[,type==-1]
      type[type==-1] <- 0
    }
  }
  group.original <- group
  group <- as.numeric(factor(group,levels=unique(group)))

  # *** model input
  if(!missing(formula)) .model.byFormula(data, formula, group, group.original,
                                         method="pan")
  if(!missing(type)) .model.byType(data, type, group, group.original,
                                   method="pan")

  # reorder colums
  cc <- which(colnames(data) %in% c(clname,grname,yvrs))
  data.ord <- cbind(data[c(clname,grname,yvrs)],data[-cc])

  # * * * * * * * * * * * * * * * * * * * *

  if(sum(is.na(y))==0) stop("Target variables do not contain any missing data.")

  if(is.null(prior)){
    prior <- list( a=ncol(y), Binv=diag(1,ncol(y)),
      c=ncol(y)*length(zcol), Dinv=diag(1,ncol(y)*length(zcol)) )
  }

  if(is.null(seed)){
    set.seed(as.integer(runif(1,0,10^6)))
  }else{
    set.seed(as.integer(seed))
  }
  rns <- sapply(unique(group), function(x,m) as.integer(runif(m+1,0,10^6)), m=m)

  # prepare output
  ind <- which(is.na(data.ord), arr.ind=TRUE, useNames=FALSE)
  ind <- ind[ ind[,2] %in% which(colnames(data.ord)%in%colnames(y)),,drop=FALSE ]
  rpm <- matrix(NA, nrow(ind), m)

  ng <- length(unique(group))
  np <- length(xcol)
  nq <- length(zcol)
  nr <- ncol(y)
  bpar <- list(beta=array( NA, c(np,nr,n.burn,ng) ),
               psi=array( NA, c(nr*nq,nr*nq,n.burn,ng) ),
               sigma=array( NA, c(nr,nr,n.burn,ng) ))
  ipar <- list(beta=array( NA, c(np,nr,n.iter*m,ng) ),
               psi=array( NA, c(nr*nq,nr*nq,n.iter*m,ng) ),
               sigma=array( NA, c(nr,nr,n.iter*m,ng) ))
  
  # burn-in
  if(!silent){
    cat("Running burn-in phase ...\n")
    flush.console()
  }
  glast <- as.list(unique(group))
  for(gg in unique(group)){

    gi <- group==gg
    gy <- y[gi,]
    gpred <- pred[gi,]
    gclus <- clus[gi]
    # sort 1, ..., k
    gclus <- match(gclus, unique(gclus)) 

    cur <- pan::pan(gy, subj=gclus, gpred, xcol, zcol, prior, seed=rns[1,gg], iter=n.burn)
    glast[[gg]] <- cur$last

    bpar[["beta"]][,,,gg] <- cur$beta
    bpar[["psi"]][,,,gg] <- cur$psi
    bpar[["sigma"]][,,,gg] <- cur$sigma

  }
  
  # imputation
  for(ii in 1:m){
    if(!silent){
      cat("Creating imputed data set (",ii,"/",m,") ...\n")
      flush.console()
    }

    gy.imp <- as.list(unique(group))
    for(gg in unique(group)){

      gi <- group==gg
      gy <- y[gi,]
      gpred <- pred[gi,]
      gclus <- clus[gi]
      # sort 1, ..., k
      gclus <- match(gclus, unique(gclus)) 
  
      cur <- pan::pan(gy, subj=gclus, gpred, xcol, zcol, prior, seed=rns[ii+1,gg], iter=n.iter, 
        start=glast[[gg]])
      glast[[gg]] <- cur$last
  
      # populate output
      gy.imp[[gg]] <- cur$y
      ipar[["beta"]][,,(n.iter*(ii-1)+1):(n.iter*ii),gg] <- cur$beta
      ipar[["psi"]][,,(n.iter*(ii-1)+1):(n.iter*ii),gg] <- cur$psi
      ipar[["sigma"]][,,(n.iter*(ii-1)+1):(n.iter*ii),gg] <- cur$sigma

    }
    y.imp <- do.call(rbind,gy.imp)
    rpm[,ii] <- y.imp[is.na(y)]

  }

  if(!silent){
    cat("Done!\n")
  }

  # clean up
  srt <- data.ord[,ncol(data.ord)]
  data.ord <- data.ord[,-ncol(data.ord)]

  # prepare output data
  if( save.pred & !missing(formula) ) data.ord <- cbind(data.ord,pred[,psave,drop=F])
  attr(data.ord,"sort") <- srt
  attr(data.ord,"group") <- group.original
  model <- list(clus=clname, yvrs=yvrs, pvrs=pvrs, qvrs=qvrs)
  attr(model,"full.names") <- list(pvrs=pnames, qvrs=qnames)

  out <- list(
    data=data.ord,
    replacement.mat=rpm,
    index.mat=ind,
    call=match.call(),
    model=model,
    random.L1="none",
    prior=prior,
    iter=list(burn=n.burn, iter=n.iter, m=m),
    par.burnin=bpar,
    par.imputation=ipar
  )
  class(out) <- c("mitml","pan")
  out
  
}

