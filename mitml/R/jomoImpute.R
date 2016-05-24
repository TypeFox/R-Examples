jomoImpute <- function(data, type, formula, random.L1=c("none","mean","full"),
                       n.burn=5000, n.iter=100, m=10, group=NULL, prior=NULL,
                       seed=NULL, save.pred=FALSE, silent=FALSE){

# wrapper function for the different samplers of the jomo package

  # .model.byFormula <- mitml:::.model.byFormula
  # .model.byType <- mitml:::.model.byType
  # .extractMatrix <- mitml:::.extractMatrix

  # empty objects to assign to
  clname <- yvrs <- y <- ycat <- zcol <- xcol <- pred <- clus <- psave <-
    pvrs <- qvrs <- pnames <- qnames <- NULL

  # *** checks
  if(!missing(type) & !missing(formula)) stop("Only one of 'type' or 'formula' may be specified.")
  if(save.pred & !missing(type)){
    warning("Option 'save.pred' is ignored if 'type' is specified")
    save.pred=FALSE 
  }

  random.L1 <- match.arg(random.L1)

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

  # *** retreive model input
  if(!missing(formula)) .model.byFormula(data, formula, group, group.original,
                                         method="jomo.matrix")
  if(!missing(type)) .model.byType(data, type, group, group.original,
                                   method="jomo.matrix")

  # more checks
  if(any(is.na(pred))) stop("Predictor variables must not contain missing data.")
  if(any(is.na(group))) stop("Grouping variable must not contain missing data.")
  if(any(!sapply(data[yvrs], function(a) is.factor(a) || is.numeric(a)))) stop("Target variables must either be numeric or factors.")
  if(sum(is.na(y))==0 & sum(is.na(ycat))==0) stop("Target variables do not contain any missing data.")

  # reorder colums
  cc <- which(colnames(data) %in% c(clname,grname,yvrs))
  data.ord <- cbind(data[c(clname,grname,yvrs)],data[-cc])

  # *** jomo setup

  ycat.labels <- lapply(data[,colnames(ycat),drop=F],levels)
  # select function
  func <- if(ncol(ycat)==0) "con" else if(ncol(y)==0) "cat" else "mix"
  func <- get(paste0("jomo1ran", func, if(random.L1!="none") "hr", ".MCMCchain"),
              asNamespace("jomo"))

  # data properties per group
  ng <- length(unique(group))
  np <- length(xcol)
  nq <- length(zcol)
  ncon <- ncol(y)
  ncat <- ncol(ycat)
  nr <- ncon + ncat   # combined con + cat (variables)
  y.numcat <- matrix(0,ng,ncat)
  nc <- nr2 <- integer(ng)
  for(gg in unique(group)){
    y.numcat[gg,] <- apply(ycat[group==gg,,drop=F], 2, function(x) length(unique(x[!is.na(x)])))
    nc[gg] <- length(unique(clus[group==gg]))
    nr2[gg] <- ncon+sum(y.numcat[gg,])-length(y.numcat[gg,]) # combined con + cat (indicators)
  }

  # * * * * * * * * * * * * * * * * * * * *

  # seed
  original.seed <- NULL
  if(!is.null(seed)){
    if(exists(".Random.seed", .GlobalEnv)) original.seed <- .Random.seed
    set.seed(seed)
  }

  # priors
  if(is.null(prior)){
    prior <- as.list(unique(group))
    for(gg in unique(group)){
      prior[[gg]] <- list( Binv=diag(1,nr2[gg]),
                           Dinv=diag(1,nq*nr2[gg]) )
      if(random.L1!="none") prior[[gg]]$a <- nr2[gg]
    }
  }

  # prepare output
  ind <- which(is.na(data.ord), arr.ind=TRUE, useNames=FALSE)
  ind <- ind[ ind[,2] %in% which(colnames(data.ord)%in%yvrs),,drop=FALSE ]
  rpm <- matrix(NA, nrow(ind), m)

  bpar <- list(beta=array( NA, c(np,max(nr2),n.burn,ng) ),
               psi=array( NA, c(max(nr2)*nq,max(nr2)*nq,n.burn,ng) ),
               sigma=array( NA, c(ifelse(random.L1=="full",max(nr2)*max(nc),max(nr2)),
                                  max(nr2),n.burn,ng) ))
  ipar <- list(beta=array( NA, c(np,max(nr2),n.iter*m,ng) ),
               psi=array( NA, c(max(nr2)*nq,max(nr2)*nq,n.iter*m,ng) ),
               sigma=array( NA, c(ifelse(random.L1=="full",max(nr2)*max(nc),max(nr2)),
                                  max(nr2),n.iter*m,ng) ))
  
  # burn-in
  if(!silent){
    cat("Running burn-in phase ...\n")
    flush.console()
  }
  glast <- as.list(unique(group))
  for(gg in unique(group)){

    gi <- group==gg
    gprior <- prior[[gg]]

    # jomo1ran...: function arguments (group specific)
    gclus <- clus[gi]
    gclus <- matrix( match(gclus, unique(gclus))-1, ncol=1 )
    func.args <- list( Y=if(ncon>0 & ncat==0) y[gi,,drop=F] else NULL,
                       Y.con=if(ncon>0 & ncat>0) y[gi,,drop=F] else NULL,
                       Y.cat=if(ncat>0) ycat[gi,,drop=F] else NULL,
                       Y.numcat=if(ncat>0) y.numcat[gg,] else NULL,
                       X=pred[gi,xcol,drop=F],
                       Z=pred[gi,zcol,drop=F],
                       clus=gclus,
                       beta.start=matrix(0,np,nr2[gg]),
                       u.start=matrix(0,nc[gg],nq*nr2[gg]),
                       l1cov.start=if(random.L1!="none"){
                         matrix(diag(1,nr2[gg]),nr2[gg]*nc[gg],nr2[gg],byrow=T)
                       }else{
                         diag(1,nr2[gg])
                       },
                       l2cov.start=diag(1,nq*nr2[gg]),
                       start.imp=NULL,
                       l1cov.prior=gprior$Binv,
                       l2cov.prior=gprior$Dinv,
                       a=gprior$a,
                       meth=if(random.L1!="none") "random" else NULL,
                       nburn=n.burn,
                       output=0
    )
    func.args <- func.args[!sapply(func.args,is.null)]
    # imp <- do.call( func, c(func.args,list(nburn=n.burn)) )

    cur <- do.call( func, func.args )
    glast[[gg]] <- cur

    # populate output
    bdim <- dim(cur$collectbeta)[1:2]
    pdim <- dim(cur$collectcovu)[1:2]
    sdim <- dim(cur$collectomega)[1:2]
    bpar[["beta"]][1:bdim[1],1:bdim[2],,gg] <- cur$collectbeta
    bpar[["psi"]][1:pdim[1],1:pdim[2],,gg] <- cur$collectcovu
    if(random.L1=="mean"){
      tmp <- cur$collectomega
      dim(tmp) <- c(nr2[gg],nc[gg],nr2[gg],n.burn)
      bpar[["sigma"]][1:sdim[2],1:sdim[2],,gg] <- apply(tmp,c(1,3,4),mean)
    } else {
      bpar[["sigma"]][1:sdim[1],1:sdim[2],,gg] <- cur$collectomega
    }

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
      gprior <- prior[[gg]]

      # jomo1ran...: function arguments (group specific)
      gclus <- clus[gi]
      gclus <- matrix( match(gclus, unique(gclus))-1, ncol=1 )
      it <- dim(glast[[gg]]$collectbeta)[3]
      func.args <- list( Y=if(ncon>0 & ncat==0) y[gi,,drop=F] else NULL,
                         Y.con=if(ncon>0 & ncat>0) y[gi,,drop=F] else NULL,
                         Y.cat=if(ncat>0) ycat[gi,,drop=F] else NULL,
                         Y.numcat=if(ncat>0) y.numcat[gg,] else NULL,
                         X=pred[gi,xcol,drop=F],
                         Z=pred[gi,zcol,drop=F],
                         clus=gclus,
                         #beta.start=glast[[gg]]$collectbeta[,,it],
                         #u.start=glast[[gg]]$collectu[,,it],
                         #l1cov.start=glast[[gg]]$collectomega[,,it],
                         #l2cov.start=glast[[gg]]$collectcovu[,,it],
                         beta.start=.extractMatrix(glast[[gg]]$collectbeta,it),
                         u.start=.extractMatrix(glast[[gg]]$collectu,it),
                         l1cov.start=.extractMatrix(glast[[gg]]$collectomega,it),
                         l2cov.start=.extractMatrix(glast[[gg]]$collectcovu,it),
                         start.imp=glast[[gg]]$finimp.latnorm,
                         l1cov.prior=gprior$Binv,
                         l2cov.prior=gprior$Dinv,
                         a=gprior$a,
                         meth=if(random.L1!="none") "random" else NULL,
                         nburn=n.iter,
                         output=0
      )
      func.args <- func.args[!sapply(func.args,is.null)]
      # imp <- do.call( func, c(func.args,list(nburn=n.burn)) )

      cur <- do.call( func, func.args )
      glast[[gg]] <- cur

      # populate output
      ri <- (nrow(gclus)+1):nrow(cur$finimp)
      ci <- which(colnames(cur$finimp) %in% yvrs)
      gy.imp[[gg]] <- cur$finimp[ri,ci,drop=F]
      bdim <- dim(cur$collectbeta)[1:2]
      pdim <- dim(cur$collectcovu)[1:2]
      sdim <- dim(cur$collectomega)[1:2]
      iind <- (n.iter*(ii-1)+1):(n.iter*ii)
      ipar[["beta"]][1:bdim[1],1:bdim[2],iind,gg] <- cur$collectbeta
      ipar[["psi"]][1:pdim[1],1:pdim[2],iind,gg] <- cur$collectcovu
      if(random.L1=="mean"){
        tmp <- cur$collectomega
        dim(tmp) <- c(nr2[gg],nc[gg],nr2[gg],n.iter)
        ipar[["sigma"]][1:sdim[2],1:sdim[2],iind,gg] <- apply(tmp,c(1,3,4),mean)
      } else {
        ipar[["sigma"]][1:sdim[1],1:sdim[2],iind,gg] <- cur$collectomega
      }

    }
    y.imp <- data.matrix(do.call(rbind,gy.imp))
    rpm[,ii] <- y.imp[,yvrs][is.na(data.ord[,yvrs,drop=F])]

  }

  if(!silent){
    cat("Done!\n")
  }

  # clean up
  srt <- data.ord[,ncol(data.ord)]
  data.ord <- data.ord[,-ncol(data.ord)]

  # restore seed
  if(is.null(original.seed)){
    rm(".Random.seed", envir = .GlobalEnv)
  }else{
    assign(".Random.seed", original.seed, envir=.GlobalEnv)
  }

  # prepare output data
  if( save.pred & !missing(formula) ) data.ord <- cbind(data.ord,pred[,psave,drop=F])
  attr(data.ord,"sort") <- srt
  attr(data.ord,"group") <- group.original
  if(ncat>0){
    attr(data.ord,"cvrs") <- names(ycat.labels)
    attr(data.ord,"levels") <- y.numcat
    attr(data.ord,"labels") <- ycat.labels
  }
  model <- list(clus=clname, yvrs=yvrs, pvrs=pvrs, qvrs=qvrs)
  attr(model,"full.names") <- list(pvrs=pnames, qvrs=qnames)

  out <- list(
    data=data.ord,
    replacement.mat=rpm,
    index.mat=ind,
    call=match.call(),
    model=model,
    random.L1=random.L1,
    prior=prior,
    iter=list(burn=n.burn, iter=n.iter, m=m),
    par.burnin=bpar,
    par.imputation=ipar
  )
  class(out) <- c("mitml","jomo")
  out
  
}

