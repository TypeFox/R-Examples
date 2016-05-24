###{{{ missingModel

missingModel <- function(model,data,var=endogenous(model),fix=FALSE,type=2,keep=NULL,weight=NULL,weight2=NULL,cluster=NULL,...) {
  if (!inherits(model,"lvm")) stop("Needs a lvm-object")
  if (type==3) {
    var <- manifest(model)
  }

  data.mis <- is.na(data[,var,drop=FALSE])
  colnames(data.mis) <- var
  patterns <- unique(data.mis,MARGIN=1)

  mis.type <- apply(data.mis,1,
                  function(x) which(apply(patterns,1,function(y) identical(x,y))))
  pattern.allmis <- which(apply(patterns,1,all)) ## Remove entry with all missing

  models <- datasets <- weights <- weight2 <- clusters <- c()
  mymodel <- baptize(model)
  pattern.compl <- 0
  count <- 0
  A <- index(model)$A
  topendo <- endogenous(model,top=TRUE)
  exo <- exogenous(model)
  exclude <- c()

  warned <- FALSE
  for (i in setdiff(seq_len(nrow(patterns)),pattern.allmis)) {
    exoremove <- c()
    includemodel <- TRUE
    count <- count+1
    mypattern <- patterns[i,]
    m0 <- mymodel;
    if (any(mypattern)) {
      latent(m0,zero=FALSE) <- colnames(data.mis)[mypattern]
      if (type>1) {
        mytop <- intersect(topendo,colnames(data.mis)[mypattern])
        if (!is.null(mytop)) {
          kill(m0) <- mytop     
          for (xx in exo) {
          ## If exogenous variable only have effect on missing variables,
          ##  then remove it from the model
              if (all(c(rownames(A)[A[xx,]==1])%in%mytop) &&
                  !(xx%in%m0$par)
                  ##&& !(xx%in%names(index(m0))$parval)
                  ) {
                  exoremove <- c(exoremove,xx)
                  kill(m0) <- xx
              }
          }
        }
      }
    } else
    pattern.compl <- count
##    d0 <- data[mis.type==i,manifest(m0),drop=FALSE];
    d0 <- data[which(mis.type==i),c(manifest(m0),keep),drop=FALSE];
    w0.var <- intersect(manifest(m0),colnames(weight))
    w0 <- weight[which(mis.type==i),w0.var,drop=FALSE];
    if (!is.list(weight2)) {
      w02.var <- intersect(manifest(m0),colnames(weight2))
      w02 <- weight2[which(mis.type==i),w02.var,drop=FALSE];
    } else {
      weight2 <- weight2
    }

    clust0 <- cluster[which(mis.type==i)]
    ex0 <- exogenous(m0) <- setdiff(exo,exoremove)
    xmis <- which(apply(d0[,ex0,drop=FALSE],1,function(x) any(is.na(x))))

    if (length(xmis)>0) {
      misx <- ex0[apply(d0[xmis,ex0,drop=FALSE],2,function(x) any(is.na(x)))]
      if (!warned)
          warning("Missing exogenous variables: ", paste(misx,collapse=","),
                  ". Removing rows...")
      warned <- TRUE
      d0 <- d0[-xmis,,drop=FALSE]
      w0 <- w0[-xmis,,drop=FALSE]
      clust0 <- clust0[-xmis]
      w02 <- w02[-xmis,,drop=FALSE]
    }
    if (length(misx <- intersect(ex0,latent(m0)))>0) {
      warning("Missing exogenous variables:", paste(misx,collapse=","),
              "! Remove manually!.")
    }
##    else
    {
      if( sum(unlist(index(m0)[c("npar","npar.mean")]))>0 ) {
        models <- c(models, list(m0))
        datasets <- c(datasets, list(d0))
        weights <- c(weights, list(w0))
        if (!is.list(weight2))
          weight2 <- c(weight2, list(w02))
        clusters <- c(clusters, list(clust0))
      } else {
        exclude <- c(exclude,count)
      }
    }
  }

  rmset <- c()
  for (i in seq_len(length(datasets))) {
    if (nrow(datasets[[i]])==0) rmset <- c(rmset,i)
  }
  if (length(rmset)>0) {
    models[[rmset]] <- NULL
    datasets[[rmset]] <- NULL
    weights[[rmset]] <- NULL
    weight2[[rmset]] <- NULL
    clusters[[rmset]] <- NULL
    patterns <- patterns[-rmset,,drop=FALSE]
  }

  Patterns <- patterns
  if (length(exclude)>0)
    Patterns <- Patterns[-exclude,]
  pattern.allcomp<- which(apply(Patterns,1,function(x) all(!x))) ## Complete cases

  res <- list(models=models, datasets=datasets,
              weights=weights,
              weight2=weight2,
              clusters=clusters,
              patterns=Patterns,
              pattern.compl=pattern.compl,
              pattern.allmis=pattern.allmis,
              pattern.allcomp=pattern.allcomp,
              mis.type=mis.type)
  return(res)
}

###}}}

###{{{ estimate.MAR.lvm

##' @export
estimate.MAR <- function(x,data,which=endogenous(x),fix,type=2,startcc=FALSE,control=list(),silent=FALSE,weight,weight2,cluster,onlymodel=FALSE,estimator="gaussian",hessian=TRUE,keep=NULL,...) {
  cl <- match.call()

  Debug("estimate.MAR")
  redvar <- intersect(intersect(parlabels(x),latent(x)),colnames(data))
  if (length(redvar)>0 & !silent)
    warning(paste("Remove latent variable colnames from dataset",redvar))

  xfix <- setdiff(colnames(data)[(colnames(data)%in%parlabels(x,exo=TRUE))],latent(x))
  if (missing(fix))
    fix <- ifelse(length(xfix)>0,FALSE,TRUE)

  S <- diag(nrow=length(manifest(x)));
  mu <- rep(0,nrow(S));
  K <- length(exogenous(x))
  vnames <- index(x)$manifest
  names(mu) <- rownames(S) <- colnames(S) <- vnames
  if (K>0) {
    xx <- subset(Model(x),exogenous(x))
    exogenous(xx) <- NULL
    covfix(xx, vars(xx)) <- NA
    xx <- covariance(xx,exogenous(x),exogenous(x))
    datax <- data[,exogenous(x),drop=FALSE]
    exo.idx <- match(exogenous(x),manifest(x))

    mu0 <- colMeans(datax,na.rm=TRUE)
    cov0 <- cov(datax,use="pairwise.complete.obs")*(nrow(datax)-1)/nrow(datax)
    cov0upper <- cov0[upper.tri(cov0,diag=TRUE)]
    exogenous(xx) <- NULL
    coefpos <- matrices(xx,seq_len(K*(K-1)/2+K))$P
    ii <- coefpos[upper.tri(coefpos,diag=TRUE)]
    start <- c(mu0, cov0upper[order(ii)])
    S[exo.idx,exo.idx] <- cov0
    mu[exo.idx] <- mu0
    ##    message("\n")
  }

  x0 <- x
  x <- fixsome(x, measurement.fix=fix, exo.fix=TRUE, S=S, mu=mu, n=1)
  if (!silent)
    message("Identifying missing patterns...")

  val <- missingModel(x,data,var=which,type=type,keep=c(keep,xfix),weight=weight,weight2=weight2,cluster=cluster,...)
  if (!silent)
    message("\n")

  if (nrow(val$patterns)==1) {
    res <- estimate(x,data=data,fix=fix,weight=weight,weight2=weight2,estimator=estimator,silent=silent,control=control,...)
    return(res)
  }

  if (startcc & is.null(control$start)) {
    if (!silent)
      message("Obtaining starting value...")
    start0 <- rep(1,sum(unlist(index(x)[c("npar","npar.mean")])))
    mystart <- tryCatch(
                        (estimate(x,data=na.omit(data),silent=TRUE,
                                     weight=weight,weight2=weight2,estimator=estimator,quick=TRUE,...
                                      )),
                        error=function(e) rep(1,sum(unlist(index(x)[c("npar","npar.mean")])))
                        )
    control$start <- mystart
    if (!silent)
      message("\n")
  }
  if (is.null(control$meanstructure))
    control$meanstructure <- TRUE

  mg0 <- with(val, suppressWarnings(multigroup(models,datasets,fix=FALSE,exo.fix=FALSE,missing=FALSE)))
  if (!is.null(names(control$start))) {
    parorder1 <- attributes(parpos(mg0,p=names(control$start)))$name
    paridx <- match(parorder1,names(control$start))
    newpos <- paridx[which(!is.na(paridx))]
    start0 <- control$start
    start0[which(!is.na(paridx))] <- control$start[na.omit(paridx)]
    names(start0)[which(!is.na(paridx))] <- names(control$start[na.omit(paridx)])
    control$start <- start0
  }


  if (onlymodel) return(list(mg=mg0,val=val,weight=val$weights,weight2=val$weight2,cluster=val$clusters))

  if (all(unlist(lapply(val$weights,is.null)))) val$weights <- NULL
  if (all(unlist(lapply(val$weight2,is.null)))) val$weight2 <- NULL
  if (all(unlist(lapply(val$clusters,is.null)))) val$clusters <- NULL

  e.mis <- estimate(mg0,control=control,silent=silent,
                    weight=val$weights,weight2=val$weight2,
                    cluster=val$clusters,estimator=estimator,...)

  cc <- coef(e.mis,level=1)
  mynames <- c()
  if (e.mis$model$npar.mean>0)
    mynames <- c(mynames,paste0("m",seq_len(e.mis$model$npar.mean)))
   if (e.mis$model$npar>0)
     mynames <- c(mynames,paste0("p",seq_len(e.mis$model$npar)))
  rownames(cc) <- mynames


  mycc <- val$pattern.allcomp ## Position of complete-case model
  nmis <- with(val, as.numeric(table(mis.type)[pattern.allmis])) ## Number of completely missing observations
  if (length(nmis)>0 & length(mycc)>0) ## Any individuals with all missing?
    if (val$pattern.allmis<mycc)
      mycc <- mycc-1

  if (length(xfix)>0) {
    nrow <- length(vars(x))
    xpos <- lapply(xfix,function(y) which(regfix(x)$labels==y))
    colpos <- lapply(xpos, function(y) ceiling(y/nrow))
    rowpos <- lapply(xpos, function(y) (y-1)%%nrow+1)
    myfix <- list(var=xfix, col=colpos, row=rowpos)
    for (i in seq_along(xfix))
      regfix(x, from=vars(x)[rowpos[[i]]],to=vars(x)[colpos[[i]]]) <-
        rep(colMeans(data[,xfix[i],drop=FALSE],na.rm=TRUE),length(rowpos[[i]]))
    x <-
      updatelvm(x,zeroones=TRUE,deriv=TRUE)
  }

  ord <- c()
  ordlist <- list()
  for (i in seq_len(nrow(val$patterns))) {
    ordlist <- c(ordlist, list(which(val$mis.type==i)))
    ord <- c(ord, ordlist[[i]])
  }

  res <- with(val, list(coef=cc,
                        patterns=patterns, table=table(mis.type),
                        mis.type=mis.type,
                        order=ord,
                        orderlist=ordlist,
                        nmis=nmis,
                        allmis=pattern.allmis,
                        cc=mycc,
                        ncc=as.numeric(table(mis.type)[pattern.allcomp]),
                        multigroup=e.mis$model,
                        estimate=e.mis,
                        model=x,
                        model0=x0,
                        vcov=e.mis$vcov, opt=e.mis$opt,
                        control=control,
                        data=list(model.frame=data),
                        estimator=estimator,
                        call=cl
                        ))
  class(res) <- c("lvm.missing","lvmfit")
  if (inherits(e.mis,"lvmfit.randomslope"))
    class(res) <- c(class(res),"lvmfit.randomslope")

  if (hessian & is.null(cluster)) {
    if (!silent)
      message("Calculating asymptotic variance...\n")
    res$vcov <- solve(information(res$estimate,type="hessian"))
    cc[] <- coef(e.mis,level=1,vcov=res$vcov)
    res$coef <- cc
  }

  return(res)
}

###}}} estimate.MAR.lvm
