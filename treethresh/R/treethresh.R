######################### TIM

subtree <- function(object,C)
  UseMethod("subtree")

prune <- function(object,v=5,sd.mult=0.5,plot=TRUE)
  UseMethod("prune")

get.w <- function(object,C)
  UseMethod("get.w")

get.t <- function(object,C)
  UseMethod("get.t")

thresh <- function(object,data,C,postmed=TRUE)
  UseMethod("thresh")

treethresh <- function(data,beta,criterion="score",control=list(),rho=sys.frame(sys.parent())) {
  use.beta <- FALSE
  use.data <- NULL
  if (!missing(beta) && !is.null(beta)) {
    use.beta <- TRUE
    use.data <- as.array(beta)
  } else {
    if (!missing(data) && !is.null(data)) {
      use.data <- as.array(data)
    }
  }
  if (is.null(use.data))
    stop("You either need to specify the original data sequence or the betas.")
  default.control <- list(max.iter=30, max.depth=10, minimum.width=3, minimum.size=5^length(dim(use.data)), tolerance.grad=1e-8, tolerance=1e-6, absolute.improvement=-Inf, relative.improvement=-Inf, absolute.criterion=0, beta.max=1e5, a=0.5, crit.signif=0.90, lr.signif=0.5, first.step=0.2)
  if (!is.list(control))
    stop("The control argument must be a list.")
  for (name in names(default.control)) {
    if (!(name %in% names(control)))
      control[[name]] <- default.control[[name]]
  }
  list.diff <- setdiff(names(control),names(default.control))
  if (length(list.diff)>0) {
    warning.text <- "The following arguments in the control list were ignored:"
    for (elt in list.diff)
      warning.text <- paste(warning.text," \"",elt,"\"",sep="")
    warning(warning.text)
  }
  if (!is.function(criterion)) {
    possible.criteria <- c("score","likelihood","heuristic")
    criterion <- possible.criteria[pmatch(criterion,possible.criteria)]
  }
  data.dim <- dim(use.data)
  result<-.External("fit_tree",data=as.numeric(use.data),dims=as.integer(length(data.dim)),size=as.integer(data.dim),use.beta=as.logical(use.beta),criterion=criterion,control=control,rho=rho,PACKAGE="treethresh")
  names(result) <- c("membership","splits","beta")
  result$membership <- array(result$membership,dim=data.dim)
  result$beta <- array(result$beta,dim=data.dim)
  dimnames(result$splits)=list(NULL,c("id","parent.id","dim","pos","left.child.id","right.child.id","crit","w","t","loglikelihood","alpha","C"))
  if (!all(is.na(result$splits[,"alpha"])))
    result$splits[,"C"] <- result$splits[,"alpha"]/max(result$splits[,"alpha"],na.rm=TRUE)
  else
    result$splits[,"C"] <- NA
  result$splits[,"t"] <- tfromw(result$splits[,"w"])
  result$control <- control
  result$criterion <- criterion
  if (!missing(data))
    result$data <- data
  class(result) <- "treethresh"
  result
}

subtree.treethresh <- function(object,C) {
  if (missing(C) || is.null(C))
    return(object)
  if (all(is.na(object$splits[,"C"])))
    return(object)
  delete.me <- rep(FALSE,nrow(object$splits))
  for (i in nrow(object$splits):2) {
    if (object$splits[i,"C"]<C) {
      delete.me[i] <- TRUE
      parent.id <- abs(object$splits[i,"parent.id"])
      object$splits[parent.id==object$splits[,"id"],3:7] <- NA
      object$membership[object$membership==object$splits[i,"id"]] <- parent.id
    }
  }
  object$splits <- object$splits[!delete.me,,drop=FALSE]
  object
}

get.w.treethresh <- function(object,C) {
  if (!missing(C))
    if (!is.null(C))
      object <- subtree.treethresh(object,C)
  map <- numeric(max(object$splits[,"id"]))
  map[object$splits[,"id"]] <- 1:nrow(object$splits)
  array(object$splits[map[as.vector(object$membership)],"w"],dim=dim(object$membership))
}

get.t.treethresh <- function(object,C) {
  if (!missing(C))
    if (!is.null(C))
      object <- subtree.treethresh(object,C)
  map <- numeric(max(object$splits[,"id"]))
  map[object$splits[,"id"]] <- 1:nrow(object$splits)
  array(dim=dim(object$membership),object$splits[map[as.vector(object$membership)],"t"])
}




# Compute best C in cross-validation
prune.treethresh <- function(object,v=5,sd.mult=0.5,plot=TRUE) {
  sizes <- function(object) {
    C <- object$splits[-1,"C"]
    C.values <- sort(unique(C))
    c(apply(outer(C,C.values,function(a,b) as.numeric(a>=b)),2,sum)/2+1,1)
  }
  if (is.na(v))
    return(object)
  if (all(is.na(object$splits[,"C"])))
    return(object)
  # Split the dataset into v portions
  permut <- sample(length(object$data))
  sample.size <- length(object$data)/v
  C.values <- sort(as.numeric(na.omit(c(unique(object$splits[,"C"]),2))))
  results <- matrix(nrow=v,ncol=length(C.values))
  to <- 0
  # Fitting v cv models
  for (j in 1:v) {
    from <- to + 1
    to <- round(sample.size*j)
    test.index <- permut[from:to]
    a.blanked <- object$data[test.index]
    if (is.function(object$criterion) || (object$criterion=="heuristic")) {
      a.data <- object$data
      a.data[test.index] <- NA
      a.tree <- treethresh(data=a.data,control=object$control,criterion=object$criterion)
    } else {
      a.beta <- object$beta
      a.beta[test.index] <- NA
      a.tree <- treethresh(beta=a.beta,control=object$control,criterion=object$criterion)
    }
    for (k in 1:length(C.values)) {
      a.tree <- subtree.treethresh(a.tree,C=C.values[k])
      w <- get.w.treethresh(a.tree,C=NULL)[test.index]
      beta <- pmin(1e5,object$beta[test.index])
      results[j,k] <- sum(log(1+w*beta))
    }
  }
  results <- results - results[,ncol(results)] 
  results.mean <- apply(results,2,mean)
  results.sd <- apply(results,2,sd)
  lower.bound <- (results.mean-sd.mult*results.sd)[which.max(results.mean)]
  if (plot) {
    results.range <- c(min(results.mean-results.sd),max(results.mean+results.sd))
    ns <- 1:length(results.mean)
    plot(ns,results.mean,xaxt="n",ylab="Loglikelihood from CV (+const)",type="o",xlab="C",ylim=results.range)
    segments(ns, results.mean - results.sd, ns, results.mean + results.sd)
    axis(1, at = ns, labels = c(format(C.values[-length(C.values)],digits=3),">1"))
    mtext("number of regions", side = 3, line = 3)
    axis(3, at = ns, labels = as.character(sizes(object)))
    abline(h=lower.bound,lty=3)
  }
  subtree.treethresh(object,C.values[max(which(results.mean>=lower.bound))])
}

thresh.treethresh <- function(object,data,C,postmed=TRUE) {
  if (missing(data) && !is.null(object$data))
    data <- object$data
  data.dim <- dim(data)
  if (postmed) {
    w <- get.w.treethresh(object,C)
    if (length(data)!=length(w))
      stop("Data does not have the right length.")
    data <- postmed.laplace(data,w)
  } else {
    t <- get.t.treethresh(object,C)
    if (length(data)!=length(t))
      stop("Data does not have the right length.")
    data[abs(data)<t] <- 0
  }
  dim(data) <- data.dim
  data
}

######################### TIM_WAVE

# Input for data or betas is a list 
wtthresh <- function(data,beta,weights,control=list()) {
  use.beta <- FALSE
  src.data <- NULL
  if (!missing(beta) && !is.null(beta)) {
    use.beta <- TRUE
    src.data <- beta
    data <- NULL
  } else {
    if (!missing(data) && !is.null(data)) {
      src.data <- data
    }
  }
  if (missing(weights))
    weights <- rep(1,length(src.data))
  if (length(weights)!=length(src.data))
    stop("There must be a weight for every array in the data list.")
  if (is.null(src.data))
    stop("You either need to specify the original data sequence or the betas.")
  use.data <- list()
  data.dim <- length(dim(src.data[[1]]))
  for (i in 1:length(src.data)) {
    use.data[[i]] <- as.array(src.data[[i]])
    cur.dim <- dim(use.data[[i]])
    if (length(unique(cur.dim))!=1)
      stop("Input data needs to be a list of squared arrays.")
    if (data.dim!=length(cur.dim))
      stop("Input data must have identical dimensions.")
  }
  default.control <- list(max.iter=30, max.depth=10, minimum.width=4, min.minimum.width=1, minimum.size=8^data.dim, min.minimum.size=4^data.dim,tolerance.grad=1e-8, tolerance=1e-6, absolute.improvement=-Inf, relative.improvement=-Inf, absolute.criterion=0, beta.max=1e5, a=0.5, rescale.quantile=0.5, lr.signif=0.5, first.step=0.2, min.width.scale.factor=1, min.size.scale.factor=1)
  if (!is.list(control))
    stop("The control argument must be a list.")
  for (name in names(default.control)) {
    if (!(name %in% names(control)))
      control[[name]] <- default.control[[name]]
  }
  list.diff <- setdiff(names(control),names(default.control))
  if (length(list.diff)>0) {
    warning.text <- "The following arguments in the control list were ignored:"
    for (elt in list.diff)
      warning.text <- paste(warning.text," \"",elt,"\"",sep="")
    warning(warning.text)
  }
  result<-.External("fit_tree_wave",data=use.data,dims=data.dim,use.beta=as.logical(use.beta),weights=as.double(weights),control=control, PACKAGE="treethresh")
  names(result) <- c("splits","details","w","t","membership","beta")
  dimnames(result$splits)=list(NULL,c("id","parent.id","dim","pos","left.child.id","right.child.id","crit","loglikelihood","alpha","C"))
  if (!all(is.na(result$splits[,"alpha"]))) {
    result$splits[,"C"] <- result$splits[,"alpha"]/max(result$splits[,"alpha"],na.rm=TRUE)
  } else {
    result$splits[,"C"] <- NA
  }
  result$t <- matrix(tfromw(result$w),nrow=nrow(result$w))
  result$data <- data
  result$weights <- weights
  result$control <- control
  if (is.null(result$beta)) {
    result$beta <- beta
  } else {
    for (i in 1:length(result$beta))
      result$beta[[i]] <- array(result$beta[[i]],dim=dim(result$data[[i]]))
  }
  for (i in 1:length(result$membership))
    result$membership[[i]] <- array(result$membership[[i]],dim=dim(result$beta[[i]]))
  class(result) <- "wtthresh"
  result
}

subtree.wtthresh <- function(object,C=NULL) {
  if (missing(C) || is.null(C))
    return(object)
  if (all(is.na(object$splits[,"C"])))
    return(object)
  if (is.null(object$splits))
    stop("Not the right type of object")
  nrow <- nrow(object$splits)
  old.types <- as.integer(object$splits[,"id"])
  result <- .C("prune_tree",splits=as.double(object$splits),nrow=as.integer(nrow),kill=integer(nrow),leaf=as.integer(is.na(object$splits[,"dim"])),membership=as.integer(old.types),C=as.double(C),NAOK=TRUE) # ,package="treethresh")
  result$leaf <- as.logical(result$leaf)
  result$kill <- as.logical(result$kill)
  object$splits[result$leaf,3:7] <- NA
  object$details[result$leaf] <- NA
  object$splits <- object$splits[!result$kill,,drop=FALSE]
  object$details <- object$details[!result$kill,,drop=FALSE]
  object$w <- object$w[!result$kill,,drop=FALSE]
  object$t <- object$t[!result$kill,,drop=FALSE]
  for (i in 1:length(object$membership)) {
    len <- length(object$membership[[i]])
    object$membership[[i]] <- array(.C("update_membership",old.membership=as.integer(object$membership[[i]]),new.membership=integer(len),n=as.integer(len),old.types=as.integer(old.types),new.types=as.integer(result$membership),n.types=as.integer(length(old.types)),NAOK=TRUE)$new.membership,dim=dim(object$membership[[i]]))
  }
  object
}


get.w.wtthresh <- function(object,C) {
  if (!missing(C))
    if (!is.null(C))
      object <- subtree.wtthresh(object,C)
  map <- numeric(max(object$splits[,"id"]))
  map[object$splits[,"id"]] <- 1:nrow(object$splits)
  result <- list()
  for (k in 1:length(object$beta)) {
    result[[k]] <- array(object$w[map[as.vector(object$membership[[k]])],k],dim=dim(object$membership[[k]]))
  }
  result
}

get.t.wtthresh <- function(object,C) {
  if (!missing(C))
    if (!is.null(C))
      object <- subtree.wtthresh(object,C)
  map <- numeric(max(object$splits[,"id"]))
  map[object$splits[,"id"]] <- 1:nrow(object$splits)
  result <- list()
  for (k in 1:length(object$beta)) {
    result[[k]] <- array(object$t[map[as.vector(object$membership[[k]])],k],dim=dim(object$membership[[k]]))
  }
  result
}

# Compute best C in cross-validation
prune.wtthresh <- function(object,v=5,sd.mult=0.5,plot=TRUE) {
  reep <- function(max,cur,dim) {
    reep <- 1:max
    if (dim>1) {
      for (pos in 2:dim) 
        if (pos<=cur)
          reep <- rep(reep,each=max)
        else
          reep <- rep(reep,max)
    }
    reep
  }
  create.coords <- function(index,size,dim,factor) {
    coords <- matrix(nrow=length(index)*factor^dim,ncol=dim)
    index <- index-1
    multiple <- 1
    pos <- rep(1,length(index)*factor^dim)
    for (i in 1:dim) {
      cur.coords <- index%%size
      cur.coords <- factor*rep(cur.coords,each=factor^dim)+reep(factor,i,dim)
      coords[,dim+1-i] <- cur.coords
      index <- index%/%size
      pos <- pos+(cur.coords-1)*multiple
      multiple <- multiple*size*factor
    }
    list(coords=coords,pos=pos)
  }
  get.sizes <- function(object) {
    C <- object$splits[-1,"C"]
    C.values <- sort(unique(C))
    c(apply(outer(C,C.values,function(a,b) as.numeric(a>=b)),2,sum)/2+1,1)
  }
  if (is.na(v))
    return(object)
  if (all(is.na(object$splits[,"C"])))
    return(object)
  sizes <- sapply(object$beta,length)
  data.dim <- length(dim(object$beta[[1]]))
  # Split the dataset into v portions
  permut <- sample(min(sizes))
  sample.size <- min(sizes)/v
  C.values <- sort(as.numeric(na.omit(c(unique(object$splits[,"C"]),2))))
  results <- matrix(0,nrow=v,ncol=length(C.values))
  to <- 0
  for (j in 1:v) {
    from <- to + 1
    to <- round(sample.size*j)
    test.index <- permut[from:to]
    a.beta <- object$beta
    for (k in 1:length(a.beta)) {
      a.beta[[k]][create.coords(test.index,min(sizes)^(1/data.dim),data.dim,(sizes[k]/min(sizes))^(1/data.dim))$pos] <- NA
    }
    a.tree <- wtthresh(beta=a.beta,control=object$control,weights=object$weights)
    for (l in 1:length(C.values)) {
      a.tree <- subtree.wtthresh(a.tree,C=C.values[l])
      w <- get.w.wtthresh(a.tree)
      for (k in 1:length(a.beta)) {
        test.pos <- create.coords(test.index,min(sizes)^(1/data.dim),data.dim,(sizes[k]/min(sizes))^(1/data.dim))$pos
        beta <- pmin(1e5,object$beta[[k]][test.pos])
        results[j,l] <- results[j,l] + sum(log(1+w[[k]][test.pos]*beta))
      }
    }
  }
  results <- results - results[,ncol(results)] 
  results.mean <- apply(results,2,mean)
  results.sd <- apply(results,2,sd)
  lower.bound <- (results.mean-sd.mult*results.sd)[which.max(results.mean)]
  if (plot) {
    results.range <- c(min(results.mean-results.sd),max(results.mean+results.sd))
    ns <- 1:length(results.mean)
    plot(ns,results.mean,xaxt="n",ylab="Loglikelihood from CV (+const)",type="o",xlab="C",ylim=results.range)
    segments(ns, results.mean - results.sd, ns, results.mean + results.sd)
    axis(1, at = ns, labels = c(format(C.values[-length(C.values)],digits=3),">1"))
    mtext("number of regions", side = 3, line = 3)
    axis(3, at = ns, labels = as.character(get.sizes(object)))
    abline(h=lower.bound,lty=3)
  }
  subtree.wtthresh(object,C.values[max(which(results.mean>=lower.bound))])
}

thresh.wtthresh <- function(object,data,C,postmed=TRUE) {
  if (missing(data) && !is.null(object$data))
    data <- object$data 
  if (postmed) {
    w <- get.w.wtthresh(object,C)
    if (length(data)!=length(w))
      stop("Data does not have the right length.")
    for (i in 1:length(data)) 
      data[[i]] <- postmed.laplace(data[[i]],w[[i]])
  } else {
    t <- get.t.wtthresh(object,C)
    if (length(data)!=length(t))
      stop("Data does not have the right length.")
    for (i in 1:length(data)) 
      data[[i]][abs(data[[i]])<t[[i]]] <- 0
  }
  data
}


########################### UTILITY FUNCTIONS

extract.coefficients.imwd <- function(object,start.level=5) {
  if (start.level>=object$nlevels)
    return(list())
  types <- c("CD","DD","DC")
  result <- list()
  for (level in start.level:object$nlevels) 
    for (type in types) {
      nm <- lt.to.name(level=level-1,type)
      result[[nm]] <- array(object[[nm]],dim=rep(2^(level-1),2))
    }
  result
}

estimate.sdev.imwd <- function(object,dev=mad) {
  types <- c("CD","DD","DC")
  mads <- c()
  for (type in types) {
    nm <- lt.to.name(level=object$nlevels-1,type)
    elt <- object[[nm]]
    our.mads <- dev(elt)
    mads <- c(mads,as.vector(our.mads))
  }
  median(mads)
}
    
insert.coefficients.imwd <- function(object,update) {
  for (name in names(update))
    object[[name]] <- update[[name]]
  object
}

extract.coefficients.wd <- function(object,start.level=5) {
  if (start.level>=object$nlevels)
    return(list())
  result <- list()
  for (level in start.level:object$nlevels) {
    foo <- accessD(object,level=level-1)
    nm <- as.character(level)
    result[[nm]] <- array(foo,dim=length(foo))
  }
  result
}

estimate.sdev.wd <- function(object,dev=mad) {
  elt <- accessD(object,level=object$nlevels-1)
  dev(elt)
}
    
insert.coefficients.wd <- function(object,update) {
  for (name in names(update))
    object <- putD(object,level=as.numeric(name)-1,update[[name]])
  object
}

extract.coefficients <- function(object,start.level=5)
  UseMethod("extract.coefficients")

estimate.sdev <- function(object,dev=mad)
UseMethod("estimate.sdev")

estimate.sdev.numeric <- function(object, block=8, dev=mad) {
  if (is.na(block))
    return(dev(object))
  block <- min(c(block,length(object)))
  n.blocks <- floor(length(object)/block)
  our.mads <- numeric(n.blocks)
  for (i in 1:n.blocks) {
    dta <- object[block*(i-1)+1:block]
    our.mads[i] <- dev(dta)
  }
  median(our.mads)
}

insert.coefficients <- function(object, update)
  UseMethod("insert.coefficients")

wavelet.treethresh <- function(object,sdev=NA,dev=mad,start.level=5,levelwise=FALSE,v=5,sd.mult=0.5,postmed=TRUE,...) {
  # If no sdev is provided estimate it
  if (is.na(sdev)) {
    sdev <- estimate.sdev(object,dev=dev)   
  }
  # Extract the coefficients
  coefs <- extract.coefficients(object,start.level=start.level)
  # Rescale the coefficients
  for (nm in names(coefs))
    coefs[[nm]] <- coefs[[nm]] / sdev
  if (!levelwise) {
    # JOINT THRESHOLDING
    # Compute the thresholding tree
    coefs.tree <- wtthresh(coefs,...)
    # Prune the tree
    coefs.pruned.tree <- prune.wtthresh(coefs.tree,v=v,sd.mult=sd.mult,plot=FALSE)
    # Threshold according to the pruned tree
    coefs.threshed <- thresh(coefs.pruned.tree,postmed=postmed)
  } else {
    # LEVELWISE THRESHOLDING
    coefs.threshed <- list()
    for (nm in names(coefs)) {
      coefs.tree <- treethresh(coefs[[nm]],...)
      coefs.pruned.tree <- prune(coefs.tree,v=v,sd.mult=sd.mult,plot=FALSE)
      coefs.threshed[[nm]] <- thresh(coefs.pruned.tree,postmed=postmed)
    }
  }     
  # Undo the resclaing
  for (nm in names(coefs))
    coefs.threshed[[nm]] <- coefs.threshed[[nm]] * sdev
  # Update coefficients
  insert.coefficients(object,coefs.threshed)
}

