# Copyright 2011-, Tobias Abenius, Rebecka Jornsten, Sven Nelander
#
# S4 classes not used because of ignorance
#setClass("EPOCA",contains="list")
#setClass("EPOCG",contains="list")
nodiag <- function(x) { 
  if (is(x[1,1], 'logical')) {
    diag(x) <- FALSE
  } else {
    diag(x) <- 0
  }
  x
}
c.lambda = "\u03BB" #small lambda
c.square = "\u00B2" #superscript square
c.infty = "\u221E"  #infinity
progressbar.width = 50
progressbar <- function(i,k,p,q,progress) {
  progress.old <- round(progress,digits=0)
  progress <- round(progressbar.width * (p * (k-1) + i) / (p*q),digits=0)
  if (progress.old == 0 | progress > progress.old) {
    cat(paste('\r|',paste(rep('=', progress),collapse=''), '>',paste(rep(' ', progressbar.width-progress),collapse=''),'|',sep=''))
  }
  return(progress)
}
plapply <- function(X1,X2,FUN, ...) {
  FUN <- match.fun(FUN)
  if (length(X1) != length(X2)) stop("x1 and x2 are not of same length")
  if (!is.vector(X1) || is.object(X1)) X1 <- as.list(X1)
  if (!is.vector(X2) || is.object(X2)) X2 <- as.list(X2)
  l.new <- list()
  for (k in 1:length(X1))
    l.new[[k]] <- FUN(X1[[k]], X2[[k]], ...)
  return(l.new)
}
reg <- function(y,x) {
#  mode <- 3
#  if (mode==1) {
#    require('corpcor')
#    pinvu <- pseudoinverse(as.matrix(t(x)))
#    d <- as.matrix(t(y)) %*% pinvu
#  } else if (mode==2) {
#    d <- coef(lm(y~x))[1]
#  } else if (mode==3) {
    d <- coef(lsfit(x,y,intercept=T))[2]
#  }
  d
}
coef.EPOCA <- function(object, k=1, ...) {
  object$coefficients[[k]]
}
predict.EPOCG <- function(object,newdata,k=1,trace=0, ...) {
  if (typeof(newdata) == "list") # works for data.frame also
    U <- newdata$U
  else
    U <- newdata
  if (is.null(U)) stop("predict require U")
  m <- dim(U)[1]
  p <- dim(U)[2]
  if (trace > 0) cat("Predicting for p =",p,"variables, in",m,"points\n")
  muYresM <- t(array(rep(object$Yres.mean, m),dim=c(p,m)))
#  YonU <- array(0,dim=c(N,p))
#  for (i in 1:p) YonU[,i] <- d[i]*U[,i]
  (U - object$U.mean) %*% coef(object,k=k, ...) + muYresM  #+ YonU
}
predict.EPOCA <- function(object,newdata,k=1,trace=0, ...) {
  if (typeof(newdata) == "list") { # works for data.frame also
    Y <- newdata$Y
    if (is.null(Y))
      stop("Y is required")
    N <- dim(Y)[1]
    p <- dim(Y)[2]
    U <- newdata$U
    if (is.null(U))
      U <- array(object$U.mean,dim=c(N,p))
    else
      U <- U - object$U.mean
  } else {
    Y <- newdata
    N <- dim(Y)[1]
    p <- dim(Y)[2]
    U <- array(object$U.mean,dim=c(N,p))
  }
  Y <- Y - object$Y.mean
  if (is.null(Y)) stop("predict require Y")
  if (trace > 0)
    cat("\nUsing direct effects from object\n")
  YonU <- U%*%diag(object$d) 
 
  if (trace > 0) cat("Predicting for p =",p,"variables, in",N,"points\n")
  Y %*% coef(object,k=k, ...) + YonU
}
summary.EPOCA <- function(object, k=NULL,...) {
  if (is.null(k)) {
    K <- length(object$lambdas)
    q <- length(object$coefficients)
    ks <- 1:q
  } else {
    K = 1
    ks <- k
  }
  sp <- array(0,dim=c(K,5))
  dimnames(sp) <- list(paste("lambda=",round(object$lambdas[ks],digits=4),sep=''),c(paste("R2",sep=''),"Cp",'BIC',"RSS","links"))
  sp[,1] <- round(object$R2[ks],digits=4)
  sp[,2] <- round(object$Cp[ks],digits=4)
  sp[,3] <- round(object$BIC[ks],digits=4)
  sp[,4] <- round(object$RSS[ks],digits=4)
  p <- dim(coef(object))[1]
  links <- object$links
  #links <- function (B) as.integer(sum( (B * (1 - diag(rep(1,p)))) != 0))
  #links <- function (B) as.integer(sum(B != 0))
  q <- length(object$coefficients)
  i = 1
  for(k in ks)  {
    sp[i,5] <- links[k] #links(object$coefficients[[k]])
    i = i + 1
  }
  ans <- list(call=object$call,models=sp,SS.tot=object$SS.tot,d=object$d)
  class(ans) <- "summary.EPOCA"
  return(ans)
}
print.summary.EPOCA <- function(x, ...) {
  digits = max(3, getOption("digits") - 3)
  cat("\nCall:\n",deparse(x$call), "\n", sep="")
#  cat("\nDirect effects: \n")
#  print.default(x$d, ...)
  cat("\nModels:\n")
  print.default(x$models, print.gap=2, quote=F, ...)
  cat("\nSStot:",x$SS.tot,"\n")
  cat("\n")
}
print.EPOCA <- function(x, ...) {
  #require(methods)
  digits = max(3, getOption("digits") - 3)
  cat("\nCall:\n",deparse(x$call), "\n\n", sep="")
  cat("Coefficients:\n")
  K <- length(x$lambdas)
  for(k in 1:K) {
    print(paste("For ", c.lambda, "=", x$lambdas[k], sep=''))
    #print(format(coef(x,k=k), digits=digits), print.gap=2, quote=FALSE, ...)
    print(coef(x,k=k))#, note.dropping.colnames=T,...)
  }
  cat("\nDirect effects: \n")
  print.default(x$d, ...)
  s <- summary(x)
  cat("\nModels: \n")
  print.default(s$models, print.gap=2, ...)
  cat("\nSStot:",x$SS.tot,"\n")
  cat("\n")
  invisible(x)
}
print.EPOCG <- function(x, ...) {
  #require(methods)
  digits = max(3, getOption("digits") - 3)
  cat("\nCall:\n",deparse(x$call), "\n", sep="")
  cat("\nCoefficients:\n")
  K <- length(x$lambdas)
  for(k in 1:K) {
    print(paste("For ", c.lambda, "=", x$lambdas[k], sep=''))
    print(coef(x,k=k)) #, note.dropping.colnames=T, ...)
    #print.default(format(coef(x,k=k), digits=digits), print.gap=2, quote=FALSE, ...)
  }
  s <- summary(x)
  cat("\nDirect effects: \n")
  print.default(x$d, ...)
  cat("\nModels: \n")
  print.default(s$models, print.gap=2, ...)
  cat("\nSStot:",x$SS.tot,"\n")
  cat("\n")
  invisible(x)
}
as.graph.EPOCA <- function(model, k=1) {
  #require('graph')
  p <- dim(coef(model,k=k))[1]
  A <- nodiag(abs(coef(model,k=k))) # * (1 - diag(array(1,dim=p)))
  return( new("graphAM", adjMat=A, edgemode='directed') )
}
write.sif <- function(model, k=1, file="", append=F) {
  if (file == "") 
    file <- stdout()
  else if (is.character(file)) {
    file <- file(file, ifelse(append, "a", "w"))
    on.exit(close(file))
  }
  else if (!isOpen(file, "w")) {
    open(file, "w")
    on.exit(close(file))
  }
  if (!inherits(file, "connection")) 
    stop("'file' must be a character string or connection")

  if (class(model) == 'matrix') {
    p <- dim(model)[1]
    A <- nodiag(model)
  } else {
    p <- dim(model$coefficients[[1]])[1]
    A <- nodiag(coef(model, k=k))
  }
  newnames <- gsub(' ','_',rownames(A))
  for (i in 1:p)
    for(j in 1:p) {
      e <- zapsmall(A[i,j],digits=3)
      if (e != 0)
	if (e < 0) {
	  cat(newnames[i],'inhibits',newnames[j],"\n",file=file)
	} else {
	  cat(newnames[i],'stimulates',newnames[j],"\n",file=file)
	}
    }
}
plot.EPOCA <- function (x, layout=NULL, k = 1, showtitle=F, bthr=0, showself=F, type=c('graph','modelsel'),...) {
  typeOfPlot <- match.arg(type)
  if (typeOfPlot == 'modelsel') {
    cl <- match.call()
    modelselPlot(x, ...)
  }

  if (!is.null(attr(x,'epocboot')) | !is.null(attr(x,'epocboot1')) | class(x) == 'matrix' | inherits(x,'Matrix')) {
    if (!is.null(attr(x,'epocboot'))) {
      if (k > length(x)) stop(paste("This EPoC bootstrap object contains only",length(x),"sparsity levels"))
      x <- x[[k]]
      #print(class(x))
    }
    #class(x) <- attr(x,'origclass')
    #print(class(x))
    if (showself)
      adjm <- x
    else
      adjm <- nodiag(x)
    if (bthr!=0) {
      adjm[abs(adjm) <= bthr] = 0
    }
  } else {
    if (showself)
      adjm <- coef(x, k=k) 
    else
      adjm <- nodiag(coef(x, k=k))
  }
  p1 <- dim(adjm)[1]
  p2 <- dim(adjm)[2]
  p <- max(p1,p2)
  adjmBack <- adjm
  adjm <- as.matrix(adjm)
  adjm <- t(array(t(adjm),dim=c(p,p)))
  adjm[p1:p2,] <- 0
  if (F) { # Disabled: rownames will be wrong
    # keep only rows and columns with nonzeros
    ii1 <- apply(adjm!=0,1,sum)
    ii1 <- (1:p1)[ii1 > 0]
    ii2 <- apply(adjm!=0,2,sum)
    ii2 <- (1:p2)[ii2 > 0]
    if (length(ii1) > 0 | length(ii2) > 0) {
      adjm <- adjm[ii1,,drop=F]
      adjm <- adjm[,ii2,drop=F]
    }
  } else  {
    ii <- (1:p)[apply(adjm!=0,1,sum) + apply(adjm!=0,2,sum) > 0]
    adjm <- adjm[ii,,drop=F]
    adjm <- adjm[,ii,drop=F]
    vx <- colnames(adjmBack)[ii]
  }
  # columns are targets in graph adjacency matrices

  pp <- length(ii)
  g <- new("graphNEL", nodes=vx,edgemode="directed")
  edgecolor <- NULL
  attrib1 <- NULL
  for (i in 1:pp) {
    for (j in 1:pp) {
      if (adjm[i,j] > 0) {
	g <- addEdge(vx[i], vx[j], g, 1)
	edgecolor <- c(edgecolor,color='green')
	attrib1 <- c(attrib1,foo='normal')
      } else if (adjm[i,j] < 0) {
	g <- addEdge(vx[i], vx[j], g, 1)
	edgecolor <- c(edgecolor,color='red')
	attrib1 <- c(attrib1,foo='tee')
      }
    }
  }
  names(attrib1) <- edgeNames(g)
  names(edgecolor) <- edgeNames(g)
  edgeAttrs <- list(arrowhead=attrib1,color=edgecolor)

  N <- nodes(g); 
  shapes <- rep('box',length(N)); 
  names(shapes) <- N
  ww <- rep(1.5,length(N))
  names(ww) <-  N
  nodeAttrs <- list(width=ww,shape=shapes)

  plot(g,'fdp',edgeAttrs=edgeAttrs,nodeAttrs=nodeAttrs)
  if (showtitle) title(x$call)
}
epoc.bootplot <- plot.EPOCA
epoc.svdplot <- function(G.svd,C=1) {
  x <- G.svd
  ii2 <- x$ii[x$spload.in[,C]!=0 | x$spload.out[,C]!=0]
  newg <- t(x$m)[ii2,ii2]
  plot.EPOCA(newg,showself=F)
}
plot.EPOCBOOT1 <- plot.EPOCA
plot.EPOCG <- plot.EPOCA
coef.EPOCG <- coef.EPOCA
print.summary.EPOCG <- print.summary.EPOCA
summary.EPOCG <- summary.EPOCA
as.graph.EPOCG <- as.graph.EPOCA

epoc.lambdamax <- function(X,Y,getall=F,predictorix=NULL) {
  dims <- dim(Y)
  if (!is.null(predictorix)) {
    P <- length(predictorix)
    rix <- array(0,dim=dims[2])
    rix[predictorix] <- 1:P
  }
  if (is.null(dims)) {
    return( norm(crossprod(X,Y),'i') )
  } else {
    n <- dims[2] #number of variables = genes
    lambdamax <- array(NaN,dim=n)
    for (k in 1:n) {
      predk = X
      if(is.null(predictorix))
	predk[,k] <- 0
      else if (k %in% predictorix)
	predk[,rix[k]] <- 0
      lambdamax[k] <- norm(crossprod(predk, Y[,k]),'i')
    }
    if (getall) return (lambdamax)
    else return(max(lambdamax))
  }
}
epoc.bootstrap <- function(Y,U,nboots=100,bthr=NULL,method='epocG',...) {
  first = T
  N <- dim(Y)[1]
  for (i in 1:nboots) {
    ix <- sample(1:N, N, replace=T)
    if (method == 'epocG')
      mod.boot1 <- epocG(Y[ix,],U[ix,],...)
    else
      mod.boot1 <- epocA(Y[ix,],U[ix,],...)
    D2 <- lapply(mod.boot1$coefficients, function (A) (A != 0)*1)
    if (first) {
      D <- D2
      q <- length(mod.boot1$lambdas)
      lambdas <- mod.boot1$lambdas
      first = F
    } else {
      D <- plapply(D, D2, function(A,B) A+B)
    }
  }
  D <- lapply(D,function (A) A/nboots)
  if (!is.null(bthr))
    D <- lapply(D,function(A) A >= bthr)
  for (i in 1:q) {
#    attr(D[[i]],'epocboot1') <- T
#    attr(D[[i]],'origclass') <- class(D[[i]])
#    class(D[[i]]) <- 'EPOCBOOT1'
  }
#  attr(D,'epocboot') <- T
  attr(D, 'nboots') <- nboots
  attr(D, 'lambdas') <- lambdas
  class(D) <- 'bootsize'
  D
}
epoc.final <- function(epocboot, bthr=0.2, k) {
  if (length(epocboot) < k) 
    stop(paste("epoc.final: the requested k exceeds the number of graphs created by epoc.bootstrap. which are",length(epocboot)))
  return(epocboot[[k]] >= bthr)
}
epoc <- function(method=c('G','A'),Y,U,lambdas=NULL,inorms=NULL,thr=1e-10,trace=0, ...) {
  #require('lassoshooting')
  #require('Matrix')
  cl <- match.call()
  if (!is.null(cl[['lambda']]))
    stop("The parameter `lambda' should be called `lambdas'.")
  if (!is.null(cl[['predictorix']]))
    stop("The parameter `predictorix' is now obsolete, see manual.")
  method <- match.arg(method)

  N <- dim(Y)[1] #number of equations = experiments
  #number of variables = genes
  p <- dim(Y)[2]
  if (trace > 0) cat("Solving for p =",p,"variables,",N,"equations\n")
  mlist <- strsplit(method,'.',fixed=T)
  method = mlist[[1]][1]
  method2 = ifelse(length(mlist[[1]]) > 1, mlist[[1]][2], '')
  #if (method!='G' & method!='A') stop(paste("Unknown method ",method,"!",sep=''))

  if (is.null(lambdas)) {
    lambdas = switch(method,
		     G=c(0.99999, 1.25^(-(1:10))),
		     A=c(exp(-(0:12)/6),exp(-(5:9)/2)))
  }
  lambdas <- rev(sort(lambdas))
  q <- length(lambdas)

  hasU <- !is.null(U)
  if (hasU) {
#  if (is.null(predictorix)) {
    P <- dim(U)[2]
    predictorix <- 1:P
#    P <- p
#  } else {
#    P <- length(predictorix)
  } else {
    P <- p
  }
  reindex <- array(0,dim=p)
  reindex[predictorix] <- 1:P

  if (trace > 0) cat("Centering...")
  muY <- colMeans(Y)
  Y <- Y - muY
  if (hasU) {
    muU <- colMeans(U)
    U <- U - muU
  } else {
    U <- array(0,dim=c(N,p))
    muU <- array(0,dim=c(p))
  }
  if (trace > 0) cat("DONE\n")
  #/center
  #regressing Y on U
  if (trace > 0) cat("Regressing Y on U...")
  if (hasU) {
    regf <- function(i) { # 20111004: Fixed bug with predictorix here
	return(reg(Y[,predictorix[i]], U[,i]))
    }
    d <- sapply(1:P,regf)
    d <- pmax(d,0)
  } else {
    d <- array(0,dim=p)
  }
  if (is.null(dimnames(Y)) || is.null(dimnames(Y)[[2]])) {
    gs <- paste('V',1:p,sep='')
  } else {
    gs <- dimnames(Y)[[2]]
  }
  names(d) <- gs[predictorix]
  if (trace > 0) cat("DONE\n")
  if (trace > 3) cat("Direct effects of CNA:",d,"\n")

  if (trace > 0) cat("Correcting for direct effects...")
  YonU <- array(0,dim=c(N,p))
  for (i in 1:P) YonU[,predictorix[i]] <- d[i]*U[,i]
#  YonU <- U%*%diag(d)  # requires much more memory
  Yres <- Y - YonU
  muYres <- colMeans(Yres)
  Yres <- Yres - muYres
  if (trace > 0) cat("DONE\n")
  if(method=='G') {
    pred <- U
    resp <- Yres
  } else {
    pred <- Y
    resp <- Yres
  }

  #finding maximum lambda
  if (is.null(inorms)) {
    if (trace > 0) cat("Finding ",c.lambda,"_max...",sep='')
    inorms <- epoc.lambdamax(pred,resp,getall=T, predictorix=predictorix)
  } else {
    if (trace > 0) cat("Using provided ||X'y||",c.infty,sep='')
    if (length(inorms) != p) stop("Parameter inorms should be NULL or a p-vector")
  }

  if (trace > 3) cat("\n||.||",c.infty,": ",inorms,"\n",sep='')
  lambdamax <- max(inorms)
  extreme.gene <- colnames(Y)[which.max(inorms)]
  if (trace > 0) {
    cat(paste("DONE\nRel.",c.lambda,"s:",sep=''),paste(lambdas,sep=', '),"\n")
    cat("extreme gene: ",extreme.gene,", lambdamax = ",lambdamax,"\n",sep='')
  }
  #/finding maximum lambda

  #B <- array(NaN,dim=c(n,n,q))
  B <- list()
  s2 <- array(NaN,dim=q)
  RSS <- array(NaN,dim=q)
  R2 <- array(NaN,dim=q)
  Cp <- array(NaN,dim=q)
  RMSD <- array(NaN,dim=q)
  BIC <- array(NaN,dim=q)
  E <- list()

  if (trace > 0) cat("Gram matrix calculation of predictors...")
  #XtX <- t(pred) %*% pred
  XtX <- crossprod(pred,pred) #20111004: Should be faster
  if (trace > 0) cat("DONE\n")
  if (trace == 1) cat("Lasso regression...")

  lasso <- lassoshooting
  progress <- 0
  pix <- predictorix # orginal predictorix (kludge for epocA)
  if (method=='A') {
    P <- p #!!!
    predictorix=1:p
  }
  links <- function (B) as.integer(sum(B != 0))
  for(k in 1:q) {
    if (method=='G') {
      B1 <- Matrix(0,nrow=P, ncol=p)
      diag(B1)[1:P] <- d
      dimnames(B1) <- list(gs[predictorix], gs)
    }
    else {
      B1 <- Matrix(0,nrow=p,ncol=p,sparse=T)
      dimnames(B1) <- list(gs, gs)
    }

    lambda <- lambdas[k] * lambdamax
    for(i in 1:p) {
      if (trace == 2) progress <- progressbar(i,k,p,q,progress)
      if (inorms[i] >= lambda) {
	pred.noi <- pred 
	if (i %in% predictorix)
	  pred.noi[,match(i,predictorix)] <- 0 # if we don't update this, lambdamax is wrong, perhaps lassoshooting forcezero should set this column to 0
	#XtY <- t(pred.noi) %*% resp[,i] 
	XtY <- crossprod(pred.noi, resp[,i]) #20111004: this should be faster than %*%
	if (trace == 3) cat("for var i =",i,"lasso...\n")
	l <- lasso(xtx=XtX,xty=XtY,lambda=lambda,forcezero=i,thr=thr)
	b <- l$coefficients
	if (trace == 3) cat("for var i =",i," lasso done\n")
	nonz <- (1:P)[abs(b) >= thr] # 2011-06-13
	#nonz <- setdiff(nonz,i) # FIXME
	betas <- Matrix(0,nrow=P,ncol=1,sparse=T) # sparse M don't go well with arrays..
	dimnames(betas)[[1]] <- gs[predictorix]
	if (length(nonz)>0){
	  betas[nonz,1] <- b[nonz]
	  B1[,i] <- B1[,i] + betas
	}
	if(i %in% pix)
	  B1[reindex[i],i] <- d[reindex[i]] # without this bug occurs (diagonal disappear when updating betas)
      }
    }
    if (trace==3) cat("for lambda_k, k =",k,"\n")
    B[[k]] <- B1
    if (k==1) 
      SS.tot <- sum((Y - colMeans(Y))^2)
    muYresM <- t(array(rep(muYres,N),dim=c(p,N)))
    if (method=='G') {
      yhat <- U %*% B1 + muYresM  #+ YonU # this is on the diagonal in G
    } else {
      yhat <- Y %*% B1 + muYresM + YonU 
    }
    e <- Y - yhat
    E[[k]] <- e
    if(trace >= 3) {
      coryy <- array(0,dim=p)
      for(i in 1:p) {
	coryy[i] <- cor(Y[,i],yhat[,i])
	if (trace > 3) cat("cor(y,y^)_i=",i,":",coryy[i],"\n")
      }
    }
    if (trace >= 3) {
      cat(k,"avgcor:",mean(coryy,na.rm=T))
      cat(", mincor:",min(coryy,na.rm=T),"\n")
    }
    RSS[k] <- sum(e^2)
    R2[k] <- 1 - RSS[k] / SS.tot
    s2[k] <- RSS[k] / (N - P)
    p.subset <- links(B1)
    RMSD[k] <- sqrt(RSS[k]/N)
  }
  for (k in 1:q) {
    d.subset <- sum(d != 0)
    dense <- which.min(RSS)
    P <- links(B[[dense]])
    p.subset <- links(B[[k]])
    Cp[k] <- RSS[k]/(min(RSS) / (N*p-P)) + (2*(p.subset+d.subset))
    BIC[k] <- N*p*log(RSS[k] / (N*p-p.subset-d.subset)) + (p.subset+d.subset)*log(N*p)
  }
  if (trace > 0) cat("\rDONE",rep(' ',progressbar.width),'\n',sep='')
  links <- function (B) as.integer(sum(B!=0) - sum(diag(B[,pix])!=0))
  links <- unlist(lapply(B, links))

  obj <- list(call=cl, coefficients=B, lambdas=lambdas, lambdamax=lambdamax, d=d, Y.mean=muY, U.mean=muU, Yres.mean=muYres, R2=R2, Cp=Cp, SS.tot=SS.tot, RSS=RSS, RMSD=RMSD, s2=s2, links=links, predictorix=pix, inorms=inorms, BIC=BIC, E=E)
  class(obj) <- switch(method, G="EPOCG", A="EPOCA")
  obj
}
epocA <- function(Y,U=NULL,lambdas=NULL,thr=1e-10,trace=0, ...) {
  o <- epoc('A',Y,U,lambdas,thr=thr,trace=trace,...)
  ### hairy code to remove NULLs from parameter list
  cl <- match.call()
  nnulls <- rep(T,length(names(cl)))
  for (i in 2:length(names(cl))) {
    key <- names(cl)[i]
    v <- eval(parse(text=key))
    nnulls[i] <- !is.null(v)
  }
  ii <- (1:length(names(cl)))[nnulls]
  cl2 <- cl[ii]
  #####
  o$call <- cl2
  o
}
epocG <- function(Y,U,lambdas=NULL,thr=1e-10,trace=0, ...) {
  o <- epoc('G',Y,U,lambdas,thr=thr,trace=trace,...)
  ### hairy code to remove NULLs from parameter list
  cl <- match.call()
  nnulls <- rep(T,length(names(cl)))
  for (i in 2:length(names(cl))) {
    key <- names(cl)[i]
    v <- eval(parse(text=key))
    nnulls[i] <- !is.null(v)
  }
  ii <- (1:length(names(cl)))[nnulls]
  cl2 <- cl[ii]
  #####
  o$call <- cl2
  o
}
crossvalix <- function(N,K) {
  if (K > N) stop("crossvalix: K > N !")
  folds <- list()
  left <- 1:N
  for (k in 1:(K-1)) {
    folds[[k]] <- sample(left, round(N/K,digits=0), replace=FALSE)
    left <- setdiff(left, folds[[k]])
  }
  folds[[K]] <- left
  folds
}
epoc.validation <- function(type=c('pred','concordance'),repl,Y,U,lambdas=NULL,method='G',thr=1e-10,trace=0,...) {
  N <- dim(Y)[1]
  first <- T
  type = match.arg(type)
#  trace = as.integer(match.arg(trace))
  progress <- 0
  if (type == 'pred') {
    if (trace > 0) cat ("epoc.validation: type =",type)
    K = 10
    dnames = NULL
    for(b in 1:repl) {
      folds <- crossvalix(N,K)
      for (k in 1:K) {
	if (trace > 0) progress <- progressbar(k,b,K,repl,progress)
	Y.tr <- Y[-folds[[k]],]
	U.tr <- U[-folds[[k]],]
	Y.te <- Y[folds[[k]],]
	U.te <- U[folds[[k]],]
	G.tr <- epoc(method, Y.tr, U.tr, lambdas=lambdas,thr=thr,...)
	if (first) {
	  q = length(G.tr$links)
	  E <- array(0,dim=c(K,repl,q))
	  links <- array(0,dim=c(K,repl,q))
	  first <- F
	}
	for (i in 1:q) {
	  E[k,b,i] <- sum(sum((Y.te - predict(G.tr,list(U=U.te),k=i))^2))
	  links[k,b,i] <- G.tr$links[i]
	}
      }
      dnames <- c(dnames,paste("repl",b,"fold",1:K))
    }
#    dimnames(E)[[1]] <- dnames
#    dimnames(E)[[2]] <- paste("splevel",1:q)
#    dimnames(links)[[1]] <- dnames
#    dimnames(links)[[2]] <- paste("splevel",1:q)
    lambdas <- G.tr$lambdas
    B <- dim(E)[2]
    l <- dim(E)[3]
#    flat <- array(0,dim=c(B*l,3))
#    dimnames(flat)[[2]] <- c('e','lambda','links')
#    flat[,1] <- t(E)
#    flat[,2] <- lambdas
#    flat[,3] <- t(links)
#    flat <- data.frame(flat)
    E <- apply(E,c(2,3),mean)
    links <- apply(links,c(2,3),mean)
#    opt <- findoptimvalid.pred(list(E=E, lambdas = lambdas, links=links, K=K, B=B),fromvalid=T)
    o <- list(E=E,lambdas=lambdas,links=links,B=B,K=K)
    opt <- findoptimvalid.pred(o)
    o$sopt <- opt$sopt
    o$lopt <- opt$lopt
    o$fit.loglinkslambda <- opt$mm2
    class(o) <- "EPoC.validation.pred"
    return(o)
  } else if (type == 'concordance') {
    #require(irr)
    if (trace > 0) cat ("epoc.validation: type =",type)
    first = T
    for(b in 1:repl) {
      folds <- crossvalix(N,2)
      Y.tr <- Y[-folds[[1]],]
      U.tr <- U[-folds[[1]],]
      Y.te <- Y[folds[[1]],]
      U.te <- U[folds[[1]],]
      G1 <- epoc(method, Y.tr, U.tr, lambdas=lambdas,thr=thr,trace=trace,...)
      G2 <- epoc(method, Y.te, U.te, lambdas=lambdas,thr=thr,trace=trace,...)
      K <- length(G1$links)
      for (i in 1:length(G1$links)) {
	if (trace > 0) progress <- progressbar(i,b,K,repl,progress)
	G1.v <- as(coef(G1,k=i),'sparseVector')
	G2.v <- as(coef(G2,k=i),'sparseVector')
	nzs <- (G1.v != 0 | G2.v != 0) # according to zadjust == 2 in crossval_scores.m
	r <- Matrix(0,nrow=sum(nzs),ncol=2)
	r[,1] <- G1.v[nzs]
	r[,2] <- G2.v[nzs]
	Kend <- kendall(r, correct=TRUE)
	if (trace> 0) print(Kend)
	if (first) {
	  first = F
	  q <- length(G1$links)
	  W <- array(0,dim=c(repl,q))
	  links <- array(0,dim=c(repl,q,2))
	  lambdas <- G1$lambdas
	}
	links[b,,1] <- G1$links
	links[b,,2] <- G2$links
	W[b,i] <- Kend$value
      } 
    }
    o = list(W=W, lambdas = lambdas, links=links, B=B)
    opt <- findoptimvalid(o)
    o$sopt <- opt$sopt.hi
    o$lopt <- opt$lopt.lo
    o$fit.loglinkslambda <- opt$mm2
#    o$sopt.lo <- opt$sopt.lo
#    o$sopt.hi <- opt$sopt.hi
#    o$lopt.lo <- opt$lopt.lo
#    o$lopt.hi <- opt$lopt.hi
    class(o) <- "EPoC.validation.W"
    return(o)
  } else {
    stop(paste("epoc.validation: unknown type:",type, ", valid types are 'concordance' and 'pred'"))
  }
}
findoptimvalid <- function(object) {
#  object <- object$opt
  deg <- 1
  mm<-loess(as.vector(t(object$W))~rep(as.vector(object$lambdas),object$B),degree=deg)
  o <- order(mm$x)
  lseq<-seq(min(object$lam),max(object$lam),by=.01)
  pp<-predict(mm,newdata=lseq,se=T)
  wmax<-max(pp$fit)
  lupp<-pp$fit+pp$se
  lopt.lo<-min(lseq[lupp>=wmax])
  lopt.hi<-max(lseq[lupp>=wmax])

  links <- (object$links[,,1]+ object$links[,,2])/2
#  cat("findoptimvalid: links: ", links,"\n")
#  cat("findoptimvalid: lambdas: ",object$lambdas,"\n")
#  cat("findoptimvalid: B: ",object$B,"\n")
#  recover()
  mm2<-loess(log(as.vector(t(links)))~rep(as.vector(object$lambdas),object$B),degree=deg)
  sopt.hi<-exp(predict(mm2,newdata=lopt.lo))
  sopt.lo<-exp(predict(mm2,newdata=lopt.hi))
  lopt <- lopt.lo
  sopt <- sopt.hi
  o2<-order(mm2$x)
  return(list(sopt=sopt,sopt.lo=sopt.lo,sopt.hi=sopt.hi,sopt=sopt,lopt.hi=lopt.hi,lopt.lo=lopt.lo,lopt=lopt,lseq=lseq,pp=pp,mm=mm,o=o,o2=o2,mm2=mm2,lupp=lupp,wmax=wmax))
}
findoptimvalid.pred <- function(object) {
#  o <- findoptimvalid.pred(object)
  deg <- 1
  mm<-loess(as.vector(t(object$E))~rep(as.vector(object$lambdas),object$B),degree=deg)
  o <- order(mm$x)
  lseq<-seq(min(object$lam),max(object$lam),by=.01)
  pp<-predict(mm,newdata=lseq,se=T)
  emin<-min(pp$fit)
  llow<-pp$fit-pp$se
  lopt<-min(lseq[llow<=emin])

#  links <- (object$links[,,1]+ object$links[,,2])/2
#  cat("findoptimvalid.pred: links: ", object$links,"\n")
#  cat("findoptimvalid.pred: lambdas: ",object$lambdas,"\n")
#  cat("findoptimvalid.pred: B: ",object$B,"\n")
  mm2<-loess(log(as.vector(t(object$links)))~rep(as.vector(object$lambdas),object$B),degree=deg)
  sopt<-exp(predict(mm2,newdata=lopt))
  o2<-order(mm2$x)
  return(list(sopt=sopt,lopt=lopt,lseq=lseq,pp=pp,mm=mm,o=o,o2=o2,mm2=mm2,llow=llow,emin=emin))
}
plot.EPoC.validation.W <- function(x, ...) {
  o <- findoptimvalid(x)
  plot(o$lseq,o$pp$fit+o$pp$se,lty=2,type='l',xlab="lambda",ylab="W",ylim=c(min(o$pp$fit-o$pp$se),max(o$pp$fit+o$pp$se)), ...)
  lines(o$lseq,o$pp$fit, ...) 
  abline(h=o$wmax,lty=3, ...)
  lines(o$lseq,o$pp$fit-o$pp$se,lty=2, ...)
  points(o$mm$x[o$o],o$mm$fit[o$o],pch=2, ...)
  abline(v=o$lopt.lo,lty=4, ...)
#  abline(v=o$lopt.hi,lty=4)

  axis(side=3,o$mm2$x,labels=round(exp(o$mm2$fit),0))
  text(o$lopt.lo,quantile(o$mm$fit,.3),paste('lambda* =',round(o$lopt.lo,4)))
  #text(o$lopt.lo,quantile(o$mm$fit,.3),paste('lam* in (',round(o$lopt.lo,4),',',round(o$lopt.hi,4),')',sep=''))
  #text(o$lopt.lo,quantile(o$mm$fit,.4),paste('s* in (',round(o$sopt.lo,0),',',round(o$sopt.hi,0),')',sep=''))
  text(o$lopt.lo,quantile(o$mm$fit,.4),paste('s* =',round(o$sopt.hi,0)))
  text(0.2, max(o$pp$fit+0.6*o$pp$se), paste("s: network size"))
}
plot.EPoC.validation.pred <- function(x, ...) {
  o <- findoptimvalid.pred(x)
  plot(o$lseq,o$pp$fit,type='l',xlab="lambda",ylab="CV error", ...) #,ylim=c(.9*min(o$pp$fit),1.1*max(o$pp$fit)))
  abline(h=o$emin,lty=3, ...)
  lines(o$lseq,o$pp$fit+o$pp$se,lty=2, ...)
  lines(o$lseq,o$pp$fit-o$pp$se,lty=2, ...)
  points(o$mm$x[o$o],o$mm$fit[o$o],pch=2, ...)
  abline(v=o$lopt,lty=4, ...)

  xtext <- min(o$mm$x[o$o])
  xtext <- max(xtext, o$lopt)

  axis(side=3,o$mm2$x,labels=round(exp(o$mm2$fit),0))
  text(xtext,quantile(o$mm$fit,.5),paste('lam*=',round(o$lopt,4)),pos=4)
  text(xtext,quantile(o$mm$fit,.6),paste('s*=',round(o$sopt,0)),pos=4)
  text(0.5, max(o$pp$fit+0.5*o$pp$se), paste("s: network size"),pos=1)
}
modelselPlot <- function (x, layout=NULL, k = 1, showtitle=F, bthr=0, showself=F, type=c('graph','modelsel'),...) {
  eps <- 0.1
  newCp <- (x$Cp- min(x$Cp))/(eps + max(x$Cp) - min(x$Cp))
  #newCp <- newCp - min(newCp)
  newBIC <- (x$BIC- min(x$BIC))/(eps + max(x$BIC) - min(x$BIC))
 # newBIC <- newBIC - min(newBIC)
#  print(newCp)
#  print(newBIC)
  plot(x$lambdas,newCp,type='l',ylim=c(min(newCp,newBIC),max(newCp,newBIC)),xlab='lambdas',ylab='Selection criterion', ...)
  lines(x$lambdas,newBIC,lty=2, ...)
  axis(side=3,x$lambdas,x$links, ...)
  loptc<-x$lambdas[which.min(newCp)]
  loptb<-x$lambdas[which.min(newBIC)]
  soptc<-x$links[which.min(newCp)]
  soptb<-x$links[which.min(newBIC)]
  abline(v=loptc,lty=4, ...)
  abline(v=loptb,lty=3, ...)
  #text(x$loptb,quantile(x$mm$Cp,.5),paste('lam*=',round(x$lopt,4)))
  #text(x$loptc,,quantile(x$mm$Cp,.6),paste('s*=',round(x$sopt,0)))
  text(.5, max(newBIC), paste("s: network size"), pos=1, ...)
  legend(quantile(x$lambdas,.8),quantile(newBIC,.75),c('Cp','BIC'),lty=1:2, ...)
}
plot.bootsize <- function(x, lambda.boot=NULL, B=NULL, range=c(0,1), ...) {
  if (is.null(B)) { B <- attr(x,'nboots') }
  G.boot <- x
  sz = 1/B
  bvec<-seq(range[1],range[2], by=sz)
  f <- function(b,G) links<-sum(b<=(nodiag(G)))
  hG<-lapply(G.boot, function(G) unlist(lapply(bvec,f,G=G)))
  K <- length(hG)
  ymax <- max(unlist(hG))
  ymin <- min(unlist(hG))
#  recover()
#  plot(bvec,hG[[1]],)
  plot(c(range[1],range[2]),c(ymin,ymax),log='y',type='n',xlab='bootstrap threshold',ylab='network size',ylim=c(ymin,ymax),xlim=range, ...)
  for(k in 1:K) {
    lines(bvec,hG[[k]],lty=k)
  }
  if (is.null(lambda.boot)) {
    lambda.boot <- attr(x,'lambdas')
  }
  if (!is.null(lambda.boot)) {
#    legend(mean(range), 0.9*ymax,paste('lambda:',round(lambda.boot,3)),lty=1:K, ...)
    legendtexts <- paste('lambda:',round(lambda.boot,3))
  } else {
    # FIXME: get default lambdas from epoc functions, make epoc.bootstrap give them back
    message("Currently you have to provide the lambda values in the lambda.boot parameter to get them in the legend of the plot. Giving indices instead.")
    legendtexts <- paste('lambda idx:',1:K)
#    legend(mean(range), 0.9*ymax,paste('lambda idx:',1:K),lty=1:K, ...)
  }
  legend(mean(range), 0.9*ymax,legendtexts,lty=1:K, ...)
}
epoc.svd <- function(model, k=1, C=1,numload=NULL) {
  #require(survival)
  #require(elasticnet)
  #require(Matrix)
  #require(methods)
  # NB: this method is working with a transposed matrix
  if (is.null(numload)) {
    numload = rep(10,C)
  }
  if (class(model)=="EPOCA" | class(model) == "EPOCG") {
    p <- dim(model$coefficients[[k]])[1]
    m <- t(coef(model,k))
    type = class(model)
  } else if (inherits(model,"Matrix") | class(model) == 'matrix') {
    p <- dim(model)[1]
    m <- t(model)
    isMat = T
    type = "matrix"
  } else {
    stop("epoc.svd needs an EPOC object or a Matrix")
  }
  if (T) {
    d <- diag(m)
    rs <- apply(abs(m - diag(d)), 1, sum)
    cs <- apply(abs(m - diag(d)), 2, sum)
    nz <- (rs + cs)
    ii <- (1:p)[nz>0]
    newp <- length(ii)
    if (newp > 500) {
      warning("Matrix is very large, you may run out of memory")
    } else if (newp == 0) {
      stop("Matrix without diagonal was empty")
    }
    g <- Matrix(m[ii,ii])
  } else {
    newp <- p
    g <- m
    ii <- 1:p
  }
  sg <- svd(g)
  rownames(sg$u) <- rownames(g)
  rownames(sg$v) <- rownames(g)
  ss.in<-spca(t(g)%*%g,type='Gram',K=C,sparse="varnum",para=numload)$load
  ss.out<-spca((g)%*%t(g),type='Gram',K=C,sparse="varnum",para=numload)$load
  return (list(m=m,spload.in=ss.in,spload.out=ss.out,load.in=sg$v[,1:C],load.out=sg$u[,1:C],ii=ii, type=type))
}
epoc.survival <- function(G.svd, Y, U, surv, C=1, type=NULL) {
  #require(survival)
  #require(elasticnet)
  if (G.svd$type=="EPOCA") {
    input <- Y
    output <- U
  } else if (G.svd$type == "EPOCG"){
    input <- U
    output <- Y
  } else if (G.svd$type == "matrix") {
    if (is.null(type)) 
      stop("epoc.survival with a Matrix object need to know whether it is G or A, set type accordingly")
    if (type == 'A') {
      input <- Y
      output <- U
    } else {
      input <- U
      output <- Y
    }
  } else {
    stop("epoc.survival needs an epoc.svd object")
  }
  if (C > dim(G.svd$spload.out)[2]) stop("epoc.survival: requested C exceeds number of SVD components")
  input2 <- input[,G.svd$ii]
  output2 <- output[,G.svd$ii]
  sc.out <-output2%*%G.svd$spload.out[,C]
  sc.in <-input2%*%G.svd$spload.in[,C]
  sdin<-survdiff(Surv(surv)~sign(sc.in))
  sufin<-survfit(Surv(surv)~sign(sc.in))
  sdout<-survdiff(Surv(surv)~sign(sc.out))
  sufout<-survfit(Surv(surv)~sign(sc.out))
  o <- list(survtest.in=sdin,survtest.out=sdout,survfit.in=sufin,survfit.out=sufout,C=C)
  class(o) <- "EPoC.survival"
  return (o)
}
plot.EPoC.survival <- function (x,...) {
  cl <- match.call()
  if (!is.null(cl[['gray']]))
    colors <- 1:1
  else
    colors <- 1:2
  par(mfrow=c(1,2), ...)
  plot(x$survfit.in,xlab='Time',ylab='Survival', main = paste('input scoring comp C=',x$C,sep=''),col=colors, lty=1:2)
  plot(x$survfit.out,xlab='Time',ylab='Survival', main = paste('output scoring comp C=',x$C,sep=''), col=colors, lty=1:2)
}
summary.EPoC.survival <- function(object,...) {
  o <- list(inp=object$survtest.in,outp=object$survtest.out)
  class(o) <- "summary.EPoC.survival"
  o
}
print.summary.EPoC.survival <- function(x,...) {
  cat("\n           In\n")
  print(x$inp, ...)
  cat("\n           Out\n")
  print(x$outp, ...)
}

