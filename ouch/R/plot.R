tree.plot.internal <- function (x, regimes = NULL, labels = x@nodelabels, legend = TRUE, ...) {
  rx <- range(x@times,na.rm=T)
  rxd <- 0.1*diff(rx)
  anc <- x@anc.numbers
  root <- which(is.root.node(anc))
  if (is.null(regimes)) {
    regimes <- factor(rep('unspec',length(anc)))
    names(regimes) <- x@nodes
  } else if (length(regimes)!=x@nnodes)
    stop(sQuote("regimes")," must be a vector with one entry for each node")
  if ((is.null(names(regimes)))||(!setequal(names(regimes),x@nodes)))
    stop("regime specifications must have names corresponding to the node names")
  regimes <- regimes[x@nodes]
  levs <- levels(as.factor(regimes))
  ## if the root is the only one with a certain regime, toss that regime out
  if (sum(regimes%in%regimes[root])==1)
    levs <- setdiff(levs,regimes[root])
  palette <- rainbow(length(levs))
  for (r in 1:length(levs)) {
    yy <- arrange.tree(root,anc)
    xx <- x@times
    f <- which(!is.root.node(anc) & regimes == levs[r])
    pp <- anc[f]
    X <- array(data=c(xx[f], xx[pp], rep(NA,length(f))),dim=c(length(f),3))
    Y <- array(data=c(yy[f], yy[pp], rep(NA,length(f))),dim=c(length(f),3))
    oz <- array(data=1,dim=c(2,1))
    X <- kronecker(t(X),oz)
    Y <- kronecker(t(Y),oz)
    X <- X[2:length(X)]
    Y <- Y[1:(length(Y)-1)]
    C <- rep(palette[r],length(X))
    if (r > 1) par(new=T)
    par(yaxt='n')
    plot(X,Y,type='l',col=C,xlab='time',ylab='',xlim=rx+c(-rxd,rxd),ylim=c(0,1),...)
    if (!is.null(labels))
      text(X[seq(1,length(X),6)],Y[seq(1,length(Y),6)],labels[f],pos=4,...)
  }
  if (legend)
    legend('topleft',levs,lwd=1,col=palette,bty='n')
  invisible(NULL)
}

arrange.tree <- function (root, anc) {
  k <- which(anc==root)
  n <- length(k)
  reltree <- rep(0,length(anc))
  reltree[root] <- 0.5
  p <- list()
  if (n > 0) {
    m <- rep(0,n)
    for (j in 1:n) {
      p[[j]] <- arrange.tree(k[j],anc)
      m[j] <- length(which(p[[j]] != 0))
    }
    cm <- c(0,cumsum(m))
    for (j in 1:n) {
      reltree <- reltree + (cm[j]/sum(m))*(p[[j]] != 0) + (m[j]/sum(m))*p[[j]]
    }
  }
  reltree
}

setMethod(
          'plot',
          'ouchtree',
          function (x, y, regimes = NULL, node.names = FALSE, legend = TRUE, ..., labels) {
            labels <- x@nodelabels
            if (node.names) {
              lbld <- !is.na(labels)
              labels[lbld] <- paste(x@nodes[lbld],labels[lbld])
              labels[!lbld] <- x@nodes[!lbld]
            }
            if (is.data.frame(regimes)) {
              nm <- rownames(regimes)
              regimes <- lapply(as.list(regimes),function(x){names(x)<-nm;x})
            }
            if (is.list(regimes)) {
              if (any(sapply(regimes,length)!=x@nnodes))
                stop("each element in ",sQuote("regimes")," must be a vector with one entry per node of the tree")
            } else if (!is.null(regimes)) {
              if (length(regimes)!=x@nnodes)
                stop("there must be one entry in ",sQuote("regimes")," per node of the tree")
              nm <- deparse(substitute(regimes))[1]
              regimes <- list(regimes)
              names(regimes) <- nm
            }
            if (is.null(regimes)) {
              invisible(tree.plot.internal(x,regimes=NULL,labels=labels,legend=legend,...))
            } else {
              oldpar <- par(mfrow=c(1,length(regimes)))
              on.exit(par(oldpar))
              retval <- lapply(
                               regimes,
                               function (r) tree.plot.internal(x,regimes=r,labels=labels,...)
                               )
              invisible(retval)
            }
          }
          )

setMethod(
          'plot',
          'hansentree',
          function (x, y, regimes, ...) {
            if (missing(regimes))
              regimes <- x@regimes
            plot(x=as(x,'ouchtree'),y=y,regimes=regimes,...)
          }
          )
