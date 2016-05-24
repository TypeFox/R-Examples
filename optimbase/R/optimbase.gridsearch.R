# Copyright (C) 2010-2014 - Sebastien Bihorel
#
# This file must be used under the terms of the CeCILL.
# This source file is licensed as described in the file COPYING, which
# you should have received as part of this distribution. The terms
# are also available at
# http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt


optimbase.gridsearch <- function(fun=NULL,x0=NULL,xmin=NULL,xmax=NULL,
      npts=3,alpha=10){
  
  # Get x0, xmin and xmax and change them into column vectors
  if (!is.null(x0)) x0 <- matrix(x0,nrow=prod(size(x0)),ncol=1)
  if (!is.null(xmin)) xmin <- matrix(xmin,nrow=prod(size(xmin)),ncol=1)
  if (!is.null(xmax)) xmax <- matrix(xmax,nrow=prod(size(xmax)),ncol=1)
  
  # Configure a neldermead object
  opt <- optimbase()
  opt <- optimbase.set(this=opt,key='x0',value=x0)
  opt <- optimbase.set(this=opt,key='numberofvariables',value=prod(size(x0)))
  opt <- optimbase.set(this=opt,key='function',value=fun)
  if (!is.null(xmin))
    opt <- optimbase.set(this=opt,key='boundsmin',value=xmin)
  if (!is.null(xmax))
    opt <- optimbase.set(this=opt,key='boundsmax',value=xmax)
  
  # Check bounds
  hasbounds <- optimbase.hasbounds(this=opt)
  if (hasbounds){
    tmp <- optimbase.checkbounds(this=opt)
    if (!tmp$isok)
      stop(sprintf('optimbase.gridsearch: %s',tmp$errmsg),
          call.=FALSE)
    rm(tmp)
  }
  
  # Detect number of nonlinear inequality constraints
  neq <- length(opt$fun(x=x0,index=6)$c)
  if (neq>0)
    opt <- optimbase.set(this=opt,key='nbineqconst',value=neq)
  
  # Check x0 is in bounds
  if (optimbase.isfeasible(this=opt,x=x0)$isfeasible==0)
    stop('optimbase.gridsearch: x0 is not within specified bounds.\n',call.=FALSE)
  
  # Check npts
  if (!is.numeric(npts)){
    stop('optimbase.gridsearch: npts is not numeric.', call.=FALSE)
  } else if (length(npts)>1){
    stop('optimbase.gridsearch: npts contains more than one element.', call.=FALSE)
  } else if (npts<=2){
    stop('optimbase.gridsearch: npts must be greater or equal to 3.', call.=FALSE)
  }
  
  # Determine if alpha is to be used
  if (is.null(xmin) | is.null(xmax)){  
    # Check alpha
    if (any(!is.numeric(alpha))){
      stop('optimbase.gridsearch: alpha is not numeric.', call.=FALSE)
    } else if (any(alpha<=1)){
      stop('optimbase.gridsearch: alpha contains value(s) below 1.', call.=FALSE)
    }
    if (length(alpha)>length(x0)){
      alpha <- alpha[1:length(x0)]
    } else {
      alpha <- rep(alpha, length.out=length(x0))
    }
    gridlim <- cbind(x0/alpha, x0*alpha)
  } else {
    gridlim <- cbind(xmin,xmax)
  }
  
  # Create grid
    # Create base grid
    xgrid <- apply(gridlim, 1,
      function(x,...) seq(x[1],x[2],length.out=(npts-1)),
      npts)
    xgrid <- rbind(xgrid, transpose(x0))
    
    # Identify which columns do not need to be expanded
    iscolfix <- apply(xgrid,2,function(x) all(x[1]==x))
    fixedcols <- which(iscolfix)
    varcols <- which(!iscolfix)
    
    fixedvals <- xgrid[1,iscolfix]
    xgrid <- data.frame(xgrid[,varcols])
    
    # Expand xgrid
    xgrid <- do.call(expand.grid, xgrid)
    fixedgrid <- matrix(rep(fixedvals,each=dim(xgrid)[1]),
                        nrow=dim(xgrid)[1])
    xgrid <- cbind(xgrid,fixedgrid)
    xgrid <- xgrid[,order(c(varcols,fixedcols))]
    names(xgrid) <- paste('x', 1:length(x0), sep='')
    
    # Remove duplicates
    xgrid <- unique(xgrid)
    
  # Compute f at all vector of xgrid
  cat(sprintf('The grid contains %d unique combinations.\n',size(xgrid,1)))
  xgrid <- cbind(1:size(xgrid,1),xgrid)
  
  fgrid <- apply(xgrid, 1,
      function(x,...) {
        cat(sprintf('  Evaluating combination number: %d/%d\n',
                    x[1],size(xgrid,1)))
        x <- x[-1]
        optimbase.function(this=opt,x=cbind(x),index=2)$f
      },
      opt,xgrid)
  
  xgrid <- xgrid[,-1]
    
  feasible <- apply(xgrid, 1,
      function(x,...){
        optimbase.isfeasible(this=opt,x=cbind(x))$isfeasible},
      opt)
  
  grid <- data.frame(xgrid, f=fgrid, feasible=feasible)
  
  # Reorder grid
  grid <- grid[order(-grid$feasible,grid$f),]
  
  return(grid)
  
}
