# Copyright (C) 2008-2009 - INRIA - Michael Baudin
# Copyright (C) 2009-2010 - DIGITEO - Michael Baudin
# Copyright (C) 2010-2014 - Sebastien Bihorel
#
# This file must be used under the terms of the CeCILL.
# This source file is licensed as described in the file COPYING, which
# you should have received as part of this distribution. The terms
# are also available at
# http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt
#
# This source code is a R port of the optimsimplex component
# originally written by Michael Baudin for Scilab.

optimsimplex <- function(coords=NULL,fun=NULL,data=NULL,method=NULL,
                         x0=NULL,len=NULL,deltausual=NULL,deltazero=NULL,
                         boundsmax=NULL,boundsmin=NULL,nbve=NULL,
                         simplex0=NULL){

  newobj <- simplex()

  if (!is.null(method)){
    if(!any(method==c('axes','spendley','pfeffer','randbounds','oriented')))
      stop(sprintf('optimsimplex: Unexpected key %s',method),call.=FALSE)

    if (method=='axes'){
      #
      # optimsimplex.axes
      #
      if (size(x0,1)!=1)
        stop(sprintf('optimsimplex: The x0 vector is expected to be a row matrix, but current shape is %d x %d.',
                     size(x0,1),size(x0,2)),
             call.=FALSE)
      if (!is.null(fun))
        assert.classfunction(var=fun,varname='fun',ivar=1)
      if (is.null(len)){
        len <- 1
      }else{
        assert.classreal(var=len,varname='len',ivar=1)
        if (size(len,1)!=1)
          stop(sprintf('optimsimplex: The len vector is expected to be a row matrix, but current shape is %d x %d',
                       size(len,1),size(len,2)),
               call.=FALSE)
      }
      assert.classreal(var=x0,varname='x0',ivar=3)
      n <- length(x0)
      newobj$n <- n
      newobj$nbve <- n+1
      nl <- length(len)
      if (nl==1){
        xlen <- rep(len,n)
      } else {
        # Use the length which has been given, but check that the shape of the matrix is consistent with x0
        if (size(len,2) != size (x0,2))
          stop(sprintf(paste('optimsimplex: The len vector is not consistent with the x0 point. ',
                             'Current shape of x0 is %d x %d while current shape of len is %d x %d.',sep=''),
                       size(x0,1),size(x0,2),size(len,1),size(len,2)),
               call.=FALSE)
        xlen <- len
      }
      newobj$x <- matrix(0,nrow=newobj$nbve,ncol=n)
      newobj$fv <- matrix(0,nrow=newobj$nbve,ncol=1)
    
      #
      # Set all points
      #
      nv <- newobj$nbve
      newobj$x[1:nv,] <- matrix(rep(x0[1:n],nv),nrow=nv,byrow=TRUE)
      newobj$x[2:nv,] <- newobj$x[2:nv,] + diag(xlen)
    
      #
      # Compute Function Value
      #
      if (!is.null(fun)){
        tmp <- optimsimplex.computefv(this=newobj,fun=fun,data=data)
          newobj <- tmp$this
          if (!is.null(data)) data <- tmp$data
      }
    }  
    
    if (method=='spendley'){
      # 
      #optimsimplex.spendley
      #
      if (size(x0,1)!=1)
        stop(sprintf('optimsimplex: The x0 vector is expected to be a row matrix, but current shape is %d x %d',
                     size(x0,1),size(x0,2)),
             call.=FALSE)
      if (!is.null(fun))
        assert.classfunction(var=fun,varname='fun',ivar=2)
      if (is.null(len)){
       len <- 1
      }else{
        assert.classreal(var=len,varname='len',ivar=3)
        if (size(len,1)!=1 | size(len,2)!=1 )
          stop(sprintf('optimsimplex: The len vector is expected to be a row matrix, but current shape is %d x %d',
                       size(len,1),size(len,2)),
               call.=FALSE)
      }
      assert.classreal(var=x0,varname='x0',ivar=1)
      n <- length(x0)
      newobj$n <- n
      newobj$nbve <- n + 1
      newobj$x <- matrix(0,nrow=newobj$nbve,ncol=n)
      newobj$fv <- matrix(0,nrow=newobj$nbve,ncol=1)
    
      #
      # Compute p (diagonal term) , q (off-diagonal term)
      #
      p  <- (n - 1.0 + sqrt(n + 1))/(n * sqrt(2.0))
      q <- (sqrt(n + 1) - 1.0)/(n * sqrt(2.0))
    
      #
      # Set all points
      #
      nv <- newobj$nbve
      newobj$x[1:nv,] <- matrix(rep(x0[1:n],nv),nrow=nv,byrow=TRUE)
      newobj$x[2:nv,] <- newobj$x[2:nv,,drop=FALSE] + diag(rep(p,n)) + 
                         matrix(q,nrow=n,ncol=n) - diag(rep(q,n))
    
      # Compute Function Value
      if (!is.null(fun)){
        tmp <- optimsimplex.computefv(this=newobj,fun=fun,data=data)
          newobj <- tmp$this
          if (!is.null(data)) data <- tmp$data
      } 
    }
    
    if (method=='pfeffer'){
      #
      # optimsimplex.pfeffer
      #
      if (size(x0,1)!=1)
        stop(sprintf('optimsimplex: The x0 vector is expected to be a row matrix, but current shape is %d x %d',
                     size(x0,1),size(x0,2)),
             call.=FALSE)
      if (!is.null(fun))
        assert.classfunction(var=fun,varname='fun',ivar=2)
      if (is.null(deltausual)){
        deltausual <- 0.05
      }
      if (is.null(deltazero)){
        deltazero <- 0.0075
      }
      assert.classreal(var=x0,varname='x0',ivar=1)
      assert.classreal(var=deltausual,varname='deltausual',ivar=3)
      assert.classreal(var=deltazero,varname='deltazero',ivar=4)
    
      n <- length(x0)
      newobj$n <- n
      newobj$nbve <- n + 1
      newobj$x <- matrix(0,nrow=newobj$nbve,ncol=n)
      newobj$fv <- matrix(0,nrow=newobj$nbve,ncol=1)
    
      #
      # Set all points
      #
      newobj$x[,1:n] <- matrix(rep(x0[1:n],newobj$nbve),nrow=newobj$nbve,byrow=TRUE)
    
      #
      # Set points #2 to #n+1
      #
      for (j in 2:(newobj$n+1)){
        if (x0[j-1]==0.0){
          newobj$x[j,j-1] <- deltazero
        }else{
          newobj$x[j,j-1] <- newobj$x[j,j-1] + deltausual * x0[j-1]
        }
      }
      # Compute Function Value
      if (!is.null(fun)){
        tmp <- optimsimplex.computefv(this=newobj,fun=fun,data=data)
          newobj <- tmp$this
          if (!is.null(data)) data <- tmp$data
      } 
    }
    
    if (method=='randbounds'){
      #
      # optimsimplex.randbounds
      #
      if (size(x0,1)!=1)
        stop(sprintf('optimsimplex: The x0 vector is expected to be a row matrix, but current shape is %d x %d.',
                     size(x0,1),size(x0,2)),
             call.=FALSE)
      if (size(boundsmin,1)!=1)
        stop(sprintf('optimsimplex: The boundsmin vector is expected to be a row matrix, but current shape is %d x %d.',
                     size(boundsmin,1),size(boundsmin,2)),
             call.=FALSE)
      if (length(boundsmin)<length(x0))
        stop(sprintf('optimsimplex: The boundsmin vector is expected to have %d columns, but current shape is %d x %d.',
                     size(boundsmin,1),size(boundsmin,2)),
             call.=FALSE)
      if (size(boundsmax,1)!=1)
        stop(sprintf('optimsimplex: The boundsmax vector is expected to be a row matrix, but current shape %d x %d.',
                   size(boundsmax,1),size(boundsmax,2)),
             call.=FALSE)
      if (length(boundsmax)<length(x0))
        stop(sprintf('optimsimplex: The boundsmax vector is expected to have at least as many elements as the x0 vector, but current length is %d.',
                     size(boundsmax,1),size(boundsmax,2)),
             call.=FALSE)
    
      assert.classreal(var=x0,varname='x0',ivar=1)
      assert.classfunction(var=fun,varname='fun',ivar=2)
      assert.classreal(var=boundsmin,varname='boundsmin',ivar=3)
      assert.classreal(var=boundsmax,varname='boundsmax',ivar=4)
      n <- length (x0)
      if (is.null(nbve)){
        nbve <- n + 1
      }else{
        assert.classreal(var=nbve,varname='nbve',ivar=5)
      }
      newobj$n <- n
      newobj$nbve <- nbve
      newobj$x <- matrix(0,nrow=nbve,ncol=n)
      newobj$fv <- matrix(0,nrow=nbve,ncol=1)
    
      #
      # Set all points
      #
      newobj$x[1,1:n] <- x0[1:n,drop=FALSE]
    
      #
      # Set points #2 to #nbve, by randomizing the bounds
      #
      bminmat <- matrix(rep(boundsmin[1:n],nbve-1),nrow=nbve-1,byrow=TRUE)
      bmaxmat <- matrix(rep(boundsmax[1:n],nbve-1),nrow=nbve-1,byrow=TRUE)
      thetas <- matrix(runif(n*(nbve-1)),nrow=n,ncol=(nbve-1))
      newobj$x[2:nbve,1:n] <- bminmat + transpose(thetas) * (bmaxmat - bminmat)
    
      # Compute Function Value
      if (!is.null(fun)){
        tmp <- optimsimplex.computefv(this=newobj,fun=fun,data=data)
          newobj <- tmp$this
          if (!is.null(data)) data <- tmp$data
      }
    }
    
    if (method=='oriented'){
      #
      # optimsimplex.oriented
      #
      if (simplex0$nbve != simplex0$n+1)
        stop(sprintf(paste('optimsimplex: The oriented simplex can be computed only with a simplex',
                           'made of n+1 points, but the dimension is %d and the number of vertices is %d.',sep=''),
                     simplex0$n,simplex0$nbve),
             call.=FALSE)
    
      if (!is.null(fun))
        assert.classfunction(var=fun,varname='fun',ivar=2)
      tmp <- optimsimplex.gradientfv(this=simplex0)
      sgrad <- tmp$g
      rm(tmp)
      ssize <- optimsimplex.size(this=simplex0,method='sigmaminus')
      n <- simplex0$n
    
      # Compute the betas
      ipos <- which(sgrad >= 0.0)
      ineg <- which(sgrad < 0.0)
      betav <- c()
      betav[ipos] <- ssize
      betav[ineg] <- -ssize
      betav <- -0.5 * betav
    
      # Prepare a matrix with beta as diagonal terms
      mid <- diag(betav)
    
      # Compute simplex
      newobj <- simplex()
      newobj$n <- simplex0$n
      newobj$nbve <- simplex0$n+1
      newobj$x <- matrix(0,nrow=n+1,ncol=n)
      newobj$fv <- matrix(0,nrow=n+1,ncol=1)
    
      # Store all points
      x1 <- simplex0$x[1,1:n,drop=FALSE]
      newobj$x[1:(n+1),1:n] <- matrix(rep(x1,n+1),nrow=n+1,byrow=TRUE)
    
      # Retrieve the function value for the first simplex
      # This saves one function evaluation
      newobj$fv[1] <- simplex0$fv[1]
      newobj$x[2:(n+1),1:n] <- mid[1:n,1:n,drop=FALSE] + newobj$x[2:(n+1),1:n,drop=FALSE]
    
      # Compute Function Value
      if (!is.null(fun)){
        tmp <- optimsimplex.computefv(this=newobj,fun=fun,data=data)
          newobj <- tmp$this
          if (!is.null(data)) data <- tmp$data
      }
    }
      
  } else {
    #
    # optimsimplex.coords
    #
    if (!is.null(coords))
      assert.classreal(var=coords,varname='coords',ivar=1)
  
    if (!is.null(fun))
      assert.classfunction(var=fun,varname='fun',ivar=1)
  
    if (!is.null(coords)){
      nbve <- size(coords,1)
      n <- size(coords,2)
      if (nbve<n+1)
        stop(sprintf('optimsimplex: The numbers of columns of coords is %d but is expected to be at least %d',
                     nbve,n+1),
             call.=FALSE)
      newobj$n <- n
      newobj$nbve <- nbve
      newobj$x <- coords[1:nbve,1:n,drop=FALSE]
      
      if (!is.null(fun)){
        tmp <- optimsimplex.computefv(this=newobj,fun=fun,data=data)
          newobj <- tmp$this
          if (!is.null(data)) data <- tmp$data
      }
    }
  }

  varargout <- list(newobj=newobj,data=data)
  
  class(varargout) <- 'optimsimplex'

  return(varargout)
}

