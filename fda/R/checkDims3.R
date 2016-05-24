checkDims3 <- function(x, y=NULL, xdim=2:3, ydim=2:3, defaultNames='x',
         subset=c('xiny', 'yinx', 'neither'),
         xName=substring(deparse(substitute(x)), 1, 33),
         yName=substring(deparse(substitute(y)), 1, 33) ){
#  checkDims3 coerces arguments 'x' and 'y' to 3-d arrays, compares their 
#  dimensions, subsets one or the other or throws an error or a warning 
#  as seems appropriate.  It returns a list with components 'x' and 'y' 
#  with appropriate dimnames.

###  
###
### 1.  length(xdim) == length(ydim)?
###
###  
  mDims <- length(xdim)
  if(length(ydim) != mDims)
    stop('length(xdim) = ', mDims, ' != ', length(ydim),
         ' = length(ydim)')
###  
###
### 2.  check defaultnames
###
###  
  nDefault <- length(defaultNames)
  {
    if(nDefault<1){
      dNames0 <- NULL
      dNames  <- rep(list(NULL), mDims)
    }
    else {
      dNames0 <- defaultNames[[nDefault]]      
      dNames  <- rep(as.list(defaultNames), length=mDims)
      dNms    <- names(defaultNames)
      if (!is.null(dNms)) {
        if(is.null(dNames0)) dNames0 <- dNms[nDefault]
        for (i in 1:mDims)
          if(is.null(dNames[[i]])) dNames[[i]] <- dNms[i]
      }
    }
  }
###  
###
### 3.  loop
###
###  
  for(id in 1:mDims){
    dNi  <- c(dNames[id], dNames0)
    xNmi <- paste('subscript', id, 'of', xName)
    yNmi <- paste('subscript', id, 'of', yName) 
    out  <- checkDim3(x, y, xdim[id], ydim[id], dNi, subset,
                     xName=xNmi, yName=yNmi)
    x <- out$x
    y <- out$y
  }
###  
###
### 4.  done
###
###  
  out
}

#  -----------------------------------------------------------------

checkDim3 <- function(x, y=NULL, xdim=1, ydim=1, defaultNames='x',
         subset=c('xiny', 'yinx', 'neither'),
         xName=substring(deparse(substitute(x)), 1, 33),
         yName=substring(deparse(substitute(y)), 1, 33) ){
###
###  
### 1.  Check xdim, ydim 
###
###  
  if (xdim>3) stop('Requested matching xdim = ', xdim,
                   ' > 3 = maximum allowed.')
  if (is.null(y)){ 
    y <- x
    yName <- xName
  }
  if (ydim>3) stop('Requested matching ydim = ', ydim,
                   ' > 3 = maximum allowed.')
###
###  
### 2.  ixperm, iyperm 
###
###  
  ixperm <- list(1:3, c(2, 1, 3), c(3, 2, 1))[[xdim]]
  iyperm <- list(1:3, c(2, 1, 3), c(3, 2, 1))[[ydim]]  
###
###  
### 3.  x3 <- aperm, ... 
###
###  
  x3 <- aperm(as.array3(x), ixperm);
  y3 <- aperm(as.array3(y), iyperm)   
###
###  
### 4.  xNames, yNames 
###
###  
  xNames <- dimnames(x3)
  yNames <- dimnames(y3)
  {
    if(is.null(defaultNames))
      dNames <- NULL
    else { 
      dNames <- defaultNames[[1]]
      if(is.null(dNames)){
        if(is.null(names(defaultNames)))
          dNames <- defaultNames[[2]]
        else
          dNames <- names(defaultNames)[1]
      }
    }
  }    
###
###  
### 5.  Do it:  Subset & dimnames 
###
###  
  sbst <- match.arg(subset) 
##
##  5.1.  'xiny'
##  
  if(sbst == 'xiny'){  
    if(!is.null(dNames)){
      if(length(dNames) < dim(x3)[1])
        dNames <- paste(dNames, 1:dim(x3)[1], sep='')
      else
        dNames <- dNames[1:dim(x3)[1]]
    }
    nx3 <- dim(x3)[1]
    if(is.null(xNames)){
      if(nx3>dim(y3)[1])
        stop('Can NOT subset ', yName, ' because dim(x)[xdim=',
             xdim, '] = ', nx3, ' > ', dim(y3)[1],
             ' = dim(y)[ydim=', ydim, ']') 
      y3 <- y3[1:nx3,,, drop=FALSE]
      {
        if(is.null(yNames)){ 
          if(!is.null(dNames)){
            dNm <- rep(dNames, length=nx3) 
            dimnames(x3) <- list(dNm, NULL, NULL)
            dimnames(y3) <- list(dNm, NULL, NULL)
          }
        }
        else
          if(!is.null(yNames[[1]])){
            dimnames(x3) <- list(yNames[[1]][1:nx3], NULL, NULL)
            if(nx3<dim(y3)[1]) 
              warning(xName, ' is smaller than ', yName,
                      ' but has no dimnames while ', yName,
                      ' does.  Using the first ', nx3,
                      ' elements of dimension ', ydim,
                      ' of ', yName) 
          }
      }
    }
    else {
      if(is.null(xNames[[1]])){
        if(nx3>dim(y3)[1])
          stop('Can NOT subset ', yName, ' because dim(x)[xdim=',
               xdim, '] = ', nx3, ' > ', dim(y3)[1],
               ' = dim(y)[ydim=', ydim, ']')
        y3 <- y3[1:nx3,,, drop=FALSE]
        {
          if(is.null(yNames)){ 
            if(!is.null(dNames)){
              dNm <- rep(dNames, length=nx3) 
              dimnames(x3)[[1]] <- dNm
              dimnames(y3) <- list(dNm, NULL, NULL)
            }
          }
          else {
            if(is.null(yNames[[1]])){
              if(!is.null(dNames)){
                dNm <- rep(dNames, length=nx3) 
                dimnames(x3)[[1]] <- dNm
                dimnames(y3)[[1]] <- dNm
              }
            }
            else {
              dimnames(x3)[[1]] <- dimnames(y3)[[1]]
              if(nx3<dim(y3)[1]) 
                warning(xName, ' is smaller than ', yName,
                        ' but has no dimnames while ', yName,
                        ' does.  Using the first ', nx3,
                        ' elements of dimension ', ydim,
                        ' of ', yName)
            } 
          }
        }
      }
      else {
        if(is.null(yNames)){
          y3 <- y3[1:nx3,,, drop=FALSE]
          dimnames(y3) <- list(xNames[[1]], NULL, NULL)
          if(nx3<dim(y3)[1]) 
            warning(xName, ' is smaller than ', yName,
                    ' but has no dimnames, while ', yName,
                    ' does.  Using the first ', nx3,
                    ' elements of dimension ', ydim,
                    ' of ', yName)
        }
        else {
          if(is.null(yNames[[1]])){
            y3 <- y3[1:nx3,,, drop=FALSE]
            dimnames(y3)[[1]] <- xNames[[1]] 
            if(nx3<dim(y3)[1]) 
              warning(xName, ' is smaller than ', yName,
                      ' but has no dimnames, while ', yName,
                      ' does.  Using the first ', nx3,
                      ' elements of dimension ', ydim,
                      ' of ', yName)
          }
          else {        
            xiny <- is.element(xNames[[1]], yNames[[1]])
            if(any(!xiny))
              stop('Can NOT subset ', yName, ' because some dimnames(',
                   xName, ')[[xdim=', xdim,
                   ']] are not found in dimnames(y)[[ydim=',
                   ydim, ']];  the first one is ',
                   xNames[[1]][!xiny][1]) 
            y3 <- y3[xNames[[1]],,, drop=FALSE]
          }
        }
      }
    }
  }
  else {
##
##  5.2.  'yinx'
##  
    if(sbst == 'yinx'){
      if(!is.null(dNames)){
        if(length(dNames) < dim(y3)[1])
          dNames <- paste(dNames, 1:dim(y3)[1], sep='')
        else
          dNames <- dNames[1:dim(y3)[1]]
      }
      ny3 <- dim(y3)[1] 
      if(is.null(yNames)){
        if(ny3>dim(x3)[1])
          stop('Can NOT subset ', xName, ' because dim(y)[ydim=',
               ydim, '] = ', dim(y3)[1], ' > ', dim(x3)[1],
               ' = dim(x)[xdim=', xdim, ']') 
        x3 <- x3[1:ny3,,, drop=FALSE]
        {
          if(is.null(xNames)){ 
            if(!is.null(dNames)){
              dNm <- rep(dNames, length=nx3) 
              dimnames(x3) <- list(dNm, NULL, NULL)
              dimnames(y3) <- list(dNm, NULL, NULL)
            }
          }
          else
            if(!is.null(xNames[[1]])){ 
              dimnames(y3) <- list(xNames[[1]][1:ny3], NULL, NULL)
              if(ny3<dim(x3)[1]) 
                warning(yName, ' is smaller than ', xName,
                        ' but has no dimnames while ', xName,
                        ' does.  Using the first ', ny3,
                        ' elements of dimension ', xdim,
                        ' of ', xName) 
            }
        }
      }
      else {
        if(is.null(yNames[[1]])){
          if(ny3>dim(x3)[1]) 
            stop('Can NOT subset ', xName, ' because dim(y)[ydim=',
                 ydim, '] = ', ny3, ' > ', dim(x3)[1],
                 ' = dim(x)[xdim=', xdim, ']') 
          x3 <- x3[1:ny3,,, drop=FALSE]
          {
            if(is.null(xNames)){
              if(!is.null(dNames)){
                dNm <- rep(dNames, length=ny3) 
                dimnames(y3)[[1]] <- dNm
                dimnames(x3) <- list(dNm, NULL, NULL)
              }
            }
            else {
              if(is.null(xNames[[1]])){
                if(!is.null(dNames)){
                  dNm <- rep(dNames, length=ny3) 
                  dimnames(y3)[[1]] <- dNm
                  dimnames(x3)[[1]] <- dNm
                }
              }
              else {
                dimnames(y3)[[1]] <- dimnames(x3)[[1]]
                if(ny3<dim(x3)[1]) 
                  warning(yName, ' is smaller than ', xName,
                          ' but has no dimnames while ', xName,
                          ' does.  Using the first ', ny3,
                          ' elements of dimension ', xdim,
                          ' of ', xName)
              } 
            }
          }
        }
        else {
          if(is.null(xNames)){
            x3 <- x3[1:ny3,,, drop=FALSE]
            dimnames(x3) <- list(yNames[[1]], NULL, NULL)
            if(ny3<dim(x3)[1])
              warning(yName, ' is smaller than ', xName,
                      ' but has no dimnames, while ', xName,
                      ' does.  Using the first ', ny3,
                      ' elements of dimension ', xdim,
                      ' of ', xName)
          }
          else{
            if(is.null(xNames[[1]])){
              x3 <- x3[1:ny3,,, drop=FALSE]
              dimnames(x3)[[1]] <- yNames[[1]]
              if(ny3<dim(x3)[1])
                warning(yName, ' is smaller than ', xName,
                        ' but has no dimnames, while ', xName,
                        ' does.  Using the first ', ny3,
                        ' elements of dimension ', xdim,
                        ' of ', xName)
            }
            else {           
              yinx <- is.element(yNames[[1]], xNames[[1]]) 
              if(any(!yinx))
                stop('Can NOT subset ', xName, ' because some dimnames(', 
                     yName, ')[[ydim=', ydim,
                     ']] are not found in dimnames(x)[[xdim=',
                     xdim, ']];  the first one is ',
                     yNames[[1]][!yinx][1])
              x3 <- x3[yNames[[1]],,, drop=FALSE]
            }
          }
        }
      }
    }
    else
##
##  5.3.  'neither'
##      
      if(sbst == 'neither'){
        if(dim(x3)[1] != dim(y3)[1])
          stop('dim(x)[xdim=', xdim, '] = ', dim(x3)[1],
               ' != ', dim(y3)[1], ' = dim(y)[ydim=',
               ydim, ']')
        if(is.null(xNames)){
          if(!is.null(yNames) &&
             !is.null(yNames[[1]]))
            stop('is.null(dimnames(x)) but ',
                 '!is.null(dimnames(y)[[ydim=', ydim, ']]')
        }
        else{
          if(is.null(xNames[[1]])){ 
            if(!is.null(yNames) &&
               !is.null(yNames[[1]]))
              stop('is.null(dimnames(x)[[xdim=', xdim, ']]), but ', 
                   '!is.null(dimnames(y)[[ydim=', ydim, ']]')
          }
          else{
            if(is.null(yNames))
              stop('x has dimnames;  y has none.') 
            xiny <- (xNames[[1]] %in% yNames[[1]])
            if(any(!xiny))
              stop('Some dimnames(x)[[xdim=', xdim,
                   ']] are not in dimnames(y)[[ydim=',
                   ydim, ']]')
            yinx <- (yNames[[1]] %in% xNames[[1]])
            if(any(!yinx))
              stop('Some dimnames(y)[[ydim=', ydim,
                   ']] are not in dimnames(x)[[xdim=',
                   xdim, ']]')
          }
        }
      }
  }
###
###  
### 6.  out = list(x=aperm, ... ) 
###
###  
  list(x=aperm(x3, ixperm), y=aperm(y3, iyperm) )
}
