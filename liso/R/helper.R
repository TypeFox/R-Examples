
#some helper functions
cdiag = function(x = 1,nrow,ncol){
	diag(as.vector(x) -> y,nrow = length(y))
}

##NOTES: defining a new sparse format for step functions in many dims
#
#class <- "multistep"
#$params 
#index   number of breaks    intercept
#$values 
#coefchaincount breakpoint delta previous
#$intercept 
#$ranges
#min max


# * -> used to get Y's out of it
"*.multistep" <- function(e1, e2){
  if (!inherits(e2,"multistep")){
    return ("*.multistep"(e2,e1))
  }
  pp <- nrow(e2$params)
  if (is.null(dim(e1))){
    e1 <- matrix(e1,ncol = pp)
  }
  nn <- nrow(e1)
  addup <- function(vect){
    i = vect[1]
    if (e2$params[i,2] == 0){
      return(rep(0,length(vect)-1))
    } else {
      matchtable = e2$values[e2$params[i,1]:(e2$params[i,1]+e2$params[i,2]-1),2:3,drop = FALSE]
      res = sapply(vect[-1], function(x)  sum(matchtable[matchtable[,1] <= x,2])) + e2$params[i,3]
    }
    res
  }
  results = rowSums(matrix(apply(rbind(1:pp,e1),2,function(x) addup(x)),nrow=nn)) + e2$intercept
  results
}

predict.multistep<- function(object, newx, ...){
  object * newx
}

#the same as *, but does linear interpolation
"&.multistep" <- function(x, obj){
  if (!inherits(obj,"multistep")){
    return ("&.multistep"(obj,x))
  }
  pp <- nrow(obj$params)
  if (is.null(dim(x))){
    x <- matrix(x,ncol = pp)
  }
  nn <- nrow(x)
  addup <- function(vect){
    i = vect[1]
    if (obj$params[i,2] == 0){
      return(rep(0,length(vect)-1))
    } else {
      matchtable = rbind(c(obj$range[i,1],0),obj$values[obj$params[i,1]:(obj$params[i,1]+obj$params[i,2]-1),2:3,drop = FALSE])
      matchtable[,2] = cumsum(matchtable[,2])
      res = approx(matchtable, xout=vect[-1],rule=2)$y + obj$params[i,3]
    }
    res
  }
  results = rowSums(matrix(apply(rbind(1:pp,x),2,function(a) addup(a)),nrow=nn)) + obj$intercept
  results
}

summary.multistep <- function(object, ...){
list(max = max(object), min = min(object), totalvar=abs(object), dim=dim(object))
}

max.multistep <- function(x,..., na.rm){
  individualmax <- function(vect){
    if (vect[2] == 0){
      return(0)
    } else {
      return( max(cumsum(c(vect[3],x$values[vect[1]:(vect[1]+vect[2]-1), 3]))))
    }
  }
  res <- drop(apply(x$params,1,individualmax))
  res
}
min.multistep <- function(x, ...,na.rm){
  individualmin <- function(vect){
    if (vect[2] == 0){
      return(0)
    } else {
      return( min(cumsum(c(vect[3],x$values[vect[1]:(vect[1]+vect[2]-1), 3]))))
    }
  }
 
  res <- drop(apply(x$params,1,individualmin))
  res
}
abs.multistep <- function(x){
  individualtv <- function(vect){
    if (vect[2] == 0){
      return(0)
    } else {
      return( sum(abs(x$values[vect[1]:(vect[1]+vect[2]-1), 3])) )
    }
  }
  res <- drop(apply(x$params,1,individualtv))
}

dim.multistep <- function(x){
  res <- drop(x$params[,2])
}

#pull out a coefchain
as.vector.multistep <- function(x, mode){
  coefchain = replace(rep(0, (x$n - 1)*nrow(x$params)), x$values[,1], x$values[,3])
  coefchain
}

#turns multistep into 'fits'
as.array.multistep <- function( x, dims = 1:nrow(x$param), explode = FALSE, signs = NULL, ...){
  coefchain = replace(rep(0, (x$n - 1)*nrow(x$params)), x$values[,1], x$values[,3])
  if (!explode){
  result = apply(rbind(x$params[dims,3],matrix(coefchain, x$n - 1))[,dims], 2, cumsum)
  } else {
    if (is.null(signs)){
      result = scale(apply(rbind(0, cbind(matrix(pmax(coefchain,0), x$n-1)[,dims],apply(-matrix(pmin(coefchain, 0), x$n-1)[,dims],2,rev))),2,cumsum), center=TRUE,scale=FALSE)
    } else {
      if (length(signs) == 1) signs = rep(signs, length(dims))
      #okay, we should explode or flip the appropriate signs.
      result1 = scale(apply(rbind(0, cbind(matrix(pmax(coefchain,0), x$n-1)[,dims],apply(-matrix(pmin(coefchain, 0), x$n-1)[,dims],2,rev))),2,cumsum), center=TRUE,scale=FALSE)
      result = apply(rbind(x$params[dims,3],matrix(coefchain, x$n - 1))[,dims], 2, cumsum)
      for (i in dims){
        if (signs[i] == -1){
          result[,i] = rev(result[,i])
        } else if (signs[i] == 0){
          result[,i] = result1[,i]
          result = cbind(result, result1[,i + length(dims)])
        }
      }
    }
  }
  result
}


multistep <- function(coefchain,x=NULL,intercept=0,sortedx = apply(x,2,sort),names = NULL, pinters=NULL,...){
  pp <- ncol(sortedx)
  nn <- nrow(sortedx)
  coefchainm <- matrix(coefchain, ncol=pp)
  sparsitycounts <- drop(apply(coefchainm, 2, function(a) sum(a !=0)))
  if (is.null(pinters)) {pinters <- drop(apply(coefchainm, 2, function(a) -sum(cumsum(a))/nn))}
  params <- cbind(c(1,cumsum(sparsitycounts)+1)[1:pp], sparsitycounts, pinters)
  coefchaincounts <- which(coefchain != 0)
  breakpoints <- as.vector(sortedx[-1,])[coefchaincounts]
  previous <- as.vector(sortedx[-nn,])[coefchaincounts]
  delta <- coefchain[coefchaincounts]
  values <- cbind(coefchaincounts, breakpoints, delta, previous)
  if (!is.null(names)){
    dimnames(params)[[1]] = names
  }
  obj <- list(params = params, values = values, intercept = intercept, range = cbind(sortedx[1,] , sortedx[nn,]),n=nn, ...)
  class(obj) <- "multistep"
  obj
}

#combine two or more multistep fits, or scale
#default is average
combmulti <- function(multisteplist, scales = rep(1/length(multisteplist), length(multisteplist)), x = NULL, sortedx=apply(x,2,sort)){
  newcoefchain = apply(sapply(multisteplist,as.vector) ,1,function(a) sum(a * scales))
  newintercept = sum(scales * sapply(multisteplist, function(a) a$intercept))
  newpinters = apply(sapply(multisteplist, function(a) min(a)), 1, function(b) sum(b*scales) )
  tempmultistep = multistep(newcoefchain, intercept=newintercept, sortedx = sortedx, pinters = newpinters, names = dimnames(multisteplist[[1]]$params)[[1]])
#re-intercept to give correct weighted means against the first set of weights

#copy over to get the remaining fields in the original
  result = multisteplist[[1]]
  result$params = tempmultistep$params
  result$values = tempmultistep$values
  result$intercept = tempmultistep$intercept
  result

}
#gets 'row and column' for a value
modrow <- function(n, rows){
  c( (n -1)%/% rows + 1, (n-1) %% rows +1 )
}
  
#add um adding, allow obj = NULL
  #needs a bit smarter code using plot.window. Too busy to do this
plot.multistep <- function(x = NULL, xpoints=NULL, ypoints = NULL, dims = 1:max(nrow(x$param), ncol(xpoints)) , ylimit = cbind(min(min(x),max(x)), max(max(x), min(x))), grid = TRUE, add = FALSE, titles = !add,...){
  numplots = length(dims)
  if (!add){
  if (grid) op <- par(mar = c(2,2,1,0) + 0.1,mfrow = c((ceiling(sqrt(numplots))) -> width,ceiling(numplots/width)))
  } else {
    numberrows = par()$mfrow[2]
    op = par()
  }
  for (idim in 1:length(dims)){
    if (add &grid) par(mfg = modrow(idim, numberrows))
    if (!is.null(x)){
      if (dim(x)[dims[idim]] == 0){
        plotdata = NULL
      } else {
        plotdata = x$values[x$params[dims[idim],1]:(x$params[dims[idim],1] + x$params[dims[idim],2] - 1),2:3, drop=FALSE]
        plotdata[,2] = cumsum(plotdata[,2]) + x$params[dims[idim],3]
      }
      plotdata = rbind(c(x$range[dims[idim],1], x$params[dims[idim], 3]), plotdata)
      plotdata = rbind(plotdata, c(x$range[dims[idim],2], tail(plotdata[,2], n=1)))
      if (!add){
        plot(plotdata,type="s", ann=F,frame.plot=F,axes=T, ylim = ylimit, ...)#,ann = F, axes = F)
      } else {
        lines(plotdata, type="s",...)
      }
    }
    if (titles){
      if (is.null(dimnames(x$params)[[1]][dims[idim]])){
        title(main = dims[idim])
      } else {
        title(main = dimnames(x$params)[[1]][dims[idim]])
      }
    }
    if (!(is.null(xpoints))){
      if (!is.null(x)){
        points(xpoints[,dims[idim]],ypoints,pch=".")
      } else {
        points(xpoints[,dims[idim]],ypoints,pch=".",...)
      }
    }
  }
  if (grid) invisible(op) # return the old pars, don't reset!
}

#ave for sorted factors by Charles Berry
#plus weights
ave.sorted <- function(y, x, weights = NULL){
  reps <- rle(x)$lengths
  lens <- rep(reps,reps)
  uniqLens <- unique(lens)
  for (i in uniqLens[ uniqLens != 1]){
    if (is.null(weights)){
      y[ lens == i] <- rep( colMeans(matrix(y[ lens == i], nr=i)), each=i)
    } else {
      w = matrix(weights[lens == i], nr = i)
      y[ lens == i] <- rep( colSums(matrix(y[ lens == i], nr=i) * w)/colSums(w), each=i)
    }
  }
  y
}

#changes step function to a new set of steps
#rejig <- (stepobj, newx, newY=NULL){
removeZero <- function(obj){
  obj[which(obj!=0)]
}

#create log scale grids
seq.log <- function(from = 1, to = 1, length.out = 50, add.zero = FALSE, shifting = 0, ...){
  res = exp(seq(from = log(from + shifting), to = log(to + shifting), length=length.out - add.zero)) - shifting
  if (add.zero) {
    if (from > to) {
      res = c(res,0)
    } else {
      res = c(0,res)
    }
  }
  res
}

objsort <- function(){
  sort(sapply(ls(1), function(x) object.size(get(x))))
}

#plots a matrix as a set of lines or points with legend
repeatplot <- function(x, y, margin = 1, col = 1, type = "p", form = 1, new = TRUE, lims = TRUE, ...){
  if(missing(y)) {
    y <- x
    x <- 1:dim(x)[3 - margin]
  }
  if (margin != 1){
    x = t(x)
    y = t(y)
  }
  dd = nrow(y)
  nn = ncol(y)
  if (length(col) < dd) col = rep(col,dd)
  if (length(type) < dd) type = rep(type,dd)
  if (length(form) < dd) form = rep(form,dd)
  if (!is.matrix(drop(x))){
    x <- matrix(rep(x, dd), dd, byrow=T)
  }
  if (new){
    if (lims){
      plot(x[1,], y[1,], col = col[1], type= type[1], xlim = range(x), ylim = range(y), ...)
    } else {
      plot(x[1,], y[1,], col = col[1], type= type[1], ...)
    }
  }
  for (i in (new + 1):dd){
    if (type[i] == "l") lines(x[i,], y[i,], col= col[i], lty = form[i])
    if (type[i] == "p") points(x[i,], y[i,], col= col[i], pch = form[i])
  }
}

#constrain x to a given range
constrain <- function(x, minimum, maximum=-minimum){
  if (minimum > maximum) {
    temp <- minimum
    minimum <- maximum
    maximum <- temp
  }
  x = replace(x, which(x > maximum), maximum)
  x = replace(x, which(x < minimum), minimum)
  x
}

importance.liso = function(obj, weights = rep(1/obj$n, obj$n)){
  sqrt(colSums(as.array.multistep(obj)^2 * apply(obj$sorter, 2, function(x) weights[x])))
}

# USAGE: 
# squelch((maxnumofcols), cols, F)[oldindex] will give the new index with cols
# removed
# squelch((oldnumofcols), cols, T)[oldindex] will give the new index with cols
# added (cols are original positions)
squelch <- function(n, cols=NULL, add = F){
  if (add){
    setdiff(1:n, cols)
  } else {
    cumsum(c(1,replace(rep(1, n), cols, 0)))
  }
}

#which min amongst x > eps
which.min.pos <- function(x, eps = 1e-16){
  x[x <= eps] = NA
  which.min(x)
}

#remove all but one zero column or row, for plotting
cleanzero <- function(x, mar=1, eps=1e-15){
  if (mar == 2) return(t(cleanzero(t(x))))

  zerorows = which(rowSums(abs(x)) < eps)
  if (length(zerorows) > 1){
    newx <- x[-zerorows[-1], ]
  }
  attr(newx, "indices") = (1:nrow(x))[-zerorows[-1]]
  newx
}


cv.folds = function(n, folds = 10) 
{
    split(sample(1:n), rep(1:folds, length = n))
}

error.bars = function(x, upper, lower, width = 0.02, ...) 
{
    xlim <- range(x)
    barw <- diff(xlim) * width
    segments(x, upper, x, lower, ...)
    segments(x - barw, upper, x + barw, upper, ...)
    segments(x - barw, lower, x + barw, lower, ...)
    range(upper, lower)
}

