###############################################################
####################  USER FUNCTIONS  #########################
###############################################################

n.nodes = function(object){
  if(class(object) == "CID.Gibbs"){
    object = object$CID.object
  }
  if(class(object) == "CIDnetwork"){
    return(object$n.nodes)
  }else(stop("argument must be of CID.Gibbs or CIDnetwork class!"))
}

edge.list = function(object){
  if(class(object) == "CID.Gibbs"){
    object = object$CID.object
  }
  if(class(object) == "CIDnetwork"){
    return(object$edge.list)
  }else(stop("argument must be of CID.Gibbs or CIDnetwork class!"))
}

outcome = function(object){
  if(class(object) == "CID.Gibbs"){
    object = object$CID.object
  }
  if(class(object) == "CIDnetwork"){
    return(object$outcome)
  }else(stop("argument must be of CID.Gibbs or CIDnetwork class!"))
}

is.net.directed = function(object){
  if(class(object) == "CID.Gibbs"){
    object = object$CID.object
  }
  if(class(object) == "CIDnetwork"){
    return(object$is.directed)
  }else(stop("argument must be of CID.Gibbs or CIDnetwork class!"))
}

net.density = function(object){
  if(class(object) == "CID.Gibbs"){
    object = object$CID.object
  }
  if(class(object) == "CIDnetwork"){
    n = n.nodes(object)
    yy = outcome(object)
    den = sum(yy[!is.na(yy)])/((n)*(n-1)) #incase there are missing edges
    return(den)
  }else(stop("argument must be of CID.Gibbs or CIDnetwork class!"))
}


node.names = function(object){
  if(class(object) == "CID.Gibbs"){
    object = object$CID.object
  }
  if(class(object) == "CIDnetwork"){
    return(object$node.names)
  }else(stop("argument must be of CID.Gibbs or CIDnetwork class!"))
}


inDegree = function(object){
  if(class(object) == "CID.Gibbs"){
    object = object$CID.object
  }
  if(class(object) == "CIDnetwork"){
    socioMat = socio(object)
    return(colSums(socioMat,na.rm = TRUE))
  }else(stop("argument must be of CID.Gibbs or CIDnetwork class!"))
}

outDegree = function(object){
  if(class(object) == "CID.Gibbs"){
    object = object$CID.object
  }
  if(class(object) == "CIDnetwork"){
    socioMat = socio(object)
    return(rowSums(socioMat,na.rm = TRUE))
  }else(stop("argument must be of CID.Gibbs or CIDnetwork class!"))
}


socio = function(object){
  if(class(object) == "CID.Gibbs"){
    object = object$CID.object
  }
  if(class(object) == "CIDnetwork"){
    n = n.nodes(object)
    socio = array(NA, dim = c(n,n))
    df = data.frame(edge.list(object),outcome(object))
    if(is.net.directed(object)){
      for(i in 1:dim(df)[1]){
        socio[df[i,1],df[i,2]] = df[i,3]
      }
    } else{
      for(i in 1:dim(df)[1]){
        socio[df[i,1],df[i,2]] = socio[df[i,2],df[i,1]] =df[i,3]
      } }
    return(socio)
  }else (stop("argument must be of CID.Gibbs or CIDnetwork class!"))
}


value.mat <- function(CID.Gibbs.object,prob=TRUE){
  summ.mat <- CID.Gibbs.object$CID.object$gibbs.value(CID.Gibbs.object$results)
  vals <- apply(summ.mat,1,mean)
  n.nodes <- CID.Gibbs.object$CID.object$n.nodes
  mat <- array(NA,c(n.nodes,n.nodes))
  edge.list <- CID.Gibbs.object$CID.object$edge.list

  for(rr in 1:length(vals)){
    mat[edge.list[rr,1],edge.list[rr,2]] <- vals[rr]
  }
  if(CID.Gibbs.object$CID.object$class.outcome == "gaussian"){
    return(mat)
  }else{
    if(prob){
      return(pnorm(mat))
    }else{
      return(mat)
    }
  }
}

value.mat.mean <- function(object,prob=TRUE){
  ##    summ.mat <- CID.Gibbs.object$CID.object$gibbs.value(CID.Gibbs.object$results)
  ##    vals <- apply(summ.mat,1,mean)
  if(class(object) == "CID.Gibbs"){
    object = object$CID.mean
  }
  if(class(object) != "CIDnetwork"){
    stop("object must be of typee CID.object or CIDnetwork")
  }
  vals <- object$value()
  n.nodes <- object$n.nodes
  mat <- array(NA,c(n.nodes,n.nodes))
  edge.list <- object$edge.list

  for(rr in 1:length(vals)){
    mat[edge.list[rr,1],edge.list[rr,2]] <- vals[rr]
  }
  if(object$class.outcome == "gaussian"){
    return(mat)
  }else{
    if(prob){
      return(pnorm(mat))
    }else{
      return(mat)
    }
  }

}

log.likelihood <- function(CID.Gibbs.object){
  return(CID.Gibbs.object$CID.mean$log.likelihood)
}

switcheroo <- function(CID.Gibbs.object){
  return(CID.Gibbs.object$CID.object$gibbs.switcheroo(CID.Gibbs.object$results))
}

get.blocks <- function(object){
  if(!is(object,"summary.CID.Gibbs")){
    if(is(object,"CID.Gibbs")){
      object <- summary(object)
    }else{
      stop("object must be of class CID.Gibbs or summary.CID.Gibbs")
    }
  }
  sbm.cid <- object$SBMcid
  if(is.null(sbm.cid)){
    stop("Object contains no SBM component")
  }
  return(sbm.cid$modal.membership)
}



#####  Function Required by netCV
CID.metric <- function(graph,
                       ...){
  fit <- suppressMessages(CID.Gibbs(input=graph,
                                    ...))
  return(value.mat(fit))
}



####################################################################
######################  Plotting Functions  ########################
####################################################################
network.plot <- function (x, fitted.values=FALSE, ...) {

  if (class(x) == "CID.Gibbs") {
    if (fitted.values) {
      values <- apply(x$CID.object$gibbs.value(x$results), 1, mean)  #fitted linear predictor.
    }
    x <- x$CID.object
  }
  if(class(x) == "CIDnetwork"){
    if(!fitted.values){
      values <- x$outcome   #This is the actual data outcome.
    }
    x$plot.network(values, ...)
  }else stop ("Argument 1 is not an object of class CID.Gibbs or CIDnetwork.")
}

sociogram.plot <- function (x, component.color=0, vertexcolor, add.labels=TRUE, ...){
  if(class(x) == "CID.Gibbs"){
    res <- x$results
    switched <- x$CID.object$gibbs.switcheroo(x$results)
    x <- x$CID.object
  }else{
    res <- NULL
  }
  if(class(x) == "CIDnetwork"){

    ##pull out the greater-than-zero edge list.
    edges <- x$edge.list[!is.na(x$outcome) & x$outcome > 0,]

    if(missing(vertexcolor)){
      if (component.color > length(x$components)) stop ("Invalid component number for sociogram plot.")
      if (component.color > 0 && !is.null(res)) {
        vertexcolor <- x$components[[component.color]]$gibbs.node.colors(switched[[component.color+5]])
      } else vertexcolor <- rep("#DDDDFF", x$n.nodes)
    }

    plot(graph.edgelist(edges, directed=x$is.directed),
         vertex.label=x$node.names,
         vertex.color=vertexcolor,
         ...)
  } else stop ("Argument 1 is not an object of class CID.Gibbs.")
}




                                        #sociogram.plot <- function (x, component.color=0, ...)
                                        #  if (class(x) == "CID.Gibbs") {
                                        #
                                        #    #pull out the greater-than-zero edge list.
                                        #  edges <- x$CID.object$edge.list[x$CID.object$outcome > 0,]
                                        #
                                        #  if (component.color > length(x$CID.object$components)) stop ("Invalid component number for sociogram plot.")
                                        #
                                        #  if (component.color > 0) {
                                        #    switched <- x$CID.object$gibbs.switcheroo(x$results)
                                        #    vertexcolor <- x$CID.object$components[[component.color]]$gibbs.node.colors(switched[[component.color+4]])
                                        #  } else vertexcolor <- rep("#DDDDFF", x$CID.object$n.nodes)
                                        #
                                        #  plot(graph.edgelist(edges, directed=FALSE),
                                        #       vertex.label=x$CID.object$node.names,
                                        #       vertex.color=vertexcolor,
                                        #       ...)
                                        #} else stop ("Argument 1 is not an object of class CID.Gibbs.")

##

################################################################
##################  Component Specific Plots  ##################
################################################################


likelihood.plot = function(x, ...){
  if(class(x) != "CID.Gibbs"){stop("argument must be of CID.Gibbs class!")}
  plot(x, which.plot = 3, auto.layout = FALSE, ...)
}


###
intercept.plot = function(x, mode=c("standard","trace"), ...){
  if(class(x) != "CID.Gibbs"){
    stop("argument must be of CID.Gibbs class!")
  }else{
    mode <- match.arg(mode)
    if(mode == "trace"){
      plot(x,which.plot = 1, auto.layout = FALSE)
    }else{
      switched = switcheroo(x)
      intercept = unlist(switched[[1]])
      n = length(intercept)
      plot(intercept,seq(1,10,length = n), type = 'n', main = '95% credible interval of estimated intercept from Gibbs Sampler', ylab = NA,xlab = 'Estimated intercept', yaxt = 'n', ...)
      X0 = mean(intercept)
      q1 = quantile(intercept,0.025)
      q2 = quantile(intercept, 0.975)
      points(X0,5, col = 'red',pch = 20, cex = 3.5)
      arrows( q1, 5, q2, 5,code=3,angle=90,length=0.2,col='red',lwd = 3);
      abline(v = 0, lty = 3)

    }
  }
}


mat.cov.to.edge.list.cov.old <- function(Xmat){
  kk <- nrow(Xmat)
  out.1 <- apply(Xmat,3,FUN=function(x){
    return(as.vector(t(x)))
  }    )
  out.final <- array(NA,c(kk*kk - kk,dim(Xmat)[3]))
  iter <- outer <- 0
  for(ii in 1:kk ){
    for(jj in 1:kk ){
      outer <- outer + 1
      if(ii != jj){
        iter <- iter + 1
        out.final[iter,] <- out.1[outer,]
      }
    }
  }
  return(out.final)
}

mat.cov.to.edge.list.cov <- function (Xmat, n.nodes=dim(Xmat)[1], arc.list=make.arc.list(n.nodes)) {

  if (!(is(Xmat, "array") | is(Xmat, "matrix"))) stop ("mat.cov.to.edge.list.cov: input is not a matrix or array.")
  if (dim(Xmat)[1] != dim(Xmat)[2]) stop ("mat.cov.to.edge.list.cov: input is not square.")

  if (length(dim(Xmat)) == 2) Xmat <- array(c(Xmat), c(dim(Xmat), 1))
  orig.arcs <- make.arc.list(n.nodes)
  covs <- sapply(1:dim(Xmat)[3], function (qq) t(Xmat[,,qq][non.diag(n.nodes)])) #if (symmetric) else sapply(1:dim(Xmat)[3], function (qq) Xmat[,,qq][u.diag(n.nodes)])
  covs <- covs[match(arc.list[,1]*n.nodes+arc.list[,2],
                     orig.arcs[,1]*n.nodes+orig.arcs[,2]),]

  return(covs)

}








####
COV.plot = function(x, mode=c("standard","trace","scatterplot"), ...){
  if(class(x) != "CID.Gibbs"){stop("argument must be of CID.Gibbs class!")}
  switched = switcheroo(x)
  ix = which(names(switched) == "COVout")
  mode <- match.arg(mode)
  if(length(ix) > 0 ){
    if(mode == "standard"){
      plot(x, which.plot = ix, auto.layout = FALSE, ...)
    }else{
      cov.chain <- switched$COVout
      n.draws <- length(cov.chain)
      n.coef <- length(cov.chain[[1]]$coef.cov)
      chain.table <- matrix(unlist(cov.chain),
                            nrow=n.draws,ncol=n.coef,
                            byrow=TRUE)
      p.cols <- ceiling(sqrt(n.coef))
      p.rows <- ceiling(n.coef/p.cols)
      par(mfrow=c(p.rows,p.cols))
      if(mode == "trace"){
        for(ii in 1:n.coef){
          plot(chain.table[,ii],type="l")
        }
      }else if(mode == "scatterplot"){
        pairs(chain.table)
      }
    }
  }else(stop("Wrong model specified"))
}



##
LSM.plot = function(x, ...){
  if(class(x) != "CID.Gibbs"){stop("argument must be of CID.Gibbs class!")}
  switched = switcheroo(x)
  ix = which(names(switched) == "LSMout")
  if(length(ix) > 0){
    plot(x, which.plot = ix, auto.layout = FALSE, ...)
  }else(stop("Wrong model specified"))
}

SBM.plot = function(x, ...){
  if(class(x) != "CID.Gibbs"){
    stop("argument must be of CID.Gibbs class!")
  }
  switched = switcheroo(x)
  ix = which(names(switched) == "SBMout")
  if(length(ix) > 0){
    plot(x, which.plot = ix, auto.layout = FALSE, ...)
  }else(stop("Wrong model specified"))
}

MMSBM.plot = function(x, ...){
  if(class(x) != "CID.Gibbs"){
    stop("argument must be of CID.Gibbs class!")
  }
  switched = switcheroo(x)
  ix = which(names(switched) == "MMSBMout")
  if(length(ix) > 0){
    plot(x, which.plot = ix, auto.layout = FALSE, ...)
  }else(stop("Wrong model specified"))
}



SR.plot = function(x, ...){
  if(class(x) != "CID.Gibbs"){
    stop("argument must be of CID.Gibbs class!")
  }
  switched = switcheroo(x)
  ix = which(names(switched) == "SRout")
  if(length(ix) > 0){
    plot(x, which.plot = ix, auto.layout = FALSE, ...)
  }else(stop("Wrong model specified"))
}




##Function for computing the posterioir mode of the positions
#using 2dimensional KDE

#1. Find np density.
#2. Find the position (x,y) that maximize the density
#required library(MASS)

getModeLS = function(object,gridsize = 50){
    if(!is(object,"CID.Gibbs")){
        stop('Object should be of class CID.Gibbs')
    }
    switched = switcheroo(object)
    if(is.null(switched$LSMout)){
        stop('Wrong model specified')
    }
    draws = length(switched$LSMout)
    dd = dim(switched$LSMout[[1]]$latent.space.pos)[2]
    if(dd > 2){
	warning("Current code finds posterior mode for 2-dimensional latent positions only")}
    nn = n.nodes(object)
    npmode = array(NA,dim=c(nn,dd))
    for(ii in 1:nn){
        XX = sapply(1:draws,function(y){
            switched$LSMout[[y]]$latent.space.pos[ii,1]})
        YY = sapply(1:draws,function(y){
            switched$LSMout[[y]]$latent.space.pos[ii,2]})
        npfit = kde2d(x=XX,y=YY, n = gridsize,
                      lims = c(range(XX),range(YY)))
        yind = sapply(1:gridsize,function(x){
                which(npfit$z[x,]==max(npfit$z[x,]))})
        xmax = sapply(1:gridsize,function(x) (npfit$z[x,yind[x]]))
        xi = which(xmax == max(xmax))
        npmode[ii,] = c(npfit$x[xi],npfit$y[yind[xi]])
    }
    return(npmode)
}


