BatchSOM <- function(data,grid = somgrid(),min.radius=0.0001,
                     max.radius=0.002,maxit=1000,
                     init=c("random","sample","linear"),
                     radius.type=c("gaussian","bubble","cutgauss","ep"))
{

  set.seed(10)

  if(max.radius<=min.radius)
    stop("max.radius must be larger than min.radius")

  if(min.radius<=0)
    stop("min.radius must be positive")

  if(is.null(init))
    stop("initialization methods is needed")
  initt <- pmatch(init, c("random","sample","linear"))

  if (is.null(radius.type))
    stop("radius type is needed")

  radius.type <- pmatch(radius.type,c("gaussian","bubble","cutgauss","ep"))

  data <- as.matrix(data)
  nd <- nrow(data)
  ng <- nrow(grid$pts)

  xdim<- grid$xdim
  ydim<- grid$ydim


  if(initt == 1){
    init <- matrix(NA, grid$xdim * grid$ydim, dim(data)[2])
    mi <- apply(data, 2, min)
    ma <- apply(data, 2, max)
    for (i in 1:(xdim*ydim)){
      init[i,] <- mi + (ma - mi) * runif(ncol(init))
    }
  }
  else if (initt == 2){
    init <- data[sample(1L:nd, ng, replace = FALSE), , drop = FALSE]
  }
  else {
    ## get the first two principle components
    pcm <- prcomp(data)
    pc <- pcm$rotation[,1:2]
    sd <- pcm$sdev[1:2]
    mn <- apply(data, 2, mean)
    init <- matrix(NA, xdim * ydim, dim(data)[2])
    ## give the 1st pc to the bigger dimension
    if (xdim >= ydim) {
      xtick <- sd[1] * pc[,1]
      ytick <- sd[2] * pc[,2]
    }else {
      xtick <- sd[2] * pc[,2]
      ytick <- sd[1] * pc[,1]
    }
    if (xdim == 1) xis <- rep(0, xdim)
    else xis <- seq(-2, 2, length=xdim)
    if (ydim == 1) yis <- rep(0, ydim)
    else yis <- seq(-2, 2, length=ydim)
    for (i in 1:(xdim*ydim)) {
      xi <- (i - 1) %% xdim + 1
      yi <- (i - 1) %/% xdim + 1
      init[i, ] <- mn + xis[xi] * xtick + yis[yi] * ytick
    }
    init
  }

  nhbrdist <- as.matrix(dist(grid$pts))
  radii<- seq(max.radius,min.radius,len=maxit)
  for(i in 1:maxit)
  {
    cl <- as.numeric(knn1(init, data, 1L:ng))
    if(radius.type == 1)
      A <- exp(-nhbrdist/(2*radii[i]))[,cl]
    if (radius.type == 2)
      A <- (nhbrdist <= radii[i])[,cl]
    if (radius.type == 3)
      A <- exp(-nhbrdist/(2*radii[i])) %*% (nhbrdist <= radii[i])[,cl]
    if (radius.type == 4)
      A <-(1-nhbrdist/radii[i]) %*% (nhbrdist <= radii[i])[,cl]
    ind <- rowSums(A) > 0
    init[ind, ] <- A[ind, ] %*% data / rowSums(A)[ind]
  }

  returnValue <- list(classif= cl, codes = init,grid = grid)
  return(returnValue)
}
