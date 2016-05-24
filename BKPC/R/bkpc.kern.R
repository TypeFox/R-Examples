bkpc.kern <-
function(x, y, n.kpc = NULL, thin = 100, 
                      n.iter = 100000, std = 10,  g1 = 0.001, 
                      g2 = 0.001, g3 = 1, g4 = 1, initSigmasq = NULL, initBeta = NULL, initTau = NULL, intercept = TRUE, rotate = TRUE, ... ){

  if(!is.factor(y))stop("error: y should be a factor") 

  ntr <- dim(x)[1] 
  n.class <- length(levels(y)) - 1  
  y.ind <- matrix(0, ntr, n.class)      
  for(i in 2:(n.class + 1))y.ind[y == levels(y)[i], i - 1] <- 1
  y.ind <- matrix(y.ind, ntr * n.class, 1)
  
  
  
  if (rotate){ 
    if (is.null(n.kpc)) n.kpc <- dim(x)[2]
    objectKPCA <- kPCA(x)
    design <- objectKPCA$KPCs[ , 1 : n.kpc] 
  }
  else{
    n.kpc <- dim(x)[2]
    design <- x
  }
  
  
  if (intercept){
    design <- cbind(matrix(1, ntr, 1), design)
    n.kpc <- n.kpc + 1
  }
  
  design <- matrix(t(design), 1, n.kpc * ntr)
  
  n.samples <- n.iter/thin  
  
  if(is.null(initSigmasq))initSigmasq  <- runif(1)
  
  if(is.null(initBeta))initBeta  <- matrix(runif(n.kpc * n.class), n.kpc * n.class, 1)
  else initBeta  <- matrix(initBeta, n.kpc * n.class, 1)
  if(is.null(initTau))initTau  <- matrix(runif(n.kpc * n.class), n.kpc * n.class, 1)
  else initTau  <- matrix(initTau, n.kpc * n.class, 1) 
  
  b <- .C("mainbkpc", as.double(design),  as.integer(y.ind), as.integer(n.kpc), 
          as.integer(thin), as.integer(ntr), as.integer(n.iter),
          as.integer(n.class), as.double(std),  as.double(g1), as.double(g2), as.double(g3), as.double(g4), as.double(initSigmasq), as.double(initBeta), as.double(initTau), 
          betas = as.double(matrix(0, n.kpc * n.class * n.samples, 1) ), taus = as.double(matrix(0, n.kpc * n.class * n.samples, 1) ), 
          zs = as.double(matrix(0, ntr * n.class * n.samples, 1) ), sigmasqs = as.double(matrix(0, n.samples, 1) ))      
  
  object <- new.env()

  object$beta <- t(matrix(b$betas, n.kpc * n.class, n.samples))
  object$tau <- t(matrix(b$taus, n.kpc * n.class, n.samples)) 
  object$z <- t(matrix(b$zs, ntr * n.class, n.samples))
  object$sigmasq <- matrix(b$sigmasqs, n.samples, 1)

  object$n.class  <-  n.class
  object$n.kpc  <-  n.kpc
  object$n.iter  <-  n.iter
  object$thin <- thin
  object$intercept <- intercept
  object$rotate <- rotate 
  
  if (rotate)object$kPCA <- objectKPCA
  else object$kPCA <- NULL
  object$x <- x
  object$theta <- NULL
  
  class(object) <- c("bkpc.kern", "bkpc")   
  return(object)
}
