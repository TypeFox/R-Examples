## smacof with external constraints on the configuration (de Leeuw & Heiser, 1980; Borg & Groenen, p. 236)

smacofConstraint <- function(delta, constraint = "linear", external, ndim = 2, type = c("ratio", "interval", "ordinal", "mspline"), 
                             weightmat = NULL, init = NULL, ties = "primary", verbose = FALSE, 
                             modulus = 1, itmax = 1000, eps = 1e-6,  
                             spline.intKnots = 4, spline.degree = 2, 
                             constraint.type = c("ratio", "interval", "ordinal", "spline", "mspline"), 
                             constraint.ties = "primary", 
                             constraint.spline.intKnots = 2, 
                             constraint.spline.degree = 2)
{
  # diss ... dissimilarity matrix
  # constraint ... either "linear", "unique", "diagonal", or a user-specified function
  # external ... external data for X-decomposition (Z in paper), or list with "simplex", or "circumplex"
  # weightmat ... weight structure. if not specified, weights is 1-structure
  # ndim ... number of dimensions
  # init ... starting configuration
  # metric ... if TRUE, metric MDS, if FALSE, non-metric
  # ties ... ties for pava (primary, secondary, tertiary)
  # modulus ... modulus for nonmetric update
  # itmax ... maximum number of iterations
  # eps ... change in loss function
  # spline.intKnots ... no of interior knots for spline
  # spline.degree ... degree of spline 
  # constraint.type ... transformation of external variables, either "ratio", "interval", "ordinal", "spline", or "mspline" 
  # constraint.ties ... treatment of ties for ordinal transformation of external variables (primary, secondary, tertiary) 
  # constraint.spline.intKnots ... number of interior knots of spline transformation of external variables  
  # constraint.spline.degree ... degree of spline transformation of external variables 
  
  type <- match.arg(type, c("ratio", "interval", "ordinal", "mspline"), several.ok = FALSE) 
  ties <- match.arg(ties, c("primary", "secondary", "tertiary"), several.ok = FALSE) 
  constraint.type <- match.arg(constraint.type, c("ratio", "interval", "ordinal", "spline", "mspline"), several.ok = FALSE) 
  constraint.ties <- match.arg(constraint.ties, c("primary", "secondary", "tertiary"), several.ok = FALSE) 
  
  diss <- delta
  if ((is.matrix(diss)) || (is.data.frame(diss))) {
    diss <- strucprep(diss)  #if data are provided as dissimilarity matrix
    attr(diss, "Labels") <- rownames(delta)
  }
  checkdiss(diss)
  
  p <- ndim                                     
  n <- attr(diss,"Size")
  if (p > (n - 1)) stop("Maximum number of dimensions is n-1!")
  
  nn <- n*(n-1)/2
  m <- length(diss)
  
  
  ## --- starting values 
  startconf <- init
  if (!is.null(startconf)) startconf <- as.matrix(init)   # x as matrix with starting values   
  xstart <- startconf
  
  if (is.null(attr(diss, "Labels"))) attr(diss, "Labels") <- paste(1:n)
  
  ## sanity check external
  if (is.data.frame(external)) external <- as.matrix(external)
  if (!is.list(external)) {
    if (ncol(external) < p) stop("Number of external variables can not be smaller than the number of MDS dimensions!")
  }
  #---- external specification -----
  if (is.matrix(external)) {
    extvars <- list()
    for (s in 1:ncol(external)){
      # Prepare for optimal scaling of transformations
      constraint.trans <- constraint.type
      if (constraint.trans == "ratio"){
        constraint.trans <- "none"
      } else if (constraint.trans=="ordinal" & ties=="primary"){
        constraint.trans <- "ordinalp"
      } else if(constraint.trans=="ordinal" & ties=="secondary"){
        constraint.trans <- "ordinals"
      } else if(constraint.trans=="ordinal" & ties=="tertiary"){
        constraint.trans <- "ordinalt"
      } else if(constraint.trans=="spline"){
        constraint.trans <- "spline"
      } else if(constraint.trans=="mspline"){
        constraint.trans <- "mspline"
      }
      
      extvars[[s]] <- transPrep(external[,s]-mean(external[,s], na.rm = TRUE),
                                trans = constraint.trans, 
                                spline.intKnots = constraint.spline.intKnots, 
                                spline.degree = constraint.spline.degree,
                                missing = "multiple")
      
      external[,s] <- extvars[[s]]$xInit - mean(extvars[[s]]$xInit)
    }
  } 
  
  simpcirc <- FALSE
  if (is.list(external)) {                     
    if (external[[1]] == "simplex") {                           #simplex specification
      d2 <- external[[2]]
      if (d2 >= n) stop("Simplex dimension must be < n!")
      external <- diag(1, n)[,1:d2]
      external[lower.tri(external)] <- 1
    }
    if (external[[1]] == "circumplex") {                        #circumplex specification
      d2 <- external[[2]]
      if (d2 >= n) stop("Circumplex dimension must be <= n!")
      k1 <- external[[3]]
      k2 <- external[[4]]
      if (k2 <= k1) stop("k2 must be > k1")
      external <- matrix(0, nrow = n, ncol = d2)
      ind.all <- expand.grid(1:n, 1:d2)
      inddiff <- apply(ind.all, 1, function(xx) abs(diff(xx)))      
      ind.good <- which((inddiff >= k1) + (inddiff <= k2) == 2) #k1 <= |i-s| <= k2
      el1 <- as.matrix(ind.all[ind.good,])
      external[el1] <- 1
    }
    simpcirc <- TRUE
  }
  
  K <- dim(external)[2]  
  #-------- end external -----------
  
  if (is.null(weightmat)) {
    wgths <- initWeights(diss)
  }  else  wgths <- as.dist(weightmat)
  
  # Prepare for optimal scaling
  trans <- type
  if (trans=="ratio"){
    trans <- "none"
  } else if (trans=="ordinal" & ties=="primary"){
    trans <- "ordinalp"
  } else if(trans=="ordinal" & ties=="secondary"){
    trans <- "ordinals"
  } else if(trans=="ordinal" & ties=="tertiary"){
    trans <- "ordinalt"
  } else if(trans=="spline"){
    trans <- "mspline"
  }
  disobj <- transPrep(diss,trans = trans, spline.intKnots = spline.intKnots, spline.degree = spline.degree)  
  
  
  dhat <- normDissN(diss,wgths,1)               #normalize dissimilarities
  dhat[is.na(dhat)] <- 1     ## in case of missing dissimilarities, pseudo value for dhat 
  
  w <- vmat(wgths)                              #matrix V
  v <- myGenInv(w)                              #Moore-Penrose inverse
  itel <- 1
  
  
  #----------- pre-specified functions for constraints -----------
  # linear constraint (de Leeuw & Heiser, 1980, p.515), X=ZC 
  if (!is.function(constraint)) { 
    if (constraint == "linear") {
      constrfun <-function(x,w,external) {
        external%*%solve(crossprod(external,w%*%external),crossprod(external,w%*%x))
      }
      if (is.null(xstart)) xstart <- matrix(rnorm(n*p),n,p)                     #starting value for X (before constraints)
    }
    
    # C restricted to be diagonal   
    if (constraint == "diagonal") {
      constrfun <- function(x,w,external) {
        return(external%*%diag(colSums(external*(w%*%x))/colSums(external*(w%*%external))))  
      }
      if (is.null(xstart)) xstart <- matrix(rnorm(n*K),n,K)                     #starting value for X (before constraints) of dimension n x K
    }
    
    
    # X with uniqueness coordinates  
    if (constraint == "unique") {
      constrfun <- function(x,w,external) {
        n <- dim(x)[1]
        p <- dim(x)[2]-n 
        return(cbind(x[,1:p],diag(diag(w%*%x[,p+(1:n)])/diag(w))))
      }
      if (is.null(xstart)) xstart <- cbind(matrix(rnorm(n*p),n,p), diag(1,n))   #starting value for X (before constraints) including diagonal matrix
    }
  } else {   # user-specified
    constrfun <- constraint
    if (is.null(xstart)) stop("Starting configuration must be specified!")
  }
  
  
  #---------- end pre-specified functions for constraints -------
  
  #x <- constrfun(xstart,w,external)                    #compute X
  if (constraint %in% c("linear","diagonal") & !simpcirc){
    # First make random weight matrices
    ncol.ext <- ncol(external)
    if (constraint == "linear"){
      # Make weight matrix C from a rotation matrix out of the left singular vectors of
      # a centering matrix
      C <- svd(diag(ncol.ext) - 1/ncol.ext)$u[, 1:ndim]
      #C <- matrix(runif(ncol.ext * ndim), ncol.ext, ndim)
    } else if (constraint == "diagonal") {
      C <- diag(ncol.ext)  # Make an initial C = I
    }
    # Initialize the optimally scaled external variables
    x.unc <- xstart
    x.con <- matrix(0, n, p)
    for (s in 1:ncol.ext){  # Find initial constrained configuration
      target <- x.unc %*% C[s, ]/sum(C[s, ]^2)
      loss <- sum((x.unc - outer(external[, s], C[s, ]) )^2)
      loss.old <- loss + 2 * eps
      while (loss.old - loss > eps) {  
        loss.old <- loss
        tt.plus <- transform(target, extvars[[s]], normq = 0)     # Compute update for external variable s
        tt.min  <- transform(-target, extvars[[s]], normq = 0)    # Compute update for external variable s
        if (sum((tt.plus$res - target)^2) < sum(((tt.min$res + target))^2) ) {
          external[, s] <- tt.plus$res
        } else {
          external[, s] <- tt.min$res
        }
        #x.unc <- x.unc - outer(external[, s],C[s, ])
        if (constraint == "linear"){
          C[s, ] <- t(external[,s, drop = FALSE]) %*% x.unc / sum((external[,s])^2)
          target <- x.unc %*% C[s, ]/sum(C[s, ]^2)         
        } else if (constraint == "diagonal") {
          C[s, s] <- t(external[, s, drop = FALSE]) %*% x.unc[, s, drop = FALSE] / sum((external[,s])^2)
          target <- x.unc[, s] / C[s, s]          
        } 
        loss <- sum((x.unc - outer(external[, s], C[s, ]) )^2)
      }
      x.unc <- x.unc - outer(external[, s], C[s, ])
    }  
    
    # Set external to column sum of squares n
    #external <- apply(external, 2, FUN = function(x){x <- x*(length(x)/sum(x^2))^.5})
    # Do an extra round of updates to get the column length restriction fine.
    updext.result <- updext(xstart,w,external,extvars,constraint)
    external <- updext.result$external
  }
  x <- constrfun(xstart,w,external)  
  
  d <- dist(x)                                         #distances X
  lb <- sum(wgths*d*dhat)/sum(wgths*d^2)               #denominator: normalization tr(X'VX) 
  x <- lb*x                                            #modify x with lb-factor
  d <- lb*d                                            #modify d with lb-factor
  sold <- sum(wgths*(dhat-d)^2)/nn                     #initial stress
  
  #------- begin majorization -----------
  repeat {                                         #majorization iterations
    b <- bmat(dhat,wgths,d)                        # B matrix
    y <- v %*% (b %*% x)                           # Y computation
    if (constraint %in% c("linear","diagonal") & !simpcirc){   # Update transformation of external variables
      updext.result <- updext(x,w,external,extvars,constraint)
      external <- updext.result$external
    }
    y <- constrfun(y,w,external)                   # update Y with corresponding constraints            
    e <- dist(y)                                   # Y distances
    ssma <- sum(wgths*(dhat-e)^2)                  # new stress
    
    dhat2 <- transform(e, disobj, w = wgths, normq = nn) # Transformation of the dissimilarities
    dhat <- dhat2$res
    
    
    snon <- sum(wgths*(dhat-e)^2)/nn               # nonmetric stress
    
    if (verbose) cat("Iteration: ",formatC(itel,width=3, format="d"),
                     " Stress (raw): ", formatC(c(snon),digits=8,width=10,format="f"),
                     " Difference: ", formatC(c(sold-snon),digits=8,width=10,format="f"),
                     "\n")
    
    if (((sold-snon)<eps) || (itel == itmax)) break()    #convergence 
    
    x <- y                                         #updates
    d <- e
    sold <- snon
    itel <- itel+1	
  }
  #------- end majorization -----------
  
  stress <- sqrt(snon)                   #stress normalization
  
  if (any(is.na(y))) {              #reduce ndim for external == simplex
    csy <- colSums(y)
    ind <- which(is.na(csy))
    y <- y[,-ind]
  }
  
  colnames(y) <- paste("D",1:(dim(y)[2]),sep="")
  rownames(y) <- labels(diss)
  dhat <- structure(dhat, Size = n, call = quote(as.dist.default(m=b)), class = "dist", Diag = FALSE, Upper = FALSE) 
  attr(dhat, "Labels") <- labels(diss)
  attr(e, "Labels") <- labels(diss)
  dhat[is.na(diss)] <- NA                 ## in case of NA's
  
  confdiss <- normDissN(e, wgths, 1)        #final normalization to n(n-1)/2
  
  ## stress-per-point 
  spoint <- spp(dhat, confdiss, wgths)
  rss <- sum(spoint$resmat[lower.tri(spoint$resmat)])  ## residual sum-of-squares
  
  if ((constraint == "diagonal") && (!simpcirc)) {
    if (p != ncol(y)) {
      warning("Number of dimensions is set equal to number of external variables!")
      p <- ncol(y)     ## adapt dimensions for diagonal restriction
    }
  }
  if (simpcirc) {
    if (p != ncol(y)) {
      warning("Number of dimensions is set equal to dimension defined by simplex/circumplex")
      p <- ncol(y)     ## adapt dimensions for simplex/circumplex fitting
    }
  }
  
  if (itel == itmax) warning("Iteration limit reached! Increase itmax argument!")
  
  ## compute C
  Z <- as.matrix(external)
  X <- y
  C <- ginv(t(Z)%*%Z)%*%t(Z)%*%X 
  if (constraint %in% c("linear","diagonal") && !simpcirc){
    for (s in 1:ncol(external)){
      extvars[[s]]$iord.prim <- updext.result$iord.prim[[s]]
      extvars[[s]]$final     <- external[,s]
      extvars[[s]]$c         <- C[s,]
    }
  } else {
    extvars = NULL
  }
  
  result <- list(delta = diss, dhat = dhat, confdiss = confdiss, conf = y, C = C, 
                 stress = stress, spp = spoint$spp, ndim = p, iord = dhat2$iord.prim, extvars = extvars,
                 external = external, weightmat = wgths, resmat = spoint$resmat, rss = rss, init = xstart, model = "SMACOF constraint", 
                 niter = itel, nobj = n, type = type, call = match.call()) 
  class(result) <- c("smacofB","smacof")
  result 
}



