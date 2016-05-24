#-------------------------------------------------------------------------------
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
## METHOD 1
## the alternate mehod
#-------------------------------------------------------------------------------
# This function allows precision, neighbour and polygons information to be inputed
# The code was taken from Simon Wood code, which differs from Tomas Kneib function bnd2gra,
# which we use in the past to create the precision matrix.
# Our function now takes account the neighbour with a single point in common.
MRFA <- function(y, x,  
                 precision = NULL,
                 neighbour = NULL,
                     polys = NULL,
                      area = NULL,
                   weights = rep(1,length(y)),  
                    lambda = NULL,
                        df = NULL,
                     start = 10) # this this staring value for lambda i.e. start.lambda                    
{
  ##-----------------------------------------------------------------------------  
  ## local functions 
  #------------------------------------------------------------------------------ 
  # local 1 :the fitting 
  regpen <- function(y, weights, lambda)  
  {
    # fv <- solve(W+lambda*G,  weights*y)
    beta <- solve(XWX + lambda*G, XWy) 
    names(beta) <- levels(x)
    fv <- beta[x]              
    H <- solve(XWX + lambda*G, XWX)
    edf <- sum(diag(H))
    fit <- list(fv=fv, beta=beta, edf=edf, var=diag(H))
    fit 
  }
  ##-----------------------------------------------------------------------------
  ## local 2 : local function to get df using eigen values
  edf1_df <- function(lambda)
  {
    edf <-  sum(1/(1+lambda*UDU$values))
    (edf-df)
  }  
  
  ##-----------------------------------------------------------------------------
  ## local 3 : creating the precision matrix from the neighbour information
  ##          same as the external marix with the same name but here takes 
  ##          local information
  nb2prec <- function(neighbour)
  {
    a.name <- names(neighbour$nb)
    if (all.equal(sort(a.name), sort(levels(k))) != TRUE) 
      stop("mismatch between neighbour/polys supplied area names and data area names")
    np <- nfv
    G <- matrix(0, np, np)
    rownames(G) <- colnames(G) <- levels(k)
    for (i in 1:np) 
    {
      ind <- neighbour$nb[[i]]
      lind <- length(ind)
      G[a.name[i], a.name[i]] <- lind
      if (lind > 0) 
        for (j in 1:lind) G[a.name[i], a.name[ind[j]]] <- -1
    }
    if (sum(G != t(G)) > 0) 
      stop("Something wrong with auto- penalty construction")
    if (all.equal(sort(a.name),a.name)!=TRUE) 
    { ## re-order penalty to match X
      G <- G[levels(x),]
      G <- G[,levels(x)]
    }
    G
  }
  ## end of local functions 
  ##------------------------------------------------------------------------------
  ##------------------------------------------------------------------------------
  ##------------------------------------------------------------------------------
  scall <- deparse(sys.call(), width.cutoff = 500L)
  if (!is(x, "factor")) stop("x must be a factor")
  N <- length(y)
  if (any(is.na(weights))) weights[is.na(weights)] <- 0
  nobs <- sum(!weights==0)
  # if (N>nobs) warning("some observations are weighted out")
  k <- area
  if (is.null(k)) k <- factor(levels(x),levels=levels(x)) # default knots = all regions are in the data
  #the line above can change the order of the precision matrix
  else{
    if (class(area)=="character") k <- as.factor(k)
    if (!(class(k)=="character"||class(k)=="factor")) 
      stop("area must be a factor or a chacacter vector")
  }
  if (length(levels(x))>length(levels(k))) 
    stop("MRF basis dimension set too high")
  if (sum(!levels(x)%in%levels(k))) 
    stop("data contain regions that are not contained in the area specification")
  x <- factor(x,levels=levels(k))
  X <- model.matrix(~x-1,) 
  nfv <- nlevels(x)
  if (is.null(precision)&&is.null(neighbour)&&is.null(polys))
    stop("precision matrix, boundary polygons and/or neighbours list must be supplied")
  if (!is.null(precision))
  { 
    if (!is.matrix(precision)||dim(precision)[1]!=nfv||dim(precision)[2]!=nfv) 
      stop("the precision matrix is not suitable")
    G <- precision 
  } 
  # check the precision matrix
  if (!is.null(neighbour)&&is.null(precision))
  { # if neighbour exits then calculate the precision
    G   <- nb2prec(neighbour)
  }
  if (!is.null(polys)&&is.null(neighbour)&&is.null(precision))  
  { # if polys exits then calculate the precision
    a.name <- names(polys)
    d.name <- unique(a.name[duplicated(a.name)])
    if (length(d.name)) 
    {
      for (i in 1:length(d.name)) 
      {
        ind <- (1:length(a.name))[a.name == d.name[i]]
        for (j in 2:length(ind)) 
          polys[[ind[1]]] <- rbind(polys[[ind[1]]], c(NA, NA), polys[[ind[j]]])
      }
      #now delete the un-wanted duplicates
      ind <- (1:length(a.name))[duplicated(a.name)]
      if (length(ind) > 0) 
        for (i in length(ind):1) polys[[ind[i]]] <- NULL
    }#polygon list in correct format
    neighbour <- polys2nb(polys)
    G   <- nb2prec(neighbour)
  }
  XWX <- diag.spam(diag(t(X)%*%diag(weights)%*%X))
  #     XWX <- diag.spam(tapply(weights, x, sum)) 
  XWy <- as.vector(t(X)%*%diag(weights)%*%matrix(y, nrow=length(y), ncol=1))
  #   XWy <- tapply(weights * y, x, sum)
  G <- spam(as.numeric(G), nrow=nfv)
  sig2e <- sig2b <- 0 # initialiase the sigmas
  ##      if lambda is given then just fit the model 
  if (is.null(df)&&!is.null(lambda)||!is.null(df)&&!is.null(lambda))
  {
    fit <- regpen(y, weights=weights, lambda=lambda)
    # fv <- X %*% fit$beta   
  }
  else if (is.null(df)&&is.null(lambda)) 
  {
    lambda <- start
    for (it in 1:200) 
    {
      fit <- regpen(y, weights=weights, lambda=lambda)
      gamma. <- fit$beta
      fv <- fit$fv          
      sig2e <- sum(weights * (y - fv) ^ 2) / (nobs - fit$edf)
      sig2b <- sum(t(gamma.)%*%G%*%gamma.)/(fit$edf)
      lambda.old <- lambda
      lambda <- sig2e/sig2b
      if (lambda < 1e-07) lambda <- 1e-07
      if (lambda > 1e+07) lambda <- 1e+07
      if (abs(lambda - lambda.old) < 1e-07 || lambda > 1e+10) break
    }
  } 
  else # case 3 : if df are required---------------------------------
{
  QR <- qr(sqrt(weights)*X)
  Rinv <- solve(qr.R(QR))
  UDU <- eigen(t(Rinv)%*%G%*%Rinv)           
  lambda <- if (sign(edf1_df(0))==sign(edf1_df(100000))) 100000  # in case they have the some sign
  else  uniroot(edf1_df, c(0,100000))$root
  # if (any(class(lambda)%in%"try-error")) {lambda<-100000}   
  fit <- regpen(y, weights, lambda)
  #  fv <- X %*% fit$beta          
}
  
  # fit <- regpen(y, weights=weights, lambda=lambda) 
  sig2e <- sum(weights * (y - fit$fv) ^ 2) / (nobs - fit$edf)
  sig2b <- sum(t(fit$beta)%*%G%*%fit$beta)/(fit$edf)   
  par <- c(log(sig2e), log(sig2b))
  names(par) <- c("log(sige^2)", "log(sigb^2)")
  # saving stuff 
  fit <- list(fitted = fit$fv, 
              df = fit$edf, 
              lambda = lambda, 
              sig2e =  sig2e, 
              sige = sqrt(sig2e),
              sig2b = sig2b,
              sigb = sqrt(sig2b),
              par = list(par=par, se=c(NA_real_, NA_real_)), # rubish ???
              y = y,
              x = x, 
              # sumW = (sum(weights*(y-fit$fv)^2)),
              nobs = nobs,  
              #gaGga = sum(t(fit$beta)%*%G%*%fit$beta), 
              beta = fit$beta, # test
              # XWX = XWX,
              # XWy = XWy,
              G = G, 
              var = fit$var[x],
              method = "altenating",
              value.of.Q = NULL,
              deviance.Q = 0,
              weights = weights,
              N = N,
              call = scall,
              rss = sum(weights*(y-fit$fv)^2),
              aic = sum(weights*(-2*dNO(y, mu=fit$fv, sigma=sqrt(sig2e), log=TRUE)))+2*(fit$edf+1) , 
              sbc = sum(weights*(-2*dNO(y, mu=fit$fv, sigma=sqrt(sig2e), log=TRUE)))+log(nobs)*(fit$edf+1),
              deviance = sum(weights*(-2*dNO(y, mu=fit$fv, sigma=sqrt(sig2e), log=TRUE))))
  class(fit) <- "MRF"
  return(fit)     
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------