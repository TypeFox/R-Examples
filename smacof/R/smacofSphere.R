smacofSphere <- function(delta, algorithm = c("dual", "primal"), ndim = 2, type = c("ratio", "interval", "ordinal"),
                         weightmat = NULL, init = "torgerson",  ties = "primary", verbose = FALSE,
                         penalty = 100, relax = FALSE, modulus = 1, itmax = 1000, eps = 1e-6)
{
  # penalty ... penalty term kappa >0, 100 is reasonable
  
  
 alg <- match.arg(algorithm, c("dual", "primal"))
 type <- match.arg(type, c("ratio", "interval", "ordinal"))
 
 if (alg == "dual") {
  diss <- delta
  if ((is.matrix(diss)) || (is.data.frame(diss))) diss <- strucprep(diss)  #if data are provided as dissimilarity matrix
  checkdiss(diss)
  
  p <- ndim
  n <- attr(diss,"Size")
  if (p > (n - 1)) stop("Maximum number of dimensions is n-1!")
  
  nn <- n*(n-1)/2
  m <- length(diss)
  
  if (is.null(attr(diss, "Labels"))) attr(diss, "Labels") <- paste(1:n)
  
  if (is.null(weightmat)) {
    wgths <- initWeights(diss)
  }  else  wgths <- as.dist(weightmat)
  
  dhat <- normDissN(diss,wgths,1)            #normalize dissimilarities
  dhat[is.na(dhat)] <- 1     ## in case of missing dissimilarities, pseudo value for dhat 
    
  
  ## --- starting values
  x <- initConf(init, diss, n, p)
  xstart <- x
  
  if (relax) relax <- 2 else relax <- 1 
  
  mn <- c(1,rep(0,n))
  diss <- as.dist(rbind(0,cbind(0,as.matrix(diss))))         #add row/column
  
  wgths1 <- as.dist(rbind(0,cbind(0,as.matrix(wgths))))      #distances weights
  wgths2 <- as.dist(outer(mn,mn,function(x,y) abs(x-y)))    
  dhat1 <- as.dist(rbind(0,cbind(0,as.matrix(dhat))))        #0's in the first column
  dhat2 <- mean(sqrt(rowSums(x^2)))*wgths2
  
  x <- rbind(0,x)
  w <- vmat(wgths1+penalty*wgths2); v<-myGenInv(w); itel<-1;
  d <- dist(x)
  lb <- sum(wgths1*d*dhat1)/sum(wgths1*d^2)
  x <- lb*x
  d <- lb*d
  sold1 <- sum(wgths1*(dhat1-d)^2)
  sold2 <- sum(wgths2*(dhat2-d)^2)
  sold <- sold1+penalty*sold2
  
  #---------------- begin majorization ---------------- 
  repeat
  {
    b <- bmat(dhat1,wgths1,d)+penalty*bmat(dhat2,wgths2,d)
    y <- v %*% (b %*% x)
    y <- x+relax*(y-x)
    e <- dist(y)
    ssma1 <- sum(wgths1*(dhat1-e)^2)                       #stress for delta
    ssma2 <- sum(wgths2*(dhat2-e)^2)                       #penalty term 
    ssma <- ssma1+penalty*ssma2                            #joint stress value
    
    #---- nonmetric MDS --------
    if ((type == "ordinal") && ((itel%%modulus) == 0)) {
      if (ties=="primary") daux<-monregP(diss,e,wgths1)
      if (ties=="secondary") daux<-monregS(diss,e,wgths1)
      if (ties=="tertiary") daux<-monregT(diss,e,wgths1)
      daux<-vecAsDist(daux); dhat1<-normDissN(daux,wgths1,1)
    }
    
    ## --- interval MDS
    if (type == "interval") {
      Amat <- cbind(1, as.vector(diss), as.vector(diss)^2) 
      daux <- nnlsPred(Amat, as.vector(e), as.vector(wgths1))$pred
      daux <- vecAsDist(daux)
      dhat1 <- normDissN(daux,wgths1,1)
    }
    
    dhat2 <- mean(e[1:n])*wgths2
    snon1 <- sum(wgths1*(dhat1-e)^2)
    snon2 <- sum(wgths2*(dhat2-e)^2)                       
    snon <- snon1+penalty*snon2                            #nonmetric joint stress
    
    if (verbose) cat("Iteration: ",formatC(itel,width=3, format="d")," Stress (not normalized): ",
                     formatC(c(snon),digits=8,width=12,format="f"),"\n")
    
    
    if (((sold-snon)<eps) || (itel == itmax)) break()      #convergence
    
    x <- y                               #updates
    d <- e
    sold <- snon
    itel <- itel+1
    
    if (itel == itmax) warning("Iteration limit reached! Increase itmax argument!")
    
  }
  #-------------- end majorization ---------------
  
  colnames(y) <- paste("D",1:(dim(y)[2]),sep="")
  rownames(y) <- labels(diss)
  attr(dhat1, "Labels") <- labels(diss)
  attr(e, "Labels") <- labels(diss)
  
  stress <- sqrt(snon/nn)                   #stress normalization
  
  e.temp <- as.dist(as.matrix(e)[,-1][-1,])      #remove dummy vector
  dummyvec <- as.matrix(e)[,1][-1]
  confdiss <- normDissN(e.temp, wgths, 1)        #final normalization to n(n-1)/2
  
  # point stress 
  dhat <- as.dist(as.matrix(dhat1)[-1,-1])
  spoint <- spp(dhat, confdiss, wgths)
  rss <- sum(spoint$resmat[lower.tri(spoint$resmat)])  ## residual sum-of-squares
  
  ss <- y[1,]
  y <- t(apply(y, 1, function(xx) xx - ss)[,-1] )  ## correct the configurations for plotting
  
  
## -------------------------------------- primal algorithm ------------------------------------------  
} else {
  diss <- delta
  if ((is.matrix(diss)) || (is.data.frame(diss))) diss <- strucprep(diss)  #if data are provided as dissimilarity matrix
  checkdiss(diss)
  
  p <- ndim
  n <- attr(diss,"Size")
  if (p > (n - 1)) stop("Maximum number of dimensions is n-1!")
  
  nn <- n*(n-1)/2
  m <- length(diss)
  if (is.null(attr(diss, "Labels"))) attr(diss, "Labels") <- paste(1:n)
  
  
  if (is.null(weightmat)) {
    wgths <- initWeights(diss)
  }  else  wgths <- as.dist(weightmat)
  
  dhat <- normDissN(diss,wgths,1)            #normalize dissimilarities
  dhat[is.na(dhat)] <- 1     ## in case of missing dissimilarities, pseudo value for dhat 
  
 
  ## --- starting values
  x <- initConf(init, diss, n, p)
  xstart <- x
  
  w <- vmat(wgths)
  v <- myGenInv(w)
  itel<-1;
  
  x <- x/sqrt(rowSums(x^2))
  #FIXME!!!
  d <- dist(x)                           #distance computation (to be extended with geodesic)
  
  lb <- sum(wgths*d*dhat)/sum(wgths*d^2)
  x <- lb*x
  d <- lb*d
  sold <- sum(wgths*(dhat-d)^2)          #initial stress
  
  #------------------------- begin majorization ---------------------------------
  repeat {
    b <- bmat(dhat,wgths,d)
    y <- v%*%b%*%x                       #Guttman transform
    y <- sphereProj(y,w)                 #projection on the sphere
    
    e <- dist(y)                         #extension: distances for Y to be enhanced with geodesics)
    ssma <- sum(wgths*(dhat-e)^2)        #metric stress
    
    if (type == "ordinal") {                       #nonmetric versions
      if ((itel%%modulus) == 0) {
        if (ties=="primary") daux <- monregP(diss,e,wgths)
        if (ties=="secondary") daux <- monregS(diss,e,wgths)
        if (ties=="tertiary") daux <- monregT(diss,e,wgths)
        daux <- vecAsDist(daux)
        dhat <- normDissN(daux,wgths,1)
      }
    }
    
    ## --- interval MDS
    if (type == "interval") {
      Amat <- cbind(1, as.vector(diss), as.vector(diss)^2) 
      daux <- nnlsPred(Amat, as.vector(e), as.vector(wgths))$pred
      daux <- vecAsDist(daux)
      dhat <- normDissN(daux,wgths,1)
    }
    
    
    snon <- sum(wgths*(dhat-e)^2)        #nonmetric stress
    if (verbose) cat("Iteration: ",formatC(itel,width=3, format="d")," Stress (not normalized): ",
                     formatC(c(snon),digits=8,width=12,format="f"),"\n")
    if (((sold-snon)<eps) || (itel == itmax)) break()
    
    x <- y                               #updates
    d <- e
    sold <- snon
    itel <- itel+1
  }
  #----------------------------- end majorization -------------------------------
  
  colnames(y) <- paste("D",1:(dim(y)[2]),sep="")
  rownames(y) <- labels(diss)
  attr(dhat, "Labels") <- labels(diss)
  attr(e, "Labels") <- labels(diss)
  dhat[is.na(diss)] <- NA                 ## in case of NA's
  
  stress <- sqrt(snon/nn)                   #stress normalization
  
  confdiss <- normDissN(e, wgths, 1)        #final normalization to n(n-1)/2
  
  # point stress 
  spoint <- spp(dhat, confdiss, wgths)
  rss <- sum(spoint$resmat[lower.tri(spoint$resmat)])  ## residual sum-of-squares
  
  dummyvec <- NA
  if (itel == itmax) warning("Iteration limit reached! Increase itmax argument!")
}
  
  
result <- list(delta = diss, dhat = dhat, confdiss = confdiss, conf = y, 
               stress = stress, spp = spoint$spp, ndim = p, weightmat = wgths, resmat = spoint$resmat, rss = rss, 
               dummyvec = dummyvec, init = xstart, 
               model = "Spherical SMACOF", niter = itel, nobj = n, type = type, 
               algorithm = alg, call = match.call())

class(result) <- c("smacofSP", "smacof")
  result
}
