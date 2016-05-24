`corAspect` <-
function(data, aspect = "aspectSum", level = "nominal", itmax = 100, eps=1e-6, ...) 
{
# aspect ... function names: either "aspectSum", "aspectAbs", "aspectSMC", "aspectSumSMC", "aspectEigen", "aspectDeterminant"
# or a user-specified function (specific argument for this funtion go into "...")

  m <- dim(data)[2]
  lev.check <- match.arg(level, c("nominal","ordinal"), several.ok = TRUE)
  if (length(lev.check) < length(level)) stop("Level argument should be one of nominal or ordinal!")
  if (length(level) == 1) {
    level <- rep(level, m)
  } else {
    if (length(level) != m) stop("Length of level vector must correspond to number of variables!")
  }
  n <- dim(data)[1]
  r <- diag(m)
  fold <- -Inf
  itel<-1

  ncat <- sapply(1:m,function(j) length(table(data[,j])))
  ccat <- c(0,cumsum(ncat))                                   #cumulated number of categories (plus 0 as first element)
  y <- list()
  for (j in 1:m) y<-c(y,list(1:ncat[j]))
  names(y) <- colnames(data)

  burt <- burtTable(data)                                      #compute Burt matrix
  nameslist <- apply(data, 2, unique)
  colnames(burt) <- rownames(burt) <- unlist(lapply(nameslist, sort))
  d <- diag(burt)

  for (j in 1:m) {                                             #initial scaling of y such that y'*burt*y = 1
    indj <-(ccat[j]+1):ccat[j+1]
    dj <- d[indj]
    y[[j]] <- y[[j]]-sum(dj*y[[j]])/n                          #category scores (theta in paper)
    y[[j]] <- y[[j]]/sqrt(sum(dj*y[[j]]*y[[j]])/n)
  }

  #----------------- aspect string/function check --------------
  #functions must return the value (1st argument) and the derivative (second argument) 
  
  if (!is.function(aspect)) {                                        
    if (aspect == "aspectSum") aspectfun <- aspectSum                 #r, pow = 1 (power for r); sum of r_ij^pow 
    if (aspect == "aspectAbs") aspectfun <- aspectAbs                 #r, pow = 1 (power for r); sum of |r_ij|^pow
    if (aspect == "aspectSMC") aspectfun <- aspectSMC                 #r, targvar = 1 (index of target variable y); squared multiple correlation 
    if (aspect == "aspectSumSMC") aspectfun <- aspectSumSMC           #r; sum of SMC's between each SMC combination  
    if (aspect == "aspectEigen") aspectfun <- aspectEigen             #r, p (number of eigenvalues); sum of first p eigenvalues
    if (aspect == "aspectDeterminant") aspectfun <- aspectDeterminant #r; -log determinant of r
  } else {
    aspectfun <- aspect                                               #r needed plus additonal arguments passed by ...
  #FIXME: implement check for user-defined aspect (list output, 2 elements function value, derivative)
  }

  #---------------------- end aspect check ---------------------
  
  #------------------------ begin optimization -------------------
  repeat {

    #updates correlation matrix  
    for (j in 1:m) {                             
       indj <- (ccat[j]+1):ccat[j+1]
	     for (l in 1:m) {
	       indl <- (ccat[l]+1):ccat[l+1]
	       r[j,l] <- sum(y[[j]]*(burt[indj,indl]%*%y[[l]]))/n          #correlation matrix R(theta)
	    }
    }

    #apply aspect to correlation matrix
    #a <- aspectfun(r)
    a <- aspectfun(r, ...)                        #call aspect as a function of the correlation matrix r (and additional parameters)
    f <- a[[1]]                                   #value of the aspect function (only needed for convergence checking)
    g <- a[[2]]                                   #first derivative (needed for score update)

    #update scores
    for (j in 1:m) {                              #variable index
      indj <- (ccat[j]+1):ccat[j+1]               #row index
      y[[j]] <- rep(0,ncat[j])
      dj <- d[indj]                               #frequency vector from Burt matrix
      for (l in 1:m) {
	     indl <- (ccat[l]+1):ccat[l+1]             #subsetting indices from Burt matrix
	     if (j != l) y[[j]] <- y[[j]] + (g[j,l]*burt[indj,indl]%*%y[[l]])  #\sum dphi/dr_jl*burt*theta_l
      }
      y[[j]] <- y[[j]]/dj                         #normalize scores
      y[[j]] <- y[[j]]-sum(dj*y[[j]])/n
      y[[j]] <- y[[j]]/sqrt(sum(dj*y[[j]]*y[[j]])/n)
      
      #------- pava -----------
      if (level[j] == "ordinal") {
        relfreq <- table(data[,j])/n
        y[[j]] <- as.matrix(pavasmacof(y[[j]], relfreq))
      }
      #------- pava ----------- 
      
    }
    if (((f-fold) < eps) || (itel == itmax)) break
    itel<-itel+1
    fold<-f
  }
  #--------------------------- end optimization -------------------
  

  #returns: function value (total discrepancy f), category scores y, correlation matrix of the scores r, eigenvalues (of r) eigen
  if (itel == itmax) warning("Maximum iteration limit reached!")

  dummy.mat <- as.matrix(expandFrame(data))
  scorevec <- unlist(y)
  scoremat <- t(apply(dummy.mat, 1, function(xx) scorevec[which(xx == 1)]))
  colnames(scoremat) <- colnames(data)
  rownames(r) <- colnames(r) <- colnames(scoremat)
  for (i in 1:length(y)) {
    rownames(y[[i]]) <- sort(unique(data[,i]))             
    #rownames(y[[i]]) <- unique(data[,i])
    colnames(y[[i]]) <- "score"
  }
 
  r <- cor(scoremat)
  
  result <- list(loss = f, catscores = y, cormat = r, eigencor = eigen(r,only.values=TRUE)$values, indmat = dummy.mat, scoremat = scoremat, data = data, burtmat = burt, niter = itel, call = match.call())
  class(result) <- "aspect"
  result
}

