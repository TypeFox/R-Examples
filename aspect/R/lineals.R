`lineals` <-
function(data, level = "nominal", itmax = 100, eps = 1e-6)
{
  
  m <- dim(data)[2]                                      #number of variables
  lev.check <- match.arg(level, c("nominal","ordinal"), several.ok = TRUE)
  if (length(lev.check) < length(level)) stop("Level argument should be one of nominal or ordinal!")
  if (length(level) == 1) {
    level <- rep(level, m)
  } else {
    if (length(level) != m) stop("Length of level vector must correspond to number of variables!")
  }
  n <- dim(data)[1]                                      #number of observations
  r <- diag(m)
  t <- diag(m)
  fold <- Inf
  itel <- 1
  
  ncat <- sapply(1:m, function(j) length(table(data[,j])))  #number of categories for each variable
  ccat <- c(0,cumsum(ncat))
  
  y <- list()
  for (j in 1:m) y <- c(y,list(1:ncat[j]))              #list with category indices (inital score values)
  names(y) <- colnames(data)

  burt <- burtTable(data)                               #create Burt matrix (each variable is taken as categorical)
  nameslist <- apply(data, 2, unique)
  colnames(burt) <- rownames(burt) <- unlist(lapply(nameslist, sort))
  d <- diag(burt)                                       #diagonal of the Burt matrix

  for (j in 1:m) {                                      #compute category scores (normalized to y'burt y = 1)
    indj <- (ccat[j]+1):ccat[j+1]
    dj <- d[indj]
    y[[j]] <- y[[j]] - sum(dj*y[[j]])/n
    y[[j]] <- y[[j]]/sqrt(sum(dj*y[[j]]*y[[j]]))        #burt normalization
  }
  
  
#------------------ start lineals iterations ---------------------
  repeat {
    f <- 0                                              
    #------------- loss update -----------------
    for (j in 1:m) {
      indj <- (ccat[j]+1):ccat[j+1]
      yj <- y[[j]]
      for (l in 1:m) {
        indl <- (ccat[l]+1):ccat[l+1]
        dl <- d[indl]
        yl <- y[[l]]    
 	      r[j,l] <- sum(burt[indj,indl]*outer(yj,yl))    #correlation matrix FIXME!!! for ordinal
        c <- burt[indj,indl]%*%diag(1/pmax(1,dl))%*%burt[indl,indj]
        t[j,l] <- sum(c*outer(yj,yj))			             #correlation ratios
	      f <- f+(t[j,l]-r[j,l]^2)                       #loss update (cf. p. 448, de Leeuw 1988)
      }
    }

    #------------ scores update ----------------
    for (j in 1:m) {                                                 #score updating  
        indj <- (ccat[j]+1):ccat[j+1]
        nc <- ncat[j]
        yj <- y[[j]]
	      dj <- d[indj]
        c <- matrix(0,nc,nc)
	      for (l in 1:m) {                                 #j fixed, run over remaining variables
	       if (j != l) {
	        indl <- (ccat[l]+1):ccat[l+1]                #index of categories of variable j
          dl <- d[indl]                                #category frequencies
          yl <- y[[l]]
	        u <- burt[indj,indl]%*%(diag(1/pmax(1,dl)) - 2*outer(yl,yl))%*%burt[indl,indj] 
	        c <- c+u                                     #sum up u's 
	       }
	    }
	    c.norm <- c/sqrt(outer(dj,dj))                       #normalize 
	    e <- eigen(c.norm)
	    y[[j]] <- e$vectors[,nc]/sqrt(dj)                              #scores update
      
      #--------- pava --------------
      if (level[j] == "ordinal") {
        relfreq <- table(data[,j])/n
        pavares1 <- -pavasmacof(-y[[j]], relfreq)          #increasing          
        sspava1 <- sum((pavares1-y[[j]])^2)
        pavares2 <- -pavasmacof(y[[j]], relfreq)           #decreasing
        sspava2 <- sum((pavares2-y[[j]])^2)
        if (sspava1 <= sspava2) y[[j]] <- pavares1 else y[[j]] <- pavares2                  
      } 
      #---------- pava --------------  
     }
     if (((fold-f) < eps) || (itel == itmax)) break
     itel <- itel+1
     fold <- f
  }
  
#------------------------ end lineals iterations ---------------------
 

  if (itel == itmax) warning("Maximum iteration limit reached!")

  dummy.mat <- as.matrix(expandFrame(data))
  scorevec <- unlist(y)
  scoremat <- t(apply(dummy.mat, 1, function(xx) scorevec[which(xx == 1)]))
  colnames(scoremat) <- colnames(data)
  rownames(r) <- colnames(r) <- colnames(t) <- rownames(t) <- colnames(scoremat)
  for (i in 1:length(y)) {
    names(y[[i]]) <- sort(unique(data[,i]))                       #label categories
    y[[i]] <- -y[[i]]                                                           
    y[[i]] <- as.matrix(y[[i]])                                   #represent as matrix
    colnames(y[[i]]) <- "score"
  }
  
  r <- cor(scoremat)
  
  result <- list(loss = f, catscores = y, cormat = r, cor.rat = t, indmat = dummy.mat, 
  scoremat = scoremat, data = data, burtmat = burt, niter = itel, call = match.call())
  class(result) <- "aspect"
  result
}

