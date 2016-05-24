fast.geomorph.compare.evol.rates <-
function (phy, A, gp, method="ML",ShowPlot = TRUE, iter = 999,censored=FALSE,force.diag=FALSE) 
{
  if(method=="REML") REML <- 1 else REML <- 0
  phy <- reorder(multi2di(phy),"postorder")
  if (length(dim(A)) == 3) {
    if (is.null(dimnames(A)[[3]])) {
      stop("Data matrix does not include taxa names as dimnames for 3rd dimension.")
    }
    x <- two.d.array(A)
  }
  if (length(dim(A)) == 2) {
    if (is.null(rownames(A))) {
      stop("Data matrix does not include taxa names as dimnames for rows.")
    }
    x <- A
  }
  if (is.vector(A) == TRUE) {
    if (is.null(names(A))) {
      stop("Data vector does not include taxa names as names.")
    }
    x <- as.matrix(A)
  }
  if (any(is.na(A)) == T) {
    stop("Data matrix contains missing values. Estimate these first (see 'estimate.missing').")
  }
  if (!is.factor(gp)) {
    stop("gp is not a factor.")
  }
  if (is.null(names(gp))) {
    stop("Factor contains no names. Use names() to assign specimen names to group factor.")
  }
  if (class(phy) != "phylo") 
    stop("tree must be of class 'phylo.'")
  ntaxa <- length(phy$tip.label)
  N <- nrow(x)
  if (N != dim(x)[1]) {
    stop("Number of taxa in data matrix and tree are not not equal.")
  }
  if (length(match(rownames(x), phy$tip.label)) != N) 
    stop("Data matrix missing some taxa present on the tree.")
  if (length(match(phy$tip.label, rownames(x))) != N) 
    stop("Tree missing some taxa in the data matrix.")
  p <- ncol(x)
  ones <- matrix(1, N)
  x <- x[phy$tip.label,,drop=FALSE]
  pY <- prep_multipic(x,phy = phy)
  pY_results <- do.call(multipic,pY)
  
  if(nlevels(gp==1)) return(list(sigma.d=mean(diag(crossprod(pY_results$contrasts)))/(length(phy$tip.label)-REML)))
  if(!censored) D.mat <- fast_transform(vcv(phy))
    
  if(censored)
  {
    sub_trees <- rep(phy,nlevels(gp))
    sub_ind <- vector("list",nlevels(gp))
    sub_Ns <- integer(nlevels(gp))
    sub_pY <- sub_pY_results <- vector("list",nlevels(gp))
    for(i in 1:nlevels(gp))
    {
      sub_trees[[i]] <- drop.tip(phy,name.check(phy = phy,data.names = names(gp)[gp==(levels(gp)[i])])$tree_not_data)
      sub_ind[[i]] <- match(sub_trees[[i]]$tip.label,phy$tip.label)
      sub_Ns[i] <- length(sub_trees[[i]]$tip.label)
      sub_pY[[i]] <- prep_multipic(x[sub_ind[[i]],,drop=FALSE],phy = sub_trees[[i]])
    }
  }
  sigmas <- numeric(nlevels(gp))
  
  sigma.d <- function(phy, x, N, gp) {
    gp <- gp[rownames(x)]
    ngps <- nlevels(gp)
    gpsz <- table(gp)
    pY$phe[1:N,] <- x
    pY_results <- do.call(multipic,pY)
    a.obs <- ones %*% pY_results$root
    if(censored & ngps>1)
    {
      for(i in 1:ngps)
      {
        sub_pY[[i]]$phe[1:sub_Ns[i],] <- x[sub_ind[[i]],,drop=FALSE]
        sub_pY_results[[i]] <- do.call(multipic,sub_pY[[i]])
        sigmas[i] <- mean(diag(crossprod(sub_pY_results[[i]]$contrasts)))/(sub_Ns[i]-REML)
      }
      sigma.d <- sigmas
      sigma.d.all <- mean(diag(crossprod(pY_results$contrasts)))/(N-REML)
    } else
    {
      dist.adj <- as.matrix(dist(rbind((D.mat %*% (x - a.obs)), 0)))
      vec.d2 <- dist.adj[N + 1, 1:N]^2
      sigma.d <- tapply(vec.d2, gp, sum)/(gpsz-REML)/p
      sigma.d.all <- sum(vec.d2)/(N-REML)/p
    }
    if (ngps == 1) {
      return(sigma.d.all)
    }
    if (ngps == 2) {
      sigma.d.rat <- max(sigma.d)/min(sigma.d)
      return(list(sigma.all = sigma.d.all, ratio = sigma.d.rat, 
                  sigma.d.all = sigma.d))
    }
    if (ngps > 2) {
      sigma.d.rat.gp <- array(0, dim = c(ngps, ngps))
      for (i in 1:(ngps - 1)) {
        for (j in 2:ngps) {
          tmp <- c(sigma.d[i], sigma.d[j])
          sigma.d.rat.gp[i, j] <- max(tmp)/min(tmp)
          diag(sigma.d.rat.gp) <- 0
          sigma.d.rat.gp[lower.tri(sigma.d.rat.gp)] <- 0
        }
      }
      sigma.d.rat <- max(sigma.d.rat.gp)
      return(list(sigma.all = sigma.d.all, ratio = sigma.d.rat, 
                  sigma.d.all = sigma.d, sigma.gp = sigma.d.rat.gp))
    }
  }
  if (nlevels(gp) == 1) {
    print("Single group. Sigma calculated for all specimens together.")
    sigmad.obs <- sigma.d(phy, x, ntaxa, gp)
    return(list(sigma.d = sigmad.obs))
  }
  if (nlevels(gp) > 1) {
    sigmad.obs <- sigma.d(phy, x, ntaxa, gp)
    a.obs <- ones %*% pY_results$root
    rate.mat <- crossprod(pY_results$contrasts)/(N-REML)
    if(force.diag & p>1) rate.mat <- diag(rep(mean(diag(rate.mat)),p))
    x.sim <- sim.char(phy, rate.mat, nsim = iter)
    sig.sim <- 1
    if (nlevels(gp) > 2) {
      gp.sig.sim <- array(1, dim = c(dim(sigmad.obs$sigma.gp)[1], 
                                     dim(sigmad.obs$sigma.gp)[1]))
    }
    rate.val <- rep(0, iter)
    sigmad.sim <- apply(x.sim,3,sigma.d,phy=phy,N=ntaxa,gp=gp)
    for (ii in 1:iter) {
      sig.sim <- ifelse(sigmad.sim[[ii]]$ratio >= sigmad.obs$ratio, 
                        sig.sim + 1, sig.sim)
      rate.val[ii] <- sigmad.sim[[ii]]$ratio
      if (nlevels(gp) > 2) {
        gp.sig.sim <- ifelse(sigmad.sim[[ii]]$sigma.gp >= sigmad.obs$sigma.gp, 
                             gp.sig.sim + 1, gp.sig.sim)
      }
    }
    sig.sim <- sig.sim/(iter + 1)
    if (nlevels(gp) > 2) {
      gp.sig.sim <- gp.sig.sim/(iter + 1)
      rownames(gp.sig.sim) <- colnames(gp.sig.sim) <- levels(gp)
      gp.sig.sim[lower.tri(gp.sig.sim)] <- NA
    }
    rate.val[iter + 1] = sigmad.obs$ratio
    if (ShowPlot == TRUE) {
      hist(rate.val, 30, freq = TRUE, col = "gray", xlab = "SigmaD ratio")
      arrows(sigmad.obs$ratio, 50, sigmad.obs$ratio, 5, 
             length = 0.1, lwd = 2)
    }
    if (nlevels(gp) > 2) {
      return(list(sigma.d = sigmad.obs$sigma.all, sigmad.all = sigmad.obs$sigma.d.all, 
                sigmad.ratio = sigmad.obs$ratio, pvalue = sig.sim, 
                pairwise.pvalue = gp.sig.sim))
    }
    else if (nlevels(gp) == 2) {
      return(list(sigma.d = sigmad.obs$sigma.all, sigmad.all = sigmad.obs$sigma.d.all, 
                  sigmad.ratio = sigmad.obs$ratio, pvalue = sig.sim))
    }
  }
}

fast.geomorph.physignal <-
function (phy, A, iter = 999, ShowPlot = TRUE, method = c("Kmult", 
                                                          "SSC")) 
{
  method <- match.arg(method)
  phy <- reorder(multi2di(phy,random=FALSE),"postorder")
  
  if (any(is.na(A)) == T) {
    stop("Data matrix contains missing values. Estimate these first (see 'estimate.missing').")
  }
  if (length(dim(A)) == 3) {
    if (is.null(dimnames(A)[[3]])) {
      stop("Data array does not include taxa names as dimnames for 3rd dimension.")
    }
    x <- two.d.array(A)
  }
  if (length(dim(A)) == 2) {
    if (is.null(rownames(A))) {
      stop("Data matrix does not include taxa names as dimnames for rows.")
    }
    x <- A
  }
  if (is.vector(A) == TRUE) {
    if (is.null(names(A))) {
      stop("Data vector does not include taxa names as names.")
    }
    x <- as.matrix(A)
  }
  if (class(phy) != "phylo") 
    stop("tree must be of class 'phylo.'")
  if (!is.binary.tree(phy)) 
    stop("tree is not fully bifurcating.")
  N <- length(phy$tip.label)
  if (N != dim(x)[1]) {
    stop("Number of taxa in data matrix and tree are not not equal.")
  }
  if (length(match(rownames(x), phy$tip.label)) != N) 
    stop("Data matrix missing some taxa present on the tree.")
  if (length(match(phy$tip.label, rownames(x))) != N) 
    stop("Tree missing some taxa in the data matrix.")
  if (any(is.na(match(sort(phy$tip.label), sort(rownames(x))))) == 
      T) {
    stop("Names do not match between tree and data matrix.")
  }
  x <- x[phy$tip.label, ]
  if (is.null(dim(x)) == TRUE) {
    x <- matrix(x, dimnames = list(names(x)))
  }
  if (method == "Kmult") {
    ones <- matrix(1,N)
    x <- as.matrix(x)
    N <- length(phy$tip.label)
    pY <- prep_multipic(x,phy = phy)
    pY_results <- do.call(multipic,pY)
    sumdiagC <- sum(pruningwise.distFromRoot(phy)[1:N])
    sum_invV <- pY_results$sum_invV
    Kmult <- function(x) {
      pY$phe[1:N,] <- x
      pY_results <- do.call(multipic,pY)
      a.obs <- ones %*% pY_results$root
      MSEobs.d <- sum(diag(crossprod(x-a.obs)))
      MSE.d <- sum(diag(crossprod(pY_results$contrasts)))
      K.denom <- (sumdiagC - N / sum_invV)/(N - 1)
      K.stat <- (MSEobs.d/MSE.d)/K.denom
      return(K.stat)
    }
    K.obs <- Kmult(x)
    P.val <- 1
    K.val <- rep(0, iter)
    for (i in 1:iter) {
      x.r <- as.matrix(x[sample(nrow(x)), ])
      rownames(x.r) <- rownames(x)
      K.rand <- Kmult(x.r)
      P.val <- ifelse(K.rand >= K.obs, P.val + 1, P.val)
      K.val[i] <- K.rand
    }
    P.val <- P.val/(iter + 1)
    K.val[iter + 1] = K.obs
    if (ShowPlot == TRUE && dim(x)[2] > 1) {
      plotGMPhyloMorphoSpace(phy, A, ancStates = FALSE)
    }
    return(list(phy.signal = K.obs, pvalue = P.val))
  }
  if (method == "SSC") {
    anc.states <- NULL
    options(warn = -1)
    anc.states <- fasterAnc(phy,x)
    colnames(anc.states) <- NULL
    dist.mat <- as.matrix(dist(rbind(as.matrix(x), as.matrix(anc.states)))^2)
    SSC.o <- sum(dist.mat[phy$edge])
    P.val <- 1
    SSC.val <- rep(0, iter)
    for (ii in 1:iter) {
      x.r <- x[sample(nrow(x)), ]
      if (is.null(dim(x.r)) == TRUE) {
        x.r <- matrix(x.r)
      }
      row.names(x.r) <- row.names(x)
      anc.states.r <- NULL
      options(warn = -1)
      anc.states.r <- fasterAnc(phy,x.r)
      colnames(anc.states.r) <- NULL
      dist.mat.r <- as.matrix(dist(rbind(as.matrix(x.r), 
                                         as.matrix(anc.states.r)))^2)
      SSC.r <- sum(dist.mat.r[phy$edge])
      P.val <- ifelse(SSC.r <= SSC.o, P.val + 1, P.val)
      SSC.val[ii] <- SSC.r
    }
    P.val <- P.val/(iter + 1)
    SSC.val[iter + 1] = SSC.o
    if (ShowPlot == TRUE && dim(x)[2] > 1) {
      plotGMPhyloMorphoSpace(phy, A, ancStates = FALSE)
    }
    return(list(phy.signal = SSC.o, pvalue = P.val))
  }
}

fast.geomorph.compare.multi.evol.rates <-
function (A, gp, phy, Subset = TRUE, method="ML",ShowPlot = TRUE, iter = 999) 
{
  if(method=="REML") REML <- 1 else REML <- 0
  phy <- reorder(multi2di(phy,random=FALSE),"postorder")
  if (any(is.na(A)) == T) {
    stop("Data matrix contains missing values. Estimate these first (see 'estimate.missing').")
  }
  gp <- as.factor(gp)
  if (length(dim(A)) == 3) {
    x <- two.d.array(A)
    p <- dim(A)[1]
    k <- dim(A)[2]
    if (length(gp) != p) {
      stop("Not all landmarks are assigned to a partition.")
    }
    gps <- as.factor(rep(gp, k, each = k, length = p * k))
  }
  if (length(dim(A)) == 2) {
    x <- A
    if (length(gp) != ncol(x)) {
      stop("Not all variables are assigned to a partition.")
    }
    gps <- gp
  }
  ngps <- nlevels(gp)
  if (ngps == 1) {
    stop("Only one shape assigned.")
  }
  ntaxa <- length(phy$tip.label)
  N <- nrow(x)
  if (class(phy) != "phylo") 
    stop("tree must be of class 'phylo.'")
  if (is.null(rownames(x))) {
    stop("Data matrix does not include taxa names.")
  }
  if (N != ntaxa) {
    stop("Number of taxa in data matrix and tree are not not equal.")
  }
  if (length(match(rownames(x), phy$tip.label)) != N) 
    stop("Data matrix missing some taxa present on the tree.")
  if (length(match(phy$tip.label, rownames(x))) != N) 
    stop("Tree missing some taxa in the data matrix.")
  x <- x[phy$tip.label,]
  pY <- prep_multipic(x,phy = phy)
  pY_results <- do.call(multipic,pY)
  sigma.d <- function(phy, x, Subset) {
    p <- ncol(x)
    pY$phe <- pY$phe[,1:p,drop=FALSE]
    pY$phe[1:N,] <- x
    pY$contr <- pY$contr[,1:p,drop=FALSE]
    pY_results <- do.call(multipic,pY)
    sigma <- diag(crossprod(pY_results$contrasts))
    sigma <- mean(sigma)/(N-REML) * if(!Subset) p else 1
    return(sigma = sigma)
  }
  rate.global <- sigma.d(phy, x, Subset)
  rate.gps <- array(NA, ngps)
  for (i in 1:nlevels(gps)) {
    rate.gps[i] <- sigma.d(phy, x[, which(gps == levels(gps)[i])], 
                           Subset)
  }
  rate.ratio <- max(rate.gps)/min(rate.gps)
  if (ngps > 2) {
    sig.rate.gps <- array(1, dim = c(ngps, ngps))
    rate.ratio.gps <- array(0, dim = c(ngps, ngps))
    for (i in 1:(ngps - 1)) {
      for (j in 2:ngps) {
        tmp <- c(rate.gps[i], rate.gps[j])
        rate.ratio.gps[i, j] <- rate.ratio.gps[j, i] <- max(tmp)/min(tmp)
        diag(rate.ratio.gps) <- 0
      }
    }
  }
  p <- ncol(x)
  rate.mat <- crossprod(pY_results$contrasts) / (N-REML)
  diag(rate.mat) <- rate.global
  rate.mat <- matrix(nearPD(rate.mat, corr = FALSE,keepDiag = TRUE)$mat, nrow = p, ncol = p)
  x.sim <- sim.char(phy, rate.mat, nsim = iter)
  sig.rate <- 1
  rate.val <- rep(0, iter)
  for (ii in 1:iter) {
    rate.gps.r <- array(NA, ngps)
    for (i in 1:nlevels(gps)) {
      rate.gps.r[i] <- sigma.d(phy, x.sim[, which(gps == levels(gps)[i]), ii], Subset)
    }
    rate.ratio.r <- max(rate.gps.r)/min(rate.gps.r)
    if (ngps > 2) {
      rate.ratio.gps.r <- array(0, dim = c(ngps, ngps))
      for (i in 1:(ngps - 1)) {
        for (j in 2:ngps) {
          tmp <- c(rate.gps.r[i], rate.gps.r[j])
          rate.ratio.gps.r[i, j] <- rate.ratio.gps.r[j,i] <- max(tmp)/min(tmp)
          diag(rate.ratio.gps.r) <- 0
        }
      }
    }
    sig.rate <- ifelse(rate.ratio.r >= rate.ratio, sig.rate + 1, sig.rate)
    if (ngps > 2) {
      sig.rate.gps <- ifelse(rate.ratio.gps.r >= rate.ratio.gps, 
                             sig.rate.gps + 1, sig.rate.gps)
    }
    rate.val[ii] <- rate.ratio.r
  }
  sig.rate <- sig.rate/(iter + 1)
  if (ngps > 2) {
    sig.rate.gps <- sig.rate.gps/(iter + 1)
  }
  rate.val[iter + 1] = rate.ratio
  if (ShowPlot == TRUE) {
    hist(rate.val, 30, freq = TRUE, col = "gray", xlab = "SigmaD ratio")
    arrows(rate.ratio, 50, rate.ratio, 5, length = 0.1, lwd = 2)
  }
  if (ngps == 2) {
    return(list(rates.all = rate.gps, rate.ratio = rate.ratio, 
                pvalue = sig.rate))
  }
  if (ngps > 2) {
    return(list(rates.all = rate.gps, rate.ratio = rate.ratio, 
                pvalue = sig.rate, pvalue.gps = sig.rate.gps))
  }
}

fasterAnc <- function(tree, x, vars = FALSE, CI = FALSE) 
{
  if (!is.binary.tree(tree)) 
    btree <- multi2di(tree) else btree <- tree
  btree <-reorder(btree,"postorder")
  pY <- prep_multipic(x,phy = btree)
  
  M <- btree$Nnode
  N <- length(btree$tip.label)
  anc <- v <- matrix(0,btree$Nnode,ncol(x))
  for (i in 1:M + N) {
    a <- reorder(multi2di(root(btree, node = i)),"postorder")
    pY$nnode <- a$Nnode
    pY$edge1 <- a$edge[,1]
    pY$edge2 <- a$edge[,2]
    pY$edge_len <- a$edge.length
    pY_results <- do.call(multipic,pY)
    anc[i - N,] <- pY_results$root
    #rownames(anc)[i - N] <- i
    if (vars || CI) {
      v[i - N,] <- 1/pY_results$sum_invV * colMeans(pY_results$contrasts^2)
      #names(v)[i - N] <- names(anc)[i - N]
    }
  }
  if(vars || CI) rownames(v) <- 1:M + N
  rownames(anc) <- 1:M + N
  if (!is.binary.tree(tree)) {
    ancNames <- matchNodes(tree, btree)
    anc <- anc[as.character(ancNames[, 2]),]
    rownames(anc) <- ancNames[, 1]
    if (vars || CI) {
      v <- v[as.character(ancNames[, 2]),]
      rownames(v) <- ancNames[, 1]
    }
  }
  result <- list(ace = anc)
  if (vars) 
    result$var <- v
  if (CI) {
    result$CI95 <- cbind(anc - 1.96 * sqrt(v), anc + 1.96 * 
                           sqrt(v))
    rownames(result$CI95) <- names(anc)
  }
  if (length(result) == 1) 
    return(result$ace)
  else return(result)
}

fast.geomorph.phylo.integration <- function (A1, A2, phy, iter = 999, label = NULL, 
                            verbose = FALSE, ShowPlot = TRUE) 
{
  phy <- reorder(multi2di(phy),"postorder")
  if (any(is.na(A1)) == T) {
    stop("Data matrix 1 contains missing values. Estimate these first(see 'estimate.missing').")
  }
  if (any(is.na(A2)) == T) {
    stop("Data matrix 2 contains missing values. Estimate these first(see 'estimate.missing').")
  }
  if (is.numeric(A1) == FALSE) {
    stop("A1 is not numeric. see ?numeric")
  }
  if (is.numeric(A2) == FALSE) {
    stop("A2 is not numeric. see ?numeric")
  }
  if (class(phy) != "phylo") 
    stop("phy must be of class 'phylo.'")
  if (length(dim(A1)) == 3) {
    x <- two.d.array(A1)
  }
  if (length(dim(A1)) == 2) {
    x <- A1
  }
  if (length(dim(A2)) == 3) {
    y <- two.d.array(A2)
  }
  if (length(dim(A2)) == 2) {
    y <- A2
  }
  num.taxa.X <- nrow(x)
  namesX <- rownames(x)
  num.taxa.Y <- nrow(y)
  namesY <- rownames(y)
  if (is.null(namesX)) {
    stop("No specimen names in data matrix 1. Please assign specimen names.")
  }
  if (length(match(phy$tip.label, namesX)) != num.taxa.X && 
      length(phy$tip.label) < num.taxa.X) 
    stop("Tree is missing some taxa present in the data matrix")
  if (length(match(phy$tip.label, namesX)) != num.taxa.X && 
      num.taxa.X < length(phy$tip.label)) 
    stop("Tree contains some taxa not present in present in the data matrix")
  if (is.null(namesY)) {
    stop("No specimen names in data matrix 2. Please assign specimen names")
  }
  if (is.null(namesX) == FALSE && is.null(namesY) == FALSE) {
    mtch.A <- namesX[is.na(match(namesX, namesY))]
    if (length(mtch.A) > 0) {
      stop("Specimen names in data sets are not the same.")
    }
  }
  mtch.B <- namesX[is.na(match(namesX, phy$tip.label))]
  if (length(mtch.B) > 0) {
    stop("Taxa labels on tree and taxa matrix are not the same.")
  }
  x <- x[phy$tip.label,,drop=FALSE]
  y <- y[phy$tip.label,,drop=FALSE]
  xdims <- 1:ncol(x)
  ydims <- (ncol(x)+1):(ncol(x)+ncol(y))
  
  
  data.all <- cbind(x, y)
  Nspec <- nrow(x)
  
  pALL <- prep_multipic(data.all,phy = phy)
  pALL_results <- do.call(multipic,pALL)
  all <- pALL_results$contrasts
  R <- crossprod(all)/(Nspec-1)
  R12 <- R[xdims, ydims,drop=FALSE]
  pls <- svd(R12)
  U <- pls$u
  V <- pls$v
  XScores <- ((all[,xdims,drop=FALSE])%*%U)[,1]
  YScores <- ((all[,ydims,drop=FALSE])%*%V)[,1]
  Yhat <- sum(XScores*YScores)/sum(XScores^2)*XScores
  pls.obs <- sqrt(sum(Yhat^2)/(sum((YScores-Yhat)^2)+sum(Yhat^2)))
  
  P.val <- 1
  pls.val <- rep(0, iter)
  for (ii in 1:iter) {
    y.r <- y[sample(nrow(y)), ,drop=FALSE]
    data.all.r <- cbind(x, y.r)
    
    pALL$phe[1:Nspec,] <- data.all.r
    pALL_results <- do.call(multipic,pALL)
    all.r <- pALL_results$contrasts
    R.r <- crossprod(all.r)/(Nspec-1)
    R12.r <- R.r[xdims, ydims,drop=FALSE]
    pls.r <- svd(R12.r)
    U.r <- pls.r$u
    V.r <- pls.r$v
    XScores.r <- ((all.r[,xdims])%*%U.r)[,1]
    YScores.r <- ((all.r[,ydims])%*%V.r)[,1]
    Yhat.r <- sum(XScores.r*YScores.r)/sum(XScores.r^2)*XScores.r
    pls.r <- sqrt(sum(Yhat.r^2)/(sum((YScores.r-Yhat.r)^2)+sum(Yhat.r^2)))
    
    
    pls.val[ii] <- pls.r
    P.val <- ifelse(pls.r >= pls.obs, P.val + 1, P.val)
  }
  pls.val[iter + 1] = pls.obs
  P.val <- P.val/(iter + 1)
  if (ShowPlot == TRUE) {
    if (length(dim(A1)) == 2 && length(dim(A2)) == 2) {
      plot(XScores, YScores, pch = 21, bg = "black", 
           main = "Contrast PLS Plot", xlab = "Contrast PLS1 Block 1", ylab = "Contrast PLS1 Block 2")
      if (length(label != 0)) {
        text(XScores, YScores, label, adj = c(-0.7, 
                                              -0.7))
      }
    }
    if (length(dim(A1)) == 3) {
      A1.ref <- mshape(A1)
      pls1.min <- A1[, , which.min(XScores)]
      pls1.max <- A1[, , which.max(XScores)]
    }
    if (length(dim(A2)) == 3) {
      A2.ref <- mshape(A2)
      pls2.min <- A2[, , which.min(XScores)]
      pls2.max <- A2[, , which.max(XScores)]
    }
    if (dim(A1)[2] == 2 || dim(A2)[2] == 2) {
      par(mar = c(1, 1, 1, 1) + 0.1)
      split.screen(matrix(c(0.22, 1, 0.22, 1, 0.19, 0.39, 
                            0, 0.19, 0.8, 1, 0, 0.19, 0, 0.19, 0.19, 0.39, 
                            0, 0.19, 0.8, 1), byrow = T, ncol = 4))
      screen(1)
      plot(XScores, YScores, pch = 21, bg = "black", 
           main = "Contrast PLS1 Plot: Block 1 (X) vs. Block 2 (Y) ", 
           xlab = "Contrast PLS1 Block 1", ylab = "Contrast PLS1 Block 2")
      if (length(label != 0)) {
        text(XScores, YScores, label, adj = c(-0.7, 
                                              -0.7))
      }
      close.screen(all.screens = TRUE)
      par(mar = c(5.1, 4.1, 4.1, 2.1))
    }
    if (length(dim(A1)) == 3 && dim(A1)[2] == 3) {
      plot(XScores, YScores, pch = 21, bg = "black", 
           main = "Contrast PLS Plot", xlab = "Contrast PLS1 Block 1", ylab = "Contrast PLS1 Block 2")
      if (length(label != 0)) {
        text(XScores, YScores, label, adj = c(-0.7, 
                                              -0.7))
      }
      open3d()
      plot3d(pls1.min, type = "s", col = "gray", main = paste("Contrast PLS Block1 negative"), 
             size = 1.25, aspect = FALSE)
      open3d()
      plot3d(pls1.max, type = "s", col = "gray", main = paste("Contrast PLS Block1 positive"), 
             size = 1.25, aspect = FALSE)
    }
    if (length(dim(A2)) == 3 && dim(A2)[2] == 3) {
      open3d()
      plot3d(pls2.min, type = "s", col = "gray", main = paste("Contrast PLS Block2 negative"), 
             size = 1.25, aspect = FALSE)
      open3d()
      plot3d(pls2.max, type = "s", col = "gray", main = paste("Contrast PLS Block2 positive"), 
             size = 1.25, aspect = FALSE)
    }
  }
  if (verbose == TRUE) {
    return(list(`PLS Correlation` = pls.obs, pvalue = P.val, 
                `Block 1 PLS Scores` = XScores, `Block 2 PLS Scores` = YScores[, 
                                                                               1]))
  }
  if (verbose == FALSE) {
    return(list(`PLS Correlation` = pls.obs, pvalue = P.val))
  }
}

fast.geomorph.procD.pgls <- 
function (f1, phy, iter = 999, int.first = FALSE,
          verbose = FALSE) 
{
  data = NULL
  RRPP <- FALSE
  form.in <- formula(f1)
  Y <- eval(form.in[[2]], parent.frame())
  if (length(dim(Y)) == 3) 
    Y <- two.d.array(Y) else Y <- as.matrix(Y)
  form.in <- as.formula(paste(c("Y", form.in[[3]]), collapse = "~"))
  if (int.first == TRUE) 
    ko = TRUE else ko = FALSE
  Terms <- terms(form.in, keep.order = ko)
  k <- length(attr(Terms, "term.labels"))
  mf <- model.frame(form.in)
  if (any(is.na(Y)) == T) {
    stop("Response data matrix (shape) contains missing values. Estimate these first (see 'estimate.missing').")
  }
  if (is.null(dimnames(Y)[[1]])) {
    stop("No species names with Y-data")
  }
  N <- length(phy$tip.label)
  if (length(match(rownames(Y), phy$tip.label)) != N) 
    stop("Data matrix missing some taxa present on the tree.")
  if (length(match(phy$tip.label, rownames(Y))) != N) 
    stop("Tree missing some taxa in the data matrix.")
  Xs = mod.mats2(form.in, mf)
  anova.parts.obs <- anova.parts2(form.in)
  anova.tab <- anova.parts.obs$table
  df <- anova.parts.obs$df[1:k]
  dfE <- anova.parts.obs$df[k + 1]
  
  pY <- prep_multipic(Y,phy)
  pY_results <- do.call(multipic,pY)
  pX_list <- vector("list",k)
  for(i in 1:k)
  {
    temp_X <- Xs[[1]][[i+1]][,-1,drop=FALSE]
    rownames(temp_X) <- rownames(Y)
    pX_list[[i]] <- prep_multipic(temp_X,phy)
  }
  y <- pY_results[[1]]
  SSres <- SStotal <- 0
  SSreg <- double(k)
  pX_results <- vector("list",k)
  for(i in 1:(k))
  {
    pX_results[[i]] <- do.call(multipic,pX_list[[i]])
    x <- pX_results[[i]][[1]]
    XANC <- pX_results[[i]]$root[1,]
    XX <- crossprod(cbind(0,x)) + tcrossprod(c(1,XANC))*pY_results$sum_invV
    betas <- solve(crossprod(x),crossprod(x,y))
    SSres <- sum(diag(crossprod(y-x%*%betas)))
    SStotal <- sum(diag(crossprod(y))) #anova.tab$SS[3]
    SSreg[i] <- SStotal-SSres
  }
  SSres <- SStotal - sum(SSreg)
  Rsq <- SSreg/SStotal
  MS <- SSreg/anova.parts.obs$table$df[1:length(SSreg)]
  MSE <- SSres/anova.parts.obs$table$df[k+1]
  Fs <- MS / MSE
  anova.tab$SS <- c(SSreg,SSres,SStotal)
  anova.tab$MS <- c(MS,MSE,NA)
  anova.tab$Rsq <- c(Rsq,NA,NA)
  anova.tab$F <- c(Fs,NA,NA)
  
  get.stats <- function(pX_list,pY,ind)
  {
    P <- array(, c(k, 1, iter + 1))
    
    for(ii in 1:length(ind))
    {
      pY$phe[1:N,] <- pY$phe[ind[[ii]],,drop=FALSE]
      pY_results <- do.call(multipic,pY)
      y <- pY_results[[1]]
      SSres <- SStotal <- 0
      SSreg <- double(k)
      for(i in 1:(k))
      {
        x <- pX_results[[i]][[1]]
        XANC <- pX_results[[i]]$root[1,]
        XX <- crossprod(cbind(0,x)) + tcrossprod(c(1,XANC))*pY_results$sum_invV
        betas <- solve(crossprod(x),crossprod(x,y))
        SSres <- sum(diag(crossprod(y-x%*%betas)))
        SStotal <- sum(diag(crossprod(y))) #anova.tab$SS[3]
        SSreg[i] <- SStotal-SSres
      }
      SSres <- SStotal - sum(SSreg)
      Rsq <- SSreg/SStotal
      MS <- SSreg/anova.parts.obs$table$df[1:length(SSreg)]
      MSE <- SSres/anova.parts.obs$table$df[k+1]
      Fs <- MS / MSE
      P[,,ii] <- Fs
    }
    P
  }
  ind <- c(list(1:nrow(Y)), (Map(function(x) sample(1:nrow(Y)), 1:iter)))
  P <- get.stats(pX_list,pY,ind)
  P.val <- Pval.matrix2(P)
  Z <- Effect.size.matrix2(P)
  anova.tab <- data.frame(anova.tab, Z = c(Z, NA, NA), P.value = c(P.val, NA, NA))
  
  anova.title = "\nRandomization of Raw Values used\n"
  attr(anova.tab, "heading") <- paste("\nType I (Sequential) Sums of Squares and Cross-products\n", 
                                      anova.title)
  class(anova.tab) <- c("anova", class(anova.tab))
  if (verbose == TRUE) {
    list(anova.table = anova.tab, call = match.call(), SS.rand = P)
  }
  else anova.tab
}

mod.mats2 <- 
function (f1, dat1, keep.order = FALSE) 
{
  Terms <- terms(f1, data = dat1, keep.order = keep.order)
  k <- length(attr(Terms, "term.labels"))
  n <- dim(dat1)[[1]]
  Xs <- as.list(array(0, k + 1))
  Xs[[1]] <- matrix(1, n)
  for (i in 1:k) {
    Xs[[i + 1]] <- model.matrix(Terms[1:i], data = dat1)
  }
  list(Xs = Xs, terms = attr(Terms, "term.labels"))
}


anova.parts2 <- 
function (f1) 
{
  form.in <- formula(f1)
  keep.order <- FALSE
  X <- NULL
  Y <- eval(form.in[[2]], parent.frame())
  Terms <- terms(form.in, keep.order = keep.order)
  mf <- model.frame(Terms)
  if (is.null(X)) {
    Xs <- mod.mats2(f1 = form.in, dat1 = mf, keep.order = keep.order)
  } else {
    Xs = X
  }
  anova.terms <- Xs$terms
  k <- length(Xs$Xs) - 1
  df <- SSEs <- array(0, k + 1)
  df[1] <- 1
  SSY <- SSEs[1] <- 0
  for (i in 1:k) {
    x <- Xs$Xs[[i + 1]]
    df[i + 1] <- qr(x)$rank
  }
  SS.tmp <- c(SSEs[-1], SSEs[k + 1])
  SS <- (SSEs - SS.tmp)[1:(k)]
  SS <- c(SS, SSY - sum(SS), SSY)
  SS <- SS[1:k]
  df.tmp <- c(df[-1], df[k + 1])
  df <- (df.tmp - df)[1:k]
  MS <- SS/df
  R2 <- SS/SSY
  SSE.model <- SSY - sum(SS)
  dfE <- nrow(Y) - (sum(df) + 1)
  MSE <- SSE.model/dfE
  Fs <- MS/MSE
  df <- c(df, dfE, nrow(Y) - 1)
  SS <- c(SS, SSE.model, SSY)
  MS <- c(MS, MSE, NA)
  R2 <- c(R2, NA, NA)
  Fs <- c(Fs, NA, NA)
  a.tab <- data.frame(df, SS, MS, Rsq = R2, F = Fs)
  rownames(a.tab) <- c(anova.terms, "Residuals", "Total")
  list(table = a.tab, B = coef(lm(Y ~ x - 1)), SS = SS, df = df, 
       R2 = R2, F = Fs, Y = Y)
}

Pval.matrix2 <-
function (M) 
{
  P = matrix(0, dim(M)[1], dim(M)[2])
  for (i in 1:dim(M)[1]) {
    for (j in 1:dim(M)[2]) {
      y = M[i, j, ]
      p = pval2(y)
      P[i, j] = p
    }
  }
  if (dim(M)[1] > 1 && dim(M)[2] > 1) 
    diag(P) = 1
  rownames(P) = dimnames(M)[[1]]
  colnames(P) = dimnames(M)[[2]]
  P
}

pval2 <-
function (s) 
{
  p = length(s)
  r = rank(s)[1] - 1
  pv = 1 - r/p
  pv
}

Effect.size.matrix2 <-
function (M, center = FALSE) 
{
  Z = matrix(0, dim(M)[1], dim(M)[2])
  for (i in 1:dim(M)[1]) {
    for (j in 1:dim(M)[2]) {
      y = M[i, j, ]
      n = length(y)
      z = effect.size2(y, center = center) * sqrt((n - 1)/n)
      Z[i, j] = z
    }
  }
  Z
}

effect.size2 <-
function (x, center = FALSE) 
{
  z = scale(x, center = center)
  z[1]
}