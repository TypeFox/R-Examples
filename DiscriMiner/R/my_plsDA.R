#' @title PLS Discriminant Analysis
#'
#' @description Perform a PLS discriminant analysis
#'
#' @param X matrix or data.frame with explanatory variables
#' @param y vector or factor with group membership
#' @param learn vector of learning observations
#' @param test vector of testing observations
#' @param autosel logical indicating automatic selection of PLS comps
#' @param comps number of PLS components (only when autosel=FALSE)
#' @param cv cross validation method. Options are \code{"LOO"} (Leave-One-Out)
#' and \code{"LKO"} (Leave-K fold-Out)
#' @param k fold left out if using LKO
#' @param retain.models whether to retain lower models (i.e. all lower component
#' results)
#' @keywords internal
my_plsDA <- function(X, y, learn, test, autosel, comps, 
                     cv = "LOO", k = NA, retain.models = FALSE)
{  
  ## prepare ingredients
  ntest = length(test)
  # binarize y
  Y = my_tdc(data.frame(y[learn]))
  glevs = levels(y[learn])
  # dimensions
  n = nrow(X[learn,])
  p = ncol(X)
  q = ncol(Y)
  #added in k for leave-k-out cross validation
  #k=10
  # determine number of PLS components to be computed
  # taking into account rank of X
  Xsvd = svd(X[learn,], nu=0, nv=0)
  rank_X = sum(Xsvd$d > 0.0001)
  if (rank_X == 0)
    stop("\nrank = 0: variables are numerically constant")
  nc = min(n, rank_X)
  if (!autosel) {
    if (comps < nc) nc = comps
  }
  if (nc == n) nc = n - 1
  # standardizing data
  X.old = scale(X[learn,])
  Y.old = scale(Y)
  # creating matrices to store results
  Wh = matrix(0, p, nc)
  Uh = matrix(0, n, nc)
  Th = matrix(0, n, nc)
  Ch = matrix(0, q, nc)
  Ph = matrix(0, p, nc)
  bh = rep(0, nc)
  RSS = rbind(rep(n-1,q), matrix(NA, nc, q))
  PRESS = matrix(NA, nc, q)
  Q2 = matrix(NA, nc, q)
  ## remove random kth out
  if (cv == "LKO") {
    fold = split(sample(1:n), rep(1:k, length=n))
  }
  
  ## PLS2 algorithm
  for (h in 1:nc)
  {
    # "arbitrary" vector (first column of Y.old)
    u.new = Y.old[,1]
    w.old = rep(1, p)
    iter = 1
    repeat
    {
      w.new = t(X.old) %*% u.new / sum(u.new^2)
      w.new = w.new / sqrt(sum(w.new^2))# normalize w.old
      t.new = X.old %*% w.new
      c.new = t(Y.old) %*% t.new / sum(t.new^2)
      u.new = Y.old %*% c.new / sum(c.new^2)
      w.dif = w.new - w.old
      w.old = w.new
      if (sum(w.dif^2)<1e-06 || iter==100) break
      iter = iter + 1
    }
    p.new = t(X.old) %*% t.new / sum(t.new^2)
    
    # Cross validation
    RSS[h+1,] = colSums((Y.old - t.new%*%t(c.new))^2)
    press = matrix(0, n, q)
    
    if (cv != "none")
    {
      ### Random leave-k-out
      if (cv == "LKO") {
        for (i in 1:k)
        { #removes row i, only column 1
          omit=fold[[i]]
          uh.si <- Y.old[-omit,1]
          wh.siold <- rep(1,p)
          itcv <- 1
          repeat
          {
            wh.si <- t(X.old[-omit,]) %*% uh.si / sum(uh.si^2)
            wh.si <- wh.si / sqrt(sum(wh.si^2))
            th.si <- X.old[-omit,] %*% wh.si
            ch.si <- t(Y.old[-omit,]) %*% th.si / sum(th.si^2)
            uh.si <- Y.old[-omit,] %*% ch.si / sum(ch.si^2)
            wsi.dif <- wh.si - wh.siold
            wh.siold <- wh.si
            if (sum(wsi.dif^2)<1e-06 || itcv==100) break
            itcv <- itcv + 1
          }
          Yhat.si = (X.old[omit,] %*% wh.si) %*% t(ch.si)
          press[omit,] = (Y.old[omit,] - Yhat.si)^2
        }
      }
      
      # Leave-One-Out
      if (cv == "LOO") {
        for (i in 1:n)
        {
          uh.si = Y.old[-i,1]
          wh.siold = rep(1,p)
          itcv = 1
          repeat
          {
            wh.si = t(X.old[-i,]) %*% uh.si / sum(uh.si^2)
            wh.si = wh.si / sqrt(sum(wh.si^2))
            th.si = X.old[-i,] %*% wh.si
            ch.si = t(Y.old[-i,]) %*% th.si / sum(th.si^2)
            uh.si = Y.old[-i,] %*% ch.si / sum(ch.si^2)
            wsi.dif = wh.si - wh.siold
            wh.siold = wh.si
            if (sum(wsi.dif^2)<1e-06 || itcv==100) break
            itcv = itcv + 1
          }
          Yhat.si = (X.old[i,] %*% wh.si) %*% t(ch.si)
          press[i,] = (Y.old[i,] - Yhat.si)^2
        }
      }
      PRESS[h,] = colSums(press)
      Q2[h,] = 1 - PRESS[h,]/RSS[h,]
    }
    
    # deflation
    X.old = X.old - (t.new %*% t(p.new))
    Y.old = Y.old - (t.new %*% t(c.new))
    # store new elements
    Wh[,h] = w.new
    Uh[,h] = u.new
    Th[,h] = t.new
    Ch[,h] = c.new
    Ph[,h] = p.new
    bh[h] = t(u.new) %*% t.new
    
    ## selection of PLS components
    # Q2 global
    Q2G = 1 - rowSums(PRESS)/rowSums(RSS[-nc,])
  } # finish PLS algorithm
  
  # automatic selection of PLS components?
  ncs = nc
  if (autosel)
  {
    # Rule 1: Q2G >= 0.05 (Perez & Tenenhaus, 2003)
    selcom = which(Q2G >= 0.05)
    # Rule 2: at least one Q2hk >= 0.095
    #aux = apply(Q2, 1, function(x) sum(x>=0.0975))
    #selcom = which(aux > 0)
    ncs = length(selcom)
    # selecting elements
    Wh = Wh[,selcom]
    Uh = Uh[,selcom]
    Ph = Ph[,selcom]
    Th = Th[,selcom]
    Ch = Ch[,selcom]
  }
  
  ## PLS results
  if (retain.models) {
    mylist.names <- c(paste(seq(nc), "Components", sep="."))
    Br <- vector("list", length(mylist.names))
    cte <- vector("list", length(mylist.names))
    Disc <- vector("list", length(mylist.names))
    coeffs <- vector("list", length(mylist.names))
    names(Br) <- mylist.names
    names(cte) <- mylist.names
    names(Disc) <- mylist.names
    names(coeffs) <- mylist.names
    
    for (i in 1:nc) {
      # weights
      Ws = Wh[,1:i] %*% solve(t(Ph[,1:i])%*%Wh[,1:i])
      # standardized regression coefficients
      Bs = Ws[,1:i] %*% t(Ch[,1:i])
      # regression coeffs non-standardized
      Br[[i]] = diag(1/apply(X[learn,],2,sd)) %*% Bs %*% diag(apply(Y,2,sd))
      cte[[i]] = as.vector((apply(Y,2,mean) - apply(X[learn,],2,mean)%*%Br[[i]]))
      Disc[[i]] = X[test,] %*% Br[[i]] + matrix(rep(cte[[i]],each=ntest), ntest, q)  
      coeffs[[i]] = rbind(INTERCEPT=cte[[i]], Br[[i]])
    }
  } else {
    # weights
    Ws = Wh %*% solve(t(Ph)%*%Wh)
    # standardized regression coefficients
    Bs = Ws %*% t(Ch)
    # regression coeffs non-standardized
    Br = diag(1/apply(X[learn,],2,sd)) %*% Bs %*% diag(apply(Y,2,sd))
    cte = as.vector((apply(Y,2,mean) - apply(X[learn,],2,mean)%*%Br))
    Disc = X[test,] %*% Br + matrix(rep(cte,each=ntest), ntest, q)
    coeffs = rbind(INTERCEPT=cte, Br)
  }
  
  # Q2 global accumulated
  Q2T = cbind(Q2, Q2G)
  q2 = c(paste(rep("Q2",q),colnames(Y),sep="."),"Q2.global")
  # correlations and redundancies
  cor_tx = cor(X[learn,], Th)
  cor_ty = cor(Y, Th)
  R2x = cor(X[learn,], Th)^2 # R2 coefficients
  R2y = cor(Y, Th)^2 # R2 coefficients
  Rdx = colMeans(R2x)
  Rdy = colMeans(R2y) # Sum of squares Y
  R2 = cbind(Rdx, cumsum(Rdx), Rdy, cumsum(Rdy))
  Rd.mat = matrix(0, ncs, ncs)
  for (j in 1:ncs)
    Rd.mat[1:j,j] = Rdy[1:j]
  # variable importance
  VIP = sqrt((Wh^2) %*% Rd.mat %*% diag(p/cumsum(Rdy), ncs, ncs))
  
  ### Weighted VIP for entire fitted model
  weighted.vip <- matrix(0, nrow = p, ncol=ncs)
  for (i in 1:length(Rdy)){
    weighted.vip[,i] <- Rdy[i] * (Wh[,i]^2)
  }
  VIP.weighted <- sqrt(rowSums(weighted.vip) * (p/sum(Rdy) ))
  # combine with individual component VIPs
  VIP <- cbind(VIP, VIP.weighted)
  
  ## adding names
  ### added Ws and Ch for loadings and Y.loadings respectively
  dimnames(Ws) = list(colnames(X), paste(rep("w*",ncs),1:ncs,sep=""))
  dimnames(Ch) = list(colnames(Y), paste(rep("c",ncs),1:ncs,sep=""))
  dimnames(Th) = list(rownames(X[learn,]), paste(rep("t",ncs),1:ncs,sep=""))
  dimnames(Ph) = list(colnames(X), paste(rep("p",ncs),1:ncs,sep=""))
  dimnames(Bs) = list(colnames(X), colnames(Y))
  #dimnames(Br) = list(colnames(X), colnames(Y))
  dimnames(cor_tx) = list(colnames(X), paste(rep("t",ncs),1:ncs,sep=""))
  dimnames(cor_ty) = list(colnames(Y), paste(rep("t",ncs),1:ncs,sep=""))
  dimnames(Q2T) = list(paste(rep("t",nc),1:nc,sep=""), q2)
  dimnames(R2) = list(paste(rep("t",ncs),1:ncs,sep=""),
                      c("R2X","R2Xcum","R2Y","R2Ycum"))
  dimnames(VIP) = list(colnames(X), c(paste(rep("Component ",ncs),1:ncs,sep=""),"Model VIP"))
  #dimnames(VIP.weighted) = list(paste("Model VIP"))
  
  
  # predicted class
  # confusion matrix
  if (retain.models) {
    # for all recursive models
    pred_class <- vector("list", length(mylist.names))
    names(pred_class) <- mylist.names
    conf <- vector("list", length(mylist.names))
    names(conf) <- mylist.names
    for (i in 1:nc) {
      pred_class[[i]] <- factor(max.col(Disc[[i]]), levels=seq_along(glevs), labels=glevs)
      conf = table(original=y[test], predicted=pred_class[[i]])
    }
  } else {
    # for just the highest chosen component
    pred_class = factor(max.col(Disc), levels=seq_along(glevs), labels=glevs) 
    conf = table(original=y[test], predicted=pred_class)
  }
  
  
  # results
  ### added loadings and y.loadings
  res = list(coeffs = coeffs, 
             conf = conf, 
             Disc = Disc, 
             pred_class = pred_class,
             components = Th, 
             loadings = round(Ws,4), 
             y.loadings = round(Ch,4), 
             Q2T = Q2T, 
             R2 = R2, 
             VIP = VIP,
             cor_tx = cor_tx, 
             cor_ty = cor_ty)
  res
}
