################
################
classif.npp<-function (group,x, h = NULL, Ker = AKer.norm, metric,
          type.CV = GCV.S,
          type.S = S.NW, 
          par.metric=list(),
          par.CV = list(trim = 0), 
          par.S = list(), ...) 
{
  lfdata<-x
  fdataobj<-lfdata[[1]]
  y <- group
  if (missing(metric))  {
    if (is.fdata(lfdata[[1]])) metric="metric.ldata"
    else metric="metric.dist"
  }
  C <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("y", "lfdata", "h", "Ker", "metric", "type.CV", 
               "type.S", "par.metric","par.CV", "par.S"), names(mf), 0L)
  x<-lfdata[[1]]$data
  n = nrow(x)
  np <- ncol(x)
# print(dim(x))  
  if (n != (length(y))) 
    stop("ERROR IN THE DATA DIMENSIONS")
  if (is.null(rownames(x))) 
    rownames(x) <- 1:n
  if (is.null(colnames(x))) 
    colnames(x) <- 1:np
  types = FALSE
  if (is.matrix(metric)) {
    mdist <- metric
    metric <- attributes(mdist)
  }
  else {
    par.metric$lfdata<-lfdata
    par.metric$lfdataref<-lfdata
# print("calculando distancias")    
    mdist = do.call(metric,par.metric)
  }
#print("sale"  )
  #si se pasa ua matriz classif.np==classif.npp
  #sino metric.ldata  y el resto debe ser igual q la funcion classif.np
  
  ty <- deparse(substitute(type.S))
  if (is.null(h)) 
    h = h.default(fdataobj, metric = mdist, Ker = Ker, type.S = ty, 
                  ...)
  else {
    if (any(h <= 0)) 
      stop("Error: Invalid range for h")
  }
  lenh <- length(h)
  gcv = (cv.error <- array(NA, dim = c(lenh)))
  par.S2 <- par.S
  lenh = length(h)
  if (!is.factor(group)) 
    group <- as.factor(group)
  group <- y <- factor(group, levels = levels(group)[which(table(group) > 
                                                             0)])
  ny <- levels(y)
  numg = nlevels(y)
  Y = array(0, dim = c(numg, n))
  group.est2 = group.est = array(0, dim = c(lenh, n))
  pgrup2 = array(Inf, dim = c(n, lenh))
  pgrup = array(0, dim = c(numg, n, lenh))
  misclassification = array(1, dim = c(1, lenh))
  pr <- 1
  if (is.null(par.S2$h)) 
    par.S$h <- h
  if (is.null(par.S$Ker) & ty != "S.KNN") 
    par.S$Ker <- Ker

  for (i in 1:lenh) {
    par.S$tt <- mdist
    par.S$h <- h[i]
#    if (is.null(par.S$cv) par.S$cv= TRUE
    H = do.call(ty, par.S)
    for (j in 1:numg) {
      Y[j, ] = as.integer(y == ny[j])
      pgrup[j, , i] <- (H %*% matrix(Y[j, ], ncol = 1))
    }
    if (ty == "S.KNN") {
      for (ii in 1:n) {
        l = seq_along(pgrup[, ii, i])[pgrup[, ii, i] == 
                                        max(pgrup[, ii, i], na.rm = T)]
        if (length(l) > 1) {
          l <- y[seq_along(mdist[ii, ])[mdist[ii, ] == 
                                          min(mdist[ii, ], na.rm = T)]]
        }
        group.est[i, ii] = ny[l[1]]
      }
    }
    else {
      group.est[i, ] <- ny[as.vector(apply(pgrup[, , i], 
                                           2, which.max))]
    }
    gcv[i] = sum(group.est[i, ] != y)/n
    if (pr > gcv[i]) {
      pr = gcv[i]
      iknn = i
      prob = 1 - pr
      prob.group2 = t(pgrup[, , i])
      group.pred = group.est[i, ]
    }
  }
  colnames(prob.group2) <- ny
  rownames(prob.group2) <- rownames(x)
  l = which.min(gcv)
  h.opt <- h[l]
  par.S$h <- h.opt
  if (is.null(par.S$cv)) par.S$cv=FALSE
  if (is.null(par.S$w)) par.S$w<-NULL
  if (h.opt == min(h) & par.fda.usc$warning) 
    cat(" Warning: h.opt is the minimum value of bandwidths\n   provided, range(h)=", 
        range(h), "\n")
  else if (h.opt == max(h) & par.fda.usc$warning) 
    cat(" Warning: h.opt is the maximum value of bandwidths\n   provided, range(h)=", 
        range(h), "\n")
  df = traza(H)
  names(gcv) <- h
  group.pred <- factor(group.pred, levels = ny)
  misclass = sum(group.pred != y)/n
  prob.classification <- diag(table(y, group.pred))/table(y)
  out <- list(C = C, group.est = group.pred, group = y, H = H, 
              df = df, y = y,x =lfdata, mdist = mdist, Ker = Ker, 
              metric = metric, type.S = type.S, par.S = par.S, gcv = gcv, 
              h.opt = h.opt, h = h, prob.group = prob.group2, m = m, 
              pgrup = pgrup, pgrup2 = pgrup2, ty = ty, prob.classification = prob.classification, 
              max.prob = 1 - misclass)
  class(out) = "classif"
  return(out)
}
################




