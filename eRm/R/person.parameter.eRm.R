`person.parameter.eRm` <- function(object){
  # estimation of the person parameters with jml
  # object of class eRm
  # se... whether standard errors should be computed

  se            <- TRUE
  interpolation <- TRUE

  X <- object$X
  #collapse X
  #X.full <- object$X

  max.it <- apply(X, 2L, max, na.rm = TRUE)                           #maximum item raw score without NA
  rp     <- rowSums(X, na.rm = TRUE)                                  #person raw scores
  maxrp  <- apply(X, 1L, function(x.i){ sum(max.it[!is.na(x.i)]) })   #maximum item raw score for person i
  TFrow  <- ((rp==maxrp) | (rp==0))

  pers.exe       <- (1L:nrow(X))[TFrow]       #persons excluded from estimation due to 0/full
  pers.exe.names <- rownames(X)[pers.exe]
  pers.in        <- (1L:nrow(X))[-pers.exe]   #persons in estimation

  if(length(pers.exe) > 0L){   #data matrix persons (full/0) excluded)
    X.ex <- object$X[-pers.exe,]                                        
  } else {
    X.ex <- object$X
  }

  if(any(is.na(X))){
    dichX   <- ifelse(is.na(X), 1, 0)
    strdata <- apply(dichX, 1L, function(x){ paste(x, collapse="") })
    gmemb.X <- as.vector(data.matrix(data.frame(strdata)))
  } else {
    gmemb.X <- rep(1L, nrow(X))
  }

  if(length(pers.exe) > 0L) X <- X[-pers.exe,]

  X.dummy <- X

  if(any(is.na(X))){
    dichX   <- ifelse(is.na(X), 1, 0)
    strdata <- apply(dichX, 1L, function(x){ paste(x, collapse="") })
    gmemb   <- as.vector(data.matrix(data.frame(strdata)))
    gmemb1  <- gmemb
  } else {
    gmemb   <- rep(1L, nrow(X))
    gmemb1  <- gmemb
  }

  mt_vek <- apply(X, 2L, max, na.rm = TRUE)   #number of categories - 1 for each item
  mt_ind <- rep(1L:length(mt_vek), mt_vek)    #index i for items

  indvec <- NULL                              #establish matrix with unique raw scores
  for(i in unique(gmemb)){
    gmemb.ind    <- which(gmemb == i)
    collapse.vec <- which(!duplicated(rowSums(rbind(X[gmemb==i,]), na.rm = TRUE)))
    indvec       <- c(indvec, gmemb.ind[collapse.vec])
  }
  #for (i in unique(gmemb)) indvec <- c(indvec,!duplicated(rowSums(rbind(X[gmemb==i,]),na.rm = TRUE)))

  indvec <- sort(indvec)
  X      <- X[indvec,]     #collapsing X

  beta.all <- object$betapar

  if(!is.null(object$ngroups)){
    if(object$ngroups > 1L) stop("Estimation of person parameters for models with group contrasts not possible!")
  }

  if(is.null(object$mpoints)){
    mpoints <- 1L
  } else {
    mpoints <- object$mpoints
  }

  r.pall <- rowSums(X, na.rm = TRUE)                 #person raw scores


  X01 <- object$X01
  if(length(pers.exe) > 0L) X01 <- X01[-pers.exe,]   #if persons excluded due to 0/full response

  X01           <- X01[indvec,]                      #collapsed version
  gmemb         <- gmemb[indvec]                     #collapsed version
  rownames(X01) <- rownames(X)

  rowvec <- 1L:nrow(X01)



  fitlist <- tapply(rowvec, gmemb, function(rind){      #list with nlm outputs

      if(length(rind) > 1L){
         ivec <- !is.na(X[rind[1],])                       #non-NA elements
         r.i  <- colSums(X[rind, ivec], na.rm = TRUE)      #item raw scores
      } else {                                          #if only one person belongs to raw score group
         ivec <- !is.na(X[rind[1],])
         r.i  <- X[rind, ivec]
        # r.i <- X[rind,]
        # r.i[is.na(r.i)] <- 0
      }
      #r.i     <- colSums(object$X[rind,],na.rm=TRUE)       #item raw scores
      r.p     <- r.pall[rind]                              #person raw scores for current NA group
      X01g    <- rbind(X01[rind,])
      beta    <- beta.all[!is.na(X01g[1L,])]
      X01beta <- rbind(X01g, beta.all)                     #matrix with unique 0/1 response patterns and beta vector in the last row
      theta   <- rep(0L, length(r.p))

      #==================== ML routines ===================================
      jml.rasch <- function(theta){         #fast ML for RM only
        ksi   <- exp(theta)
        denom <- 1/exp(-beta)               #-beta instead of beta since beta are easiness parameter
        lnL   <- sum(r.p*theta) - sum(r.i*(-beta)) - sum(log(1 + outer(ksi,denom)))
        return(-lnL)
      }

      jml <- function(theta){               #ML for all other models
        t1t2.list <- tapply(1L:ncol(X01beta), mt_ind, function(xin){
          #xb <- (t(X01beta)[xin,])
          xb     <- rbind(t(X01beta)[xin,])     #0/1 responses and beta parameters for one item
          beta.i <- c(0.0, xb[,ncol(xb)])       #item parameter with 0

          #person responses (row-wise) on each category for current item
          if((nrow(xb) > 1L) && (length(xin == 1L))){
            x01.i <- as.matrix(xb[,1L:(ncol(xb) - 1L)])
          } else {
            x01.i <- rbind(xb[,1L:(ncol(xb) - 1L)])  #0/1 matrix for item i without beta
          }

          cat0 <- rep(0L, ncol(x01.i))
          cat0[colSums(x01.i) == 0L] <- 1       #those with 0 on the 1-kth category get a 1
          x01.i0 <- rbind(cat0, x01.i)          #appending response vector for 0th category

          ind.h   <- 0L:(length(beta.i)-1L)
          theta.h <- ind.h %*% t(theta)         #n. categories times theta
          #!!!FIXME
          term1 <- (theta.h + beta.i) * x01.i0  #category-person matrix
          t1.i  <- sum(colSums(term1))          #sum over categories and persons
          #print(t1.i)

          term2 <- exp(theta.h+beta.i)
          t2.i  <- sum(log(colSums(term2)))     #sum over categories and persons
          #print(t2.i)

          return(c(t1.i,t2.i))
        })
        termlist <- matrix(unlist(t1t2.list), ncol = 2L, byrow = TRUE)
        termlist <- termlist[!is.na(rowSums(termlist)),]
        st1st2   <- colSums(termlist, na.rm = TRUE) #sum term1, term2

        lnL <- st1st2[1] - st1st2[2]
        return(-lnL)
      }
      #==================== end ML routines ================================

      #==================== call optimizer =================================
      if (object$model == "RM") {
        fit <- nlm(jml.rasch, theta, hessian = se, iterlim = 1000L)
      } else {
        fit <- nlm(jml, theta, hessian = se, iterlim = 1000)
      }
      #fit2 <- optim(theta, jml.rasch, method="BFGS", hessian=TRUE)

      #=================== end call optimizer ==============================
      loglik   <- -fit$minimum
      niter    <- fit$iterations
      thetapar <- fit$estimate
      if(se){
        se <- sqrt(diag(solve(fit$hessian)))
      } else {
        se          <- NA
        fit$hessian <- NA
      }

      return(list("loglik"   = loglik,
                  "niter"    = niter,
                  "thetapar" = thetapar,
                  "se"       = se,
                  "hessian"  = fit$hessian))
  })



  loglik <- niter <- npar <- numeric(length(fitlist))
  thetapar <- se.theta <- hessian <- vector(mode = "list", length = length(fitlist))

  for(i in seq_along(fitlist)){
    loglik[i] <- fitlist[[i]]$loglik
    niter[i]  <- fitlist[[i]]$niter
    npar[i]   <- length(fitlist[[i]]$thetapar)
    thetapar[[i]] <- fitlist[[i]]$thetapar
    se.theta[[i]] <- fitlist[[i]]$se
    hessian[[i]]  <- fitlist[[i]]$hessian
  }



  if(interpolation){   #cubic spline interpolation for 0 and full raw scores
    x     <- rowSums(X, na.rm = TRUE)
    xlist <- split(x, gmemb)

    splineMessage <- FALSE   # to display the spline-failure message only once
    
    max.rs.NAgroups <- lapply(sort(unique(gmemb.X)), function(i){  ## MM 2012-02-01
                         sum(max.it[!is.na(object$X[which(gmemb.X == i)[1],])])   # get the max. raw score per NA-group
                       })
    NAgroups.min <- tapply(rowSums(object$X, na.rm=T), gmemb.X, min, na.rm = TRUE)
    NAgroups.max <- tapply(rowSums(object$X, na.rm=T), gmemb.X, max, na.rm = TRUE)
    
    NAgroup_exclude <- numeric()
    for(i in seq_along(unique(gmemb.X))){   # delete groups with 0/full patterns from max.rs.NAgroups, NAgroups.min, NAgroups.max
      if(all(rowSums(object$X, na.rm=TRUE)[gmemb.X == i] %in% c(0, max.rs.NAgroups[[i]]))) NAgroup_exclude <- c(NAgroup_exclude, i)
    }
    if(length(NAgroup_exclude) > 0L){
      splineMessage <- TRUE   # so that the message is only printed once
      max.rs.NAgroups <- max.rs.NAgroups[-NAgroup_exclude]
      NAgroups.min <- NAgroups.min[-NAgroup_exclude]
      NAgroups.max <- NAgroups.max[-NAgroup_exclude]
      gmemb_reduced <- unique(gmemb.X[-which(gmemb.X %in% NAgroup_exclude)])
    }
    
    pred.list <- Map(function(xx, yy, rs, NAmin, NAmax){   # Map(...)  is  mapply(..., SIMPLFY = FALSE)
                   y    <- tapply(yy, xx, function(xy){ xy[1] })
                   x    <- unique(sort(xx))
                   from <- ifelse(NAmin == 0, 0, min(x))
                   to   <- ifelse(NAmax == rs, rs, max(x))
                   if((length(x) > 3) || (length(y) > 3)){   #otherwise splinereg is not admissible
                     fm1 <- interpSpline(x,y)
                     if((from == 0) | (to == rs)){
                       pred.val <- unclass(predict(fm1, unique(c(from, x, to))))
                     } else {
                       list(x=unname(x),y=unname(y))   #MM2012-02-01
                     }
                   } else {
                     splineMessage <- TRUE   # so that the message is only printed once
                     rval <- list(x=unname(x),y=unname(y))   #MM2012-02-01
                     if(from == 0){  rval$x <- c(0, rval$x)
                                     rval$y <- c(NA, rval$y) }
                     if(to == rs){  rval$x <- c(rval$x, rs)
                                    rval$y <- c(rval$y, NA) }
                     rval
                   }  
                 }, xlist, thetapar, max.rs.NAgroups, NAgroups.min, NAgroups.max)
    X.n <- object$X
  #  if (any(sapply(pred.list,is.null)))  pred.list <- NULL    #no spline interpolation applicable   #MM2012-02-01

    if(splineMessage) message("Spline interpolation in some subgroups not performed!\n  Less than 4 different person parameters estimable!\n  Perhaps NAs in subgroup(s).")
  }



  names(thetapar) <- names(se.theta) <- paste("NAgroup", 1L:length(thetapar), sep="")

  #---------expand theta and se.theta, labeling -------------------
  for(i in unique(gmemb)){
    o.r        <- rowSums(rbind(X.dummy[gmemb1==i,]), na.rm = TRUE)             #orginal raw scores
    names(o.r) <- rownames(X.dummy)[gmemb1 == i]
    c.r       <- rowSums(rbind(X[gmemb==i,]), na.rm = TRUE)                     #collapsed raw scores
    match.ind <- match(o.r, c.r)
    thetapar[[i]] <- thetapar[[i]][match.ind]           #de-collapse theta's
    se.theta[[i]] <- se.theta[[i]][match.ind]           #de-collapse se's
    names(thetapar[[i]]) <- names(se.theta[[i]]) <- names(o.r)
  }
  #--------------- end expand, labeling ---------------------------




  #---------------------- theta.table new ----------------------   ## MM2012-02-01 START
  #thetavec <- unlist(thetapar)
  #ind.orig <- as.vector(unlist(tapply(1:length(gmemb1), gmemb1, function(ind) {ind})))
  #theta.orig <- tapply(thetavec, ind.orig, function(ii) return(ii))   #original order theta parameter
  #theta.table <- data.frame(theta.orig, gmemb1)
  #colnames(theta.table) <- c("Person Parameter","NAgroup")
  #rownames(theta.table) <- rownames(X.ex)

  if(length(NAgroup_exclude) > 0L){
    selector <- gmemb.X %in% gmemb_reduced
  } else {
    selector <- rep(TRUE, length(gmemb.X))
  }
  gmemb.X_final <- gmemb.X[selector]
  while(any(diff(sort(unique(gmemb.X_final))) > 1L) | (min(gmemb.X_final) > 1L)){
    if(!any(gmemb.X_final == 1L)){ gmemb.X_final <- gmemb.X_final - 1L; next }
    for(i in 1:max(gmemb.X_final)){
      if(sum(gmemb.X_final == i) == 0) gmemb.X_final[gmemb.X_final == (i+1L)] <- gmemb.X_final[gmemb.X_final == (i+1L)] - 1L
    }
  }

  theta.table <- data.frame("Person Parameter" = rep(NA, nrow(object$X)),
                            "NAgroup" = NA,
                            "Interpolated" = FALSE,
                            row.names=rownames(object$X),
                            check.names = FALSE)

  theta.table[selector, "Person Parameter"] <- mapply(function(rs, NAgroup){
                                                 pred.list[[NAgroup]]$y[which(pred.list[[NAgroup]]$x == rs)]
                                               }, rowSums(object$X, na.rm=T)[selector], gmemb.X_final)
  theta.table[selector, "NAgroup"] <- gmemb.X_final
                            
  if(length(pers.exe) > 0) theta.table[pers.exe,"Interpolated"] <- TRUE
                                                                   ## MM2012-02-01 END
  

  result <- list(X = X.n, X01 = object$X01, X.ex = X.ex, W = object$W, model = object$model,
                 loglik = loglik, loglik.cml = object$loglik, npar = npar, iter = niter, betapar = object$betapar,
                 thetapar = thetapar, se.theta = se.theta, theta.table = theta.table,
                 pred.list = pred.list, hessian = hessian, mpoints = mpoints,
                 pers.ex = pers.exe, gmemb = gmemb1)
  class(result) <- "ppar"

  return(result)

}
