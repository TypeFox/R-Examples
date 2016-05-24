predict.classif<-function (object, new.fdataobj = NULL, type = "class", ...) 
{
  if (is.null(new.fdataobj)) 
    return(object$group.est)
  isfdata <- is.fdata(new.fdataobj)
  object$group <- factor(object$group, levels = levels(object$group)[which(table(object$group) > 
                                                                             0)])
  ny <- lev <- levels(object$group)
  if (object$C[[1]] == "classif.tree" | object$C[[1]] == "classif.tree2boost") {
    if (is.null(object)) 
      stop("No classif.tree object entered")
    else {
      beta.l = list()
      if (object$C[[1]] == "classif.tree2boost") {
        object <- object$fit
        new.fdataobj <- list(X = new.fdataobj)
      }
      newx = data = new.fdataobj
      basis.x = object$basis.x
      basis.b = object$basis.b
      formula = object$formula.ini
      tf <- terms.formula(formula)
      terms <- attr(tf, "term.labels")
      nt <- length(terms)
      if (attr(tf, "response") > 0) {
        response <- as.character(attr(tf, "variables")[2])
        pf <- rf <- paste(response, "~", sep = "")
      }
      else pf <- rf <- "~"
      if (attr(tf, "intercept") == 0) {
        print("No intecept")
        pf <- paste(pf, -1, sep = "")
      }
      vtab <- rownames(attr(tf, "factors"))
      vnf = intersect(terms, names(data$df))
      vnf2 = intersect(vtab[-1], names(data$df)[-1])
      vfunc2 = setdiff(terms, vnf)
      vint = setdiff(terms, vtab)
      vfunc = setdiff(vfunc2, vint)
      vnf = c(vnf2, vint)
      off <- attr(tf, "offset")
      kterms = 1
      if (length(vnf) > 0) {
        first = FALSE
        XX = data.frame(data[[1]][, c(vnf2)])
        names(XX) = vnf2
        for (i in 1:length(vnf)) {
          if (kterms > 1) 
            pf <- paste(pf, "+", vnf[i], sep = "")
          else pf <- paste(pf, vnf[i], sep = "")
          kterms <- kterms + 1
        }
        if (attr(tf, "intercept") == 0) {
          print("No intecept")
          pf <- paste(pf, -1, sep = "")
        }
      }
      else first = TRUE
      if (length(vfunc) > 0) {
        for (i in 1:length(vfunc)) {
          if (class(data[[vfunc[i]]])[1] == "fdata") {
            fdataobj <- data[[vfunc[i]]]
            dat <- fdataobj$data
            tt <- fdataobj[["argvals"]]
            if (is.null(rownames(dat))) 
              rownames(dat) <- 1:nrow(dat)
            fdnames = list(time = tt, reps = rownames(dat), 
                           values = "values")
            x.fd <- fdataobj[["data"]]
            tt <- fdataobj[["argvals"]]
            rtt <- fdataobj[["rangeval"]]
            if (object$basis.x[[vfunc[i]]]$type != "pc" & 
                  object$basis.x[[vfunc[i]]]$type != "pls") {
              x.fd = Data2fd(argvals = tt, y = t(fdata.cen(fdataobj, 
                                                           object$mean[[vfunc[i]]])[[1]]$data), 
                             basisobj = object$basis.x[[vfunc[i]]], 
                             fdnames = fdnames)
              r = x.fd[[2]][[3]]
              J <- object$JJ[[vfunc[i]]]
              Z = t(x.fd$coefs) %*% J
              colnames(Z) = colnames(J)
            }
            else {
              name.coef <- paste(vfunc[i], ".", rownames(object$basis.x[[vfunc[i]]]$basis$data), 
                                 sep = "")
              newXcen <- fdata.cen(fdataobj, object$mean[[vfunc[i]]])[[1]]
              if (object$basis.x[[vfunc[i]]]$type == 
                    "pls") {
                if (object$basis.x[[vfunc[i]]]$norm) {
                  sd.X <- sqrt(apply(fdataobj$data, 2, 
                                     var))
                  newXcen$data <- newXcen$data/(rep(1, 
                                                    nrow(newXcen)) %*% t(sd.X))
                }
              }
              Z <- inprod.fdata(newXcen, object$vs.list[[vfunc[i]]])
              colnames(Z) <- name.coef
              XX <- data.frame(Z)
              yp = predict(object = object$fit, newdata = XX, 
                           type = "class")
              group.pred <- (factor(yp, levels = lev))
              return(group.pred)
            }
            if (first) {
              XX = Z
              first = FALSE
            }
            else XX = cbind(XX, Z)
          }
          else {
            if (class(data[[vfunc[i]]])[1] == "fd") {
              if (class(object$basis.x[[vfunc[i]]]) != 
                    "pca.fd") {
                x.fd <- fdataobj <- data[[vfunc[i]]]
                r = x.fd[[2]][[3]]
                J <- object$JJ[[vfunc[i]]]
                x.fd$coefs <- x.fd$coefs - object$mean[[vfunc[i]]]$coefs[, 
                                                                         1]
                Z = t(x.fd$coefs) %*% J
                colnames(Z) = colnames(J)
              }
              else {
                name.coef[[vfunc[i]]] = paste(vfunc[i], 
                                              ".", colnames(object$basis.x[[vfunc[i]]]$harmonics$coefs), 
                                              sep = "")
                data[[vfunc[i]]]$coefs <- sweep(data[[vfunc[i]]]$coefs, 
                                                1, (object$basis.x[[vfunc[i]]]$meanfd$coefs), 
                                                FUN = "-")
                fd.cen <- data[[vfunc[i]]]
                Z <- inprod(fd.cen, object$basis.x[[vfunc[i]]]$harmonics)
                colnames(Z) <- name.coef[[vfunc[i]]]
              }
              if (first) {
                XX = Z
                first = FALSE
              }
              else XX = cbind(XX, Z)
            }
            else stop("Please, enter functional covariate")
          }
        }
      }
      if (!is.data.frame(XX)) 
        XX = data.frame(XX)
      yp = predict(object = object$fit, newdata = XX, type = "class")
      group.pred <- (factor(yp, levels = lev))
      return(group.pred)
    }
  }
  if (object$C[[1]] == "classif.glm") {
    prob <- ngroup <- length(table(object$group))
    prob.group <- array(NA, dim = c(nrow(new.fdataobj[[1]]), 
                                    ngroup))
    if (ngroup == 2) {
      probs <- predict.fregre.glm(object$fit[[1]], newx = new.fdataobj, 
                                  ...)
      yest <- ifelse(probs < 0.5, lev[1], lev[2])
      prob.group[, 1] <- 1 - probs
      prob.group[, 2] <- probs
    }
    else {
      for (i in 1:ngroup) {
        prob.group[, i] <- predict.fregre.glm(object$fit[[i]], 
                                              newx = new.fdataobj, ...)
      }
      yest <- apply(prob.group, 1, which.min)
    }
    group.pred <- (factor(yest, levels = lev))
    pgrup <- prob.group
  }
  if (object$C[[1]] == "classif.gsam") {
    prob <- ngroup <- length(table(object$group))
    prob.group <- array(NA, dim = c(nrow(new.fdataobj[[1]]), 
                                    ngroup))
    if (ngroup == 2) {
      probs <- predict.fregre.gsam(object$fit[[1]], newx = new.fdataobj, 
                                   ...)
      yest <- ifelse(probs < 0.5, lev[1], lev[2])
      prob.group[, 1] <- 1 - probs
      prob.group[, 2] <- probs
    }
    else {
      prob.group <- array(NA, dim = c(nrow(new.fdataobj[[1]]), 
                                      ngroup))
      for (i in 1:ngroup) {
        prob.group[, i] <- predict.fregre.gsam(object$fit[[i]], 
                                               newx = new.fdataobj, ...)
      }
      yest <- apply(prob.group, 1, which.min)
    }
    group.pred <- (factor(yest, levels = lev))
    pgrup <- prob.group
  }
  if (object$C[[1]] == "classif.gkam") {
    prob <- ngroup <- length(table(object$group))
    prob.group <- array(NA, dim = c(nrow(new.fdataobj[[1]]), 
                                    ngroup))
    if (ngroup == 2) {
      probs <- predict.fregre.gkam(object$fit[[1]], newx = new.fdataobj, 
                                   ...)
      yest <- ifelse(probs < 0.5, lev[1], lev[2])
      prob.group[, 1] <- 1 - probs
      prob.group[, 2] <- probs
    }
    else {
      prob.group <- array(NA, dim = c(nrow(new.fdataobj[[1]]), 
                                      ngroup))
      for (i in 1:ngroup) {
        prob.group[, i] <- predict.fregre.gkam(object$fit[[i]], 
                                               newx = new.fdataobj, ...)
      }
      yest <- apply(prob.group, 1, which.min)
    }
    group.pred <- (factor(yest, levels = lev))
    pgrup <- prob.group
  }
  if (object$C[[1]] == "classif.glm2boost") {
    newdata <- object$data
    if (isfdata) 
      newfdata <- list(X = new.fdataobj)
    else newfdata <- new.fdataobj
    prob <- ngroup <- length(table(object$group))
    prob.group <- array(NA, dim = c(nrow(new.fdataobj[[1]]), 
                                    ngroup))
    if (ngroup == 2) {
      probs <- predict.fregre.glm(object$fit[[1]], newx = newfdata, 
                                  ...)
      yest <- ifelse(probs < 0.5, lev[1], lev[2])
      prob.group[, 1] <- 1 - probs
      prob.group[, 2] <- probs
    }
    else {
      ny <- as.numeric(names(table(object$group)))
      prob.group <- array(NA, dim = c(nrow(new.fdataobj[[1]]), 
                                      ngroup))
      for (i in 1:ngroup) {
        obj <- object$fit[[i]]
        prob.group[, i] <- predict.fregre.glm(obj, newx = newfdata, 
                                              ...)
      }
      yest <- apply(prob.group, 1, which.min)
    }
    group.pred <- (factor(yest, levels = lev))
    pgrup <- prob.group
  }
  if (object$C[[1]] == "classif.gkam2boost") {
    newdata <- object$data
    prob <- ngroup <- length(table(object$group))
    prob.group <- array(NA, dim = c(nrow(new.fdataobj[[1]]), 
                                    ngroup))
    if (ngroup == 2) {
      probs <- predict.fregre.gkam(object$fit[[1]], newx = new.fdataobj, 
                                   ...)
      yest <- ifelse(probs < 0.5, lev[1], lev[2])
      prob.group[, 1] <- 1 - probs
      prob.group[, 2] <- probs
    }
    else {
      prob.group <- array(NA, dim = c(nrow(new.fdataobj[[1]]), 
                                      ngroup))
      for (i in 1:ngroup) {
        prob.group[, i] <- predict.fregre.gkam(object$fit[[i]], 
                                               newx = new.fdataobj, ...)
      }
      yest <- apply(prob.group, 1, which.min)
    }
    group.pred <- (factor(yest, levels = lev))
    pgrup <- prob.group
  }
  if (object$C[[1]] == "classif.gsam2boost") {
    newdata <- object$data
    if (isfdata) 
      newfdata <- list(X = new.fdataobj)
    else newfdata <- new.fdataobj
    prob <- ngroup <- length(table(object$group))
    prob.group <- array(NA, dim = c(nrow(new.fdataobj[[1]]), 
                                    ngroup))
    if (ngroup == 2) {
      probs <- predict.fregre.gsam(object$fit[[1]], newx = newfdata, 
                                   ...)
      yest <- ifelse(probs < 0.5, lev[1], lev[2])
      prob.group[, 1] <- 1 - probs
      prob.group[, 2] <- probs
    }
    else {
      prob.group <- array(NA, dim = c(nrow(new.fdataobj[[1]]), 
                                      ngroup))
      for (i in 1:ngroup) {
        prob.group[, i] <- predict.fregre.gsam(object$fit[[i]], 
                                               newx = newfdata, ...)
      }
      yest <- apply(prob.group, 1, which.min)
    }
    group.pred <- (factor(yest, levels = lev))
    pgrup <- prob.group
  }
  if (object$C[[1]] == "classif.np") {
    gg <- 1:nrow(new.fdataobj)
    if (isfdata) {
      nas <- apply(new.fdataobj$data, 1, count.na)
      if (any(nas)) {
        bb <- !nas
        cat("Warning: ", sum(nas), " curves with NA are omited\n")
        new.fdataobj$data <- new.fdataobj$data[bb, ]
        gg <- gg[bb]
      }
      newx <- new.fdataobj[["data"]]
      tt <- new.fdataobj[["argvals"]]
      rtt <- new.fdataobj[["rangeval"]]
    }
    else newx <- as.matrix(new.fdataobj)
    nn <- nrow(new.fdataobj)
    x = object$fdataobj
    y = object$y
    h = object$h.opt
    n = nrow(x)
    nn = nrow(newx)
    np <- ncol(x)
    ny <- levels(y)
    numg = nlevels(y)
    if (is.null(rownames(newx))) 
      rownames(newx) <- 1:nn
    bs = (as <- list())
    C <- object$call
    m <- object$m
    Ker = object$Ker
    par.metric <- list()
    par.metric <- attr(object$mdist, "par.metric")
    parm <- attr(object$mdist, "par.metric")
    a1 <- attr(object$mdist, "call")
    if (isfdata) {
      if (a1 == "semimetric.mplsr") {
        par.metric[["fdata1"]] <- x
        par.metric[["fdata2"]] <- new.fdataobj
        nmdist <- t(do.call(a1, par.metric))
      }
      else {
        par.metric[["fdata2"]] <- x
        par.metric[["fdata1"]] <- new.fdataobj
        nmdist <- do.call(a1, par.metric)
      }
    }
    else {
      par.metric[["x"]] <- new.fdataobj
      par.metric[["y"]] <- x
      nmdist <- do.call(a1, par.metric)
    }
#    object$par.S$cv = FALSE
     object$par.S$tt <- nmdist  
    kmdist = object$type.S(nmdist, h = h, Ker = object$Ker,w=object$par.S$w,
                           cv =FALSE) #SIEMPRE 
     
    pgrup = array(0, dim = c(numg, nn))
    l = array(0, dim = c(nn))
    group.pred = array(0, dim = nn)
    for (j in 1:numg) {
      grup = as.integer(y == lev[j])
      pgrup[j, ] <- kmdist %*% matrix(grup, ncol = 1)
    }
    group.pred <- factor(ny[apply(pgrup, 2, which.max)], 
                         levels = ny)
#print(group.pred)
    pgrup <- t(pgrup)
    group.est <- numeric(nn)
    ty <- "S.KNN"
    if (ty == "S.KNN") {
      for (ii in 1:nn) {
        l = seq_along(pgrup[ii, ])[pgrup[ii, ] == max(pgrup[ii, 
                                                            ], na.rm = T)]
        
        if (length(l) > 1) {
#cat(ii,"-",l)
          #clase del mas cercano
          abc<-which(nmdist[ii, ] == min(nmdist[ii, ], na.rm = T))

        group.est[ii] =y[abc]
        }
     else  group.est[ii] = ny[l[1]]
      }
      group.pred <-  factor(group.est,levels = ny)
    }
  }
  if (type == "class") 
    return(group.pred)
  if (type == "probs") {
    if (isfdata) 
      rownames(pgrup) <- rownames(new.fdataobj$data)
    else rownames(pgrup) <- rownames(new.fdataobj)
    colnames(pgrup) <- levels(object$group)
    return(list(group.pred = group.pred, prob.group = pgrup))
  }
  else stop("type argument should be one of class or probs")
}
##################################################
