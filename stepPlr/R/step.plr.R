step.plr <- function(x, y, weights = rep(1, length(y)), fix.subset = NULL,
                     level = NULL, lambda = 1e-4, cp = 'bic', max.terms = 5,
                     type = c('both', 'forward', 'forward.stagewise'),
                     trace = FALSE)  
  {
    this.call <- match.call()
    type <- match.arg(type)
    eps <- 1e-5
    nobs <- nrow(x)
    p <- ncol(x)
    if (is.null(dimnames(x)[[2]])) {
      dimnames(x)[[2]] <- xnames <- paste('x', seq(p), sep='')
    } else {
      xnames <- dimnames(x)[[2]]
    }
    if (cp == 'bic') {
      cp <- log(nobs)
    } else if (cp == 'aic') {
      cp <- 2
    }
    max.terms <- min(max.terms, 2^p - 1)
    nfix <- length(fix.subset)
    if (nfix > 0)
      names(fix.subset) <- xnames[fix.subset]
    if (max.terms <= nfix)
      stop('max.terms <= length(fix.subset): Increase max.terms\n')
    action <- action.name <- bestfit <- vector('list')
    dev <- df <- score <- rep(Inf, max.terms + 1)
    group <- vector('list', length = 2)
    ix <- NULL
    if (!is.numeric(x) || !is.null(level)) {
      imat <- get.imat(x, level)
      IX <- imat$imat
      grp <- imat$group
      level <- imat$level
    } else {
      IX <- x
      grp <- seq(p)
    }
    if (trace)
      cat('*** forward addition ***\n')
    ncolx <- 0
    term1 <- seq(max(grp))
    m <- 0
    offset.subset <- offset.coefficients <- NULL
    bestfit[[m + 1]] <- fit <- plr(rep(1,nobs), y, weights, offset.subset,
                                   offset.coefficients, lambda, cp)
    dev[m + 1] <- fit$deviance
    df[m + 1] <- fit$df
    score[m + 1] <- fit$score
    while (m < max.terms) {
      m <- m + 1
      add1 <- FALSE
      if (m > 1 && type == 'forward.stagewise') {
        offset.subset <- seq(ncol(ix))
        offset.coefficients <- bestfit[[m]]$coefficients[-1]
      }
      if (length(term1) > 0) {
        for (term in term1) {
          ix0 <- cbind(ix, IX[ , grp == term, drop = FALSE])
          fit <- plr(ix0, y, weights, offset.subset,
                     offset.coefficients, lambda, cp)
          if (term %in% fix.subset || fit$score < score[m + 1]) {
            bestfit[[m + 1]] <- fit
            ix1 <- ix0
            dev[m + 1] <- fit$dev
            df[m + 1] <- fit$df
            score[m + 1] <- fit$score 
            group[[2]] <- c(group[[1]], rep(m, ncol(ix1) - ncolx))
            action[[m]] <- -term
            action.name[[m]] <- xnames[term]
            add1 <- TRUE
            addterm1 <- term
          }
          if (term %in% fix.subset)
            break
          if (m > 1) {
            for (i in 1:(m - 1)) {
              if (!term.match(xnames[term], action.name[[i]],
                              action.name, m - 1)) {
                ix0 <- cbind(ix, cross.imat(IX[ , grp == term, drop = FALSE],
                                            ix[ , group[[1]] == i,
                                               drop = FALSE]))
                fit <- plr(ix0, y, weights, offset.subset,
                           offset.coefficients, lambda, cp)
                if (fit$score < score[m + 1]) {
                  bestfit[[m + 1]] <- fit
                  ix1 <- ix0             
                  dev[m + 1] <- fit$dev
                  df[m + 1] <- fit$df
                  score[m + 1] <- fit$score
                  group[[2]] <- c(group[[1]], rep(m, ncol(ix1) - ncolx))
                  action[[m]] <- c(-term, i)
                  action.name[[m]] <- c(xnames[term], action.name[[i]])
                  add1 <- FALSE
                }
              }
            }
          }
        }
      }
      if (m > (nfix + 2)) {
        for (i in (nfix + 1):(m - 2)) {
          for (ii in (i + 1):(m - 1)) {
            if (length(action.name[[i]]) == 1 ||
                length(action.name[[ii]]) == 1) {
              if (!term.match(action.name[[i]], action.name[[ii]],
                              action.name, m - 1)) {
                ix0 <- cross.imat(ix[ , group[[1]] == i, drop = FALSE],
                                  ix[ , group[[1]] == ii, drop = FALSE])
                ix0 <- cbind(ix, ix0)
                fit <- plr(ix0, y, weights, offset.subset,
                           offset.coefficients, lambda, cp)
                if (fit$score < score[m + 1]) {
                  bestfit[[m + 1]] <- fit
                  ix1 <- ix0
                  dev[m + 1] <- fit$dev
                  df[m + 1] <- fit$df                  
                  score[m + 1] <- fit$score
                  group[[2]] <- c(group[[1]], rep(m, ncol(ix1) - ncolx))
                  action[[m]] <- c(i, ii)
                  action.name[[m]] <- c(action.name[[i]], action.name[[ii]])
                  add1 <- FALSE
                }
              }
            }
          }
        }
      }                
      if (add1)
        term1 <- term1[term1 != addterm1]
      if (trace) {
        for (i in 1:length(action.name[[m]])) {
          cat(action.name[[m]][i], ' ')
        }
        cat('added\t', 'df=', df[m + 1], '\tscore=', score[m + 1], '\n')
      }
      ix <- ix1
      ncolx <- length(group[[2]])
      group[[1]] <- group[[2]]
    }
    if (type != 'both') {
      if (nfix > 0) {
        min.score <- min(score[-seq(nfix)])
        M <- min(which(score[-seq(nfix)] <= min.score + eps))
        M <- M + nfix
      } else {
        min.score <- min(score)
        M <- min(which(score <= min.score + eps))
      }
      if (M == 1) {
        action <- action.name <- group <- NULL
      } else {
        action <- action[1:(M - 1)]
        action.name <- action.name[1:(M - 1)]
        group <- group[[1]]
        group <- group[group < M]
      }
      object <- list(fit = bestfit[[M]], action = action,
                     action.name = action.name, deviance = dev[1:M],
                     df = nobs - df[1:M], score = score[1:M], group = group,
                     y = y, weights = weights, fix.subset = fix.subset,
                     level = level, lambda = lambda, cp = cp, type = type,
                     xnames = xnames, call = this.call)
      class(object) <- 'stepplr'
      return(object)
    }
    action.f <- action
    action.name.f <- action.name
    action.b <- rep(0, m)
    dev <- c(dev[1], rep(0, m - 1), dev[m + 1])
    df <- c(df[1], rep(0, m - 1), df[m + 1])    
    score <- c(score[1], rep(Inf, m - 1), score[m + 1])
    back.i <- seq(m)
    if (trace)
      cat('\n*** backward deletion ***\n')
    while (m > 1) {
      m0 <- ifelse(m > nfix, nfix + 1, 1)
      for (i in m0:m) {
        delete <- TRUE
        for (ii in m0:m) {
          if (i!=ii &&
              all(action.name.f[[back.i[i]]] %in% action.name.f[[back.i[ii]]]))
            delete <- FALSE
        }
        if (delete) {
          ix0 <- ix[ , group[[1]] != back.i[i]]
          fit <- plr(ix0, y, weights, NULL, NULL, lambda, cp)
          if (fit$score < score[m]) {
            bestfit[[m]] <- fit
            ix1 <- ix0
            dev[m] <- fit$dev
            df[m] <- fit$df
            score[m] <- fit$score
            group[[2]] <- group[[1]][group[[1]] != back.i[i]] 
            action.b[m] <- delete.i <- back.i[i]
          }
        }
      }
      if (trace) {
        for (i in 1:length(action.name.f[[delete.i]])) {
          cat(action.name.f[[delete.i]][i], ' ')
        }
        cat('deleted\t','df=', df[m], '\tscore=', score[m], '\n')
      }
      ix <- ix1
      group[[1]] <- group[[2]]
      back.i <- back.i[back.i != delete.i]
      m <- m - 1
    }
    action.b[1] <- back.i
    if (trace)
      cat(action.name.f[[back.i]], ' deleted\t', 'df=',df[1],
          '\tscore=', score[1], '\n')
    if (nfix > 0) {
      min.score <- min(score[-seq(nfix)])
      M <- min(which(score[-seq(nfix)] <= min.score + eps))
      M <- M + nfix 
    } else {
      min.score <- min(score)
      M <- min(which(score <= min.score + eps))
    }
    if (M == 1) {
      action <- action.name <- group <- NULL
    } else {
      action <- action.name <- vector('list', length = M - 1)
      group <- NULL
      nlevel <- rep(0, p)
      for (i in 1:p) nlevel[i] <- max(1, length(level[[i]]))
      for (i in 1:(M - 1)) {
        act0 <- action.f[[action.b[i]]]
        act1 <- NULL
        s <- 0
        for (a in act0) {
          if (a < 0) {
            act1 <- c(act1, a)
            s <- s + nlevel[-a]
          }
          else {
            ii <- which(action.b == a)
            act1 <- c(act1, ii)
            if (s == 0) {
              s <- sum(group == ii)
            } else {
              s <- s*sum(group == ii)
            }
          }
        }
        action[[i]] <- act1
        action.name[[i]] <- action.name.f[[action.b[i]]]
        group <- c(group, rep(i, s))
      }
    }
    object <- list(fit = bestfit[[M]], action = action,
                   action.name = action.name, deviance = dev[1:M],
                   df = nobs - df[1:M], score = score[1:M], group = group,
                   y = y, weights = weights, fix.subset = fix.subset,
                   level = level, lambda = lambda, cp = cp, type = type,
                   xnames = xnames, call = this.call)
    class(object) <- 'stepplr'
    return(object)
  }

predict.stepplr <- function(object, x = NULL, newx = NULL,
                            type = c('link', 'response', 'class'), ...)
  {
    type <- match.arg(type)
    if (is.null(newx)) {
      pred <- predict(object$fit, type=type, ...)
    } else {
      if (object$type == 'both') {
        if (is.null(x) || nrow(x) != length(object$y))
          stop('x in the training data must be provided')
        if (is.vector(newx))
          newx <- matrix(newx)
        if (is.null(object$action)) {
          ix.tr <- 1
          ix.new <- matrix(1, nrow(newx), 1)
        } else {
          ix.tr <- imat(object, x)
          ix.new <- imat(object, newx)
        }
        fit <- plr(ix.tr, object$y, object$weights,
                   lambda = object$lambda, cp = object$cp)
        pred <- predict(fit, ix.new, type, ...)
      } else {
        ix.new <- imat(object, newx)
        pred <- predict(object$fit, ix.new, type, ...)
      }
    }
    return(pred)
  }

cv.step.plr <- function(x, y, weights = rep(1, length(y)),
                        nfold = 5, folds = NULL, lambda = c(1e-4, 1e-2, 1),
                        cp = c('aic', 'bic'), cv.type=c('deviance', 'class'),
                        trace = TRUE, ...)
  {
    cv.type <- match.arg(cv.type)
    n <- length(y)
    m1 <- length(lambda)
    m2 <- length(cp)
    devresids <- binomial()$dev.resids
    if (is.null(folds)) {
      folds <- split(sample(seq(n)), rep(1:nfold, length = n))
    } else if (length(folds) != nfold) {
      stop('The number of folds must match nfold')
    }
    errorarr <- array(0, dim = c(m1, m2, nfold))
    for (i in 1:nfold) {
      omit <- folds[[i]]
      for (i1 in 1:m1) {
        for (i2 in 1:m2) {
          fit <- step.plr(x[-omit, ], y[-omit], weights[-omit],
                          lambda=lambda[i1], cp=cp[i2], trace=TRUE, ...)
          if (cv.type == 'deviance') {
            mu <- predict(fit, x[-omit, ], x[omit, ], 'response')
            pred <- devresids(y[omit], mu, weights[omit])
            errorarr[i1, i2, i] <- sum(pred)
          } else {
            pred <- predict(fit, x[-omit, ], x[omit, ], 'class')
            errorarr[i1, i2, i] <- mean(pred != y[omit])
          }
        }
      }
      if (trace)
        cat('CV fold', i, '\n')
    }
    error <- apply(errorarr, c(1,2), mean)
    se.error <- sqrt(apply(errorarr, c(1, 2), var) / nfold)
    dimnames(error) <- dimnames(se.error) <- list(lambda, cp)
    return(list(error = error, se.error = se.error,
                lambda = lambda, cp = cp, folds = folds))
  }

anova.stepplr <- function(object, ...)
  {
    cat('\nAnalysis of Deviance Table\n\n')
    if (object$type == 'both') {
      cat('Type: Forward stepwise selection\n\n')
    } else {
      cat('Type: Forward selection\n\n')
    }
    action <- object$action.name
    m <- length(action)
    if (m > 0) {
      for (i in 1:m) {
        k <- length(action[[i]])
        if (k > 1) {
          a <- action[[i]][1]
          for (j in 2:k) {
            a <- paste(a, ':', action[[i]][j], sep = '')
          }
          action[[i]] <- a
        }
      }
      df <- round(object$df, digits=2)
      dev <- round(object$dev, digits=2)
      ddf <- c(NA, df[1:m] - df[2:(m+1)])
      ddev <- c(NA, dev[1:m] - dev[2:(m+1)])
      junk <- cbind(ddf, ddev, df, dev)
    } else {
      junk <- cbind(NA, NA, object$df, object$dev)
    }
    dimnames(junk) <- list(c('Intercept', action),
                           c('Df', 'Deviance', 'Resid.Df', 'Resid.Dev'))
    print(junk)
    return(invisible(junk))
  }

summary.stepplr <- function(object, ...)
  {
    return(invisible(summary(object$fit)))
  }

print.stepplr <- function(x, ...)
  {
    if (x$type == 'both') {
      cat('\nType: Forward stepwise selection\n\n')
    } else {
      cat('\nType: Forward selection\n\n')
    }
    action <- x$action.name
    m <- length(action)
    if (m > 0) {
      for (i in 1:m) {
        k <- length(action[[i]])
        if (k > 1) {
          a <- action[[i]][1]
          for (j in 2:k) {
            a <- paste(a,':', action[[i]][j], sep = '')
          }
          action[[i]] <- a
        }
      }
      cat('Variables with nonzero coefficients:\n')
      for (i in 1:m) {
        cat(action[[i]],'\n')
      }
    } else {
      cat('Null model selected\n')
    }
    fit <- x$fit
    cat('\n    Null deviance:', round(fit$null.deviance, digits = 2),
        'on', fit$nobs - 1, 'degrees of freedom\n')
    cat('Residual deviance:', round(fit$deviance, digits = 2),
        'on', round(fit$nobs - fit$df, digits = 2), 'degrees of freedom\n')
    cat('            Score: deviance +', round(fit$cp, digits = 1),
        '* df =', round(fit$score, digits = 2), '\n')    
  }

