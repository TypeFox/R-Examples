term.match <- function(term1, term2, termlist, len)
  {
    if (any(term1 %in% term2))
      return(TRUE)
    if (any(term2 %in% term1))
      return(TRUE)
    term <- c(term1, term2)
    for (m in 1:len) {
      if (length(term) == length(termlist[[m]]) &&
          all(term %in% termlist[[m]]))
        return(TRUE)
    }
    return(FALSE)
  }

get.imat <- function(x, level = NULL)
  {
    n <- nrow(x)
    p <- ncol(x)
    xnames <- dimnames(x)[[2]]
    imat <- label <- group <- NULL
    if (is.null(level)) {
      level <- vector('list',length=p)
      for (j in 1:p) {
        if (is.factor(x[ ,j]))
          level[[j]] <- sort(unique(x[ ,j]))
      }
    }
    names(level) <- xnames    
    for (j in 1:p) {
      m <- max(length(level[[j]]), 1)
      if (m > 1) {
        missing <- !x[ ,j] %in% level[[j]]
        nmiss <- sum(missing)
        iimat <- matrix(0, nrow=n, ncol=m)
        for (jj in 1:m) {
          iimat[x[ ,j] == level[[j]][jj], jj] <- 1
        }
        if (nmiss > 0) {          
          prior <- apply(iimat[!missing, ], 2, sum) / (n - nmiss)
          prior <- rep(prior, rep(nmiss, m))
          iimat[missing, ] <- prior
        }
        lab <- paste(xnames[j], '.', level[[j]], sep = '')
      }
      else {
        if (sum(is.na(x[ ,j])) > 0)
          stop('A continuous factor contains missing observations')
        iimat <- x[ ,j]
        lab <- xnames[j]
      }
      imat <- cbind(imat, iimat)
      label <- c(label, lab)
      group <- c(group, rep(j,m))
    }
    dimnames(imat) <- list(seq(nrow(imat)), label)
    return(list(imat = imat, group = group, level = level))
  }

imat <- function(object, x)
  {
    dimnames(x) <- list(seq(nrow(x)), object$xnames)
    vars <- NULL
    action <- object$action
    level <- object$level
    m <- length(action)
    for (i in 1:m) {
      vars <- c(vars, action[[i]][action[[i]] < 0])
    }
    vars <- unique(vars)
    imat <- get.imat(x[ , -vars, drop = FALSE], level[-vars])
    IX <- imat$imat
    grp <- imat$group
    ix <- group <- NULL
    s <- 0
    for (i in 1:m) {
      act <- action[[i]]
      if (any(act < 0)) {
        a <- act[act < 0]
        ii <- grp == which(vars == a)
        ix0 <- IX[ , ii, drop=FALSE]
      }
      if (any(act > 0)) {
        if (act[1] < 0) {
          ix1 <- ix0
        } else {
          ix1 <- ix[ , group == act[1], drop = FALSE]
        }
        ix2 <- ix[ , group == act[2], drop = FALSE]
        ix0 <- cross.imat(ix1, ix2)
      }
      ix <- cbind(ix, ix0)
      group <- c(group, rep(i, ncol(ix0)))        
    }
    return(ix)
  }

cross.imat <- function(imat1, imat2)
  {
    imat <- label <- NULL
    lab1 <- dimnames(imat1)[[2]]
    lab2 <- dimnames(imat2)[[2]]
    m1 <- ncol(imat1)
    m2 <- ncol(imat2)
    for (i1 in 1:m1) {
      for (i2 in 1:m2) {
        imat <- cbind(imat, imat1[ ,i1] * imat2[ ,i2])
        label <- c(label, paste(lab1[i1], ':', lab2[i2], sep = ''))
      }
    }
    dimnames(imat) <- list(seq(nrow(imat)), label)
    return(imat)
  }
