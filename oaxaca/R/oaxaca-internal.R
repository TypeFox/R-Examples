# oaxaca - Blinder-Oaxaca Decomposition in R
# Author: Marek Hlavac

.onAttach <- 
function(libname, pkgname) {
  packageStartupMessage("\nPlease cite as: \n")
  packageStartupMessage(" Hlavac, Marek (2016). oaxaca: Blinder-Oaxaca Decomposition in R.")
  packageStartupMessage(" R package version 0.1.3. http://CRAN.R-project.org/package=oaxaca \n")
}

.is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  {
  if (!is.numeric(x)) { return(FALSE) }
  else { return(abs(x - round(x)) < tol) }
}

.is.all.integers <-
  function(x) {
    if (!is.numeric(x)) { return(FALSE) }
    return (all(.is.wholenumber(x))) 
  }

.is.all.integers.at.least <-
  function(x,at.least) {
    if (!is.numeric(x)) { return(FALSE) }
    if (!is.na(match(NA,x))) { return(FALSE) }
    if (.is.all.integers(x)) {
      return(min(x)>=at.least)
    }
    else {
      return(FALSE)
    }
  }

.is.all.between <-
  function(x,lower,upper) {
    if (!is.numeric(x)) { return(FALSE) }
    if (!is.na(match(NA,x))) { return(FALSE) }
    if (is.numeric(x)) {
      return(all((x>=lower) & (x<=upper)))
    }
    else {
      return(FALSE)
    }
  }

.confidence.interval <-
  function(beta.vector, se.vector, level) {
    if (is.null(names(beta.vector)) || is.null(names(se.vector))) { 
      out <- as.data.frame(cbind(beta.vector, se.vector))
    }
    else {
      out <- as.data.frame(cbind(beta.vector, se.vector[match(names(beta.vector), names(se.vector))]))
    }
    names(out) <- c("beta","se")
    
    # use the asymptotic normal distribution
    error <- qnorm((level+1)/2) * out$se
    out$ci.lower <- out$beta - error
    out$ci.upper <- out$beta + error
    
    return(out)
  }

################################## oaxaca ##################################
.oaxaca.wrap <-
function(formula, data, weights = NULL, R = 1000, reg.fun = lm, cl, ...) {
  
  create.object <-
  function(coef, se, cl) {
    object <- list()
    
    object$beta <- coef$beta
    object$call <- cl   # matched call
    object$n <- coef$n  # number of observations
    object$R <- coef$R  # bootstrap replicates
    object$reg <- coef$reg # regressions
    
    ## coefficients + std errors
    w <- nrow(coef$twofold$overall) # how many weights
    
    if (is.null(se)) {  # if se is NULL, no standard errors (instead just NA)
      se <- list()
      
      se$threefold$overall <- rep(NA, times = 3)
      se$threefold$variables <- matrix(NA, nrow = nrow(coef$threefold$variables), ncol = 3)
      names(se$threefold$overall) <- colnames(se$threefold$variables) <- c("se(endowments)", "se(coefficients)", "se(interaction)")
      
      se$twofold$overall <- matrix(NA, nrow = w, ncol = 5)
      colnames(se$twofold$overall) <- c("weight", "se(explained)", "se(unexplained)", "se(unexplained A)", "se(unexplained B)")
      for (i in 1:w) { 
        se$twofold$variables[[i]] <- matrix(NA, nrow = nrow(coef$threefold$variables), ncol = 5) 
        colnames(se$twofold$variables[[i]]) <- c("weight", "se(explained)", "se(unexplained)", "se(unexplained A)", "se(unexplained B)")
      }
    }
    
    # threefold
    order.threefold <- c("coef(endowments)", "se(endowments)", "coef(coefficients)", "se(coefficients)", "coef(interaction)", "se(interaction)")
    object$threefold$overall <- c(coef$threefold$overall, se$threefold$overall)[order.threefold]
    object$threefold$variables <- cbind(coef$threefold$variables, se$threefold$variables)[,order.threefold]
    
    # twofold
    order.twofold <- c("weight", "coef(explained)", "se(explained)", "coef(unexplained)", "se(unexplained)", "coef(unexplained A)", "se(unexplained A)", "coef(unexplained B)", "se(unexplained B)")
    object$twofold$overall <- cbind(coef$twofold$overall, se$twofold$overall[,-1,drop=FALSE])[,order.twofold]
    
    object$twofold$variables <- list()
    for (i in 1:w) {
        object$twofold$variables[[i]] <- cbind(coef$twofold$variables[[i]], se$twofold$variables[[i]][,-1,drop=FALSE])[,order.twofold]
    }
    
    object$x <- coef$x
    object$y <- coef$y   # dependent variable
  
    class(object) <- "oaxaca"
    return(object)
  }
  
  adjust.dummies <-   # assuming first coefficient is on intercept
  function(beta, dummies) {
    
    out <- beta
    
    c <- sum(beta[dummies]) / (length(dummies) + 1)
    out[dummies] <- beta[dummies] - c
    out[1] <- beta[1] + c
    
    return(out)
  }
  
  run.oaxaca.coef <- 
  function(data, indices, formula.reg, formula.reg.collapse, formula.dummy, group.var, y, w.extra, reg.fun, R, ...) {
    
    out <- list()
    
    n <- nrow(data)
    d <- data
  
    # subset data set based on group variable; 
    data.pooled <- d[(d[, group.var] == 0) || (d[, group.var] == 1), ]   # only retain observations where group.var is 0 or 1
    data.A <- data.pooled[d[, group.var] == 0,]    
    data.B <- data.pooled[d[, group.var] == 1,]   
    
    # regress on both groups
    reg.object.A <- suppressMessages(do.call(reg.fun, list(formula=formula.reg, data=data.A, ...)))
    reg.object.B <- suppressMessages(do.call(reg.fun, list(formula=formula.reg, data=data.B, ...)))
    reg.object.pooled <- suppressMessages(do.call(reg.fun, list(formula=formula.reg, data=data.pooled, ...)))
    reg.object.pooled.collapse <- suppressMessages(do.call(reg.fun, list(formula=formula.reg.collapse, data=data.pooled, ...)))
    
    # model frames  
    model.matrix.A <- model.matrix(reg.object.A)
    model.matrix.B <- model.matrix(reg.object.B)
    
    # difference in Y means
    used.y.A <- as.vector(model.frame(reg.object.A)[,y])
    used.y.B <- as.vector(model.frame(reg.object.B)[,y])
    used.y.pooled <- as.vector(model.frame(reg.object.pooled)[,y])
    
    y.mean.A <- mean(used.y.A, na.rm = TRUE)
    y.mean.B <- mean(used.y.B, na.rm = TRUE)
    y.diff <- y.mean.A - y.mean.B
    
    # number of observations used
    n.A <- length(used.y.A); n.B <- length(used.y.B); n.pooled <- length(used.y.pooled)
    proportion.B <- n.B / n.pooled  # proportion of observations that have group.var = 1
      
    # difference in coefficients
    beta.A <- reg.object.A$coefficients
    beta.B <- reg.object.B$coefficients
    beta.pooled <- reg.object.pooled$coefficients
    beta.pooled.collapse <- reg.object.pooled.collapse$coefficients
    
    # adjust coefficients for reference categories
    if (!is.null(formula.dummy)) {
      for (i in 1:length(formula.dummy)) {
    
        reg.object.dummy <- suppressMessages(do.call(reg.fun, list(formula=formula.dummy[[i]], data=data.pooled, ...))) 
        included.dummy.names <- names(reg.object.dummy$coefficients[-1])
      
        beta.A.adj <- adjust.dummies(beta.A, included.dummy.names)
        beta.B.adj <- adjust.dummies(beta.B, included.dummy.names)
        beta.pooled.adj <- adjust.dummies(beta.pooled, included.dummy.names)
        beta.pooled.collapse.adj <- adjust.dummies(beta.pooled.collapse, included.dummy.names)
      }
    }
      
    # difference in mean endowments
    x.mean.A <- apply(model.matrix.A, MARGIN = 2, FUN = mean)
    x.mean.B <- apply(model.matrix.B, MARGIN = 2, FUN = mean)
    
    ###  three-fold decompositon
    out.vars <- list()
    
    # overall
    E <- as.numeric( t(x.mean.A - x.mean.B) %*% beta.B )              # endowment
    C <- as.numeric( t(x.mean.B) %*% (beta.A - beta.B) )              # coefficients
    I <- as.numeric( t(x.mean.A - x.mean.B) %*% (beta.A - beta.B) )   # interaction
    
    # by variable
    if (is.null(formula.dummy)) {
      E.vars <- (x.mean.A - x.mean.B) * beta.B
      C.vars <- x.mean.B * (beta.A - beta.B)
      I.vars <- (x.mean.A - x.mean.B) * (beta.A - beta.B)
    }
    else {
      E.vars <- (x.mean.A - x.mean.B) * beta.B.adj
      C.vars <- x.mean.B * (beta.A.adj - beta.B.adj)
      I.vars <- (x.mean.A - x.mean.B) * (beta.A.adj - beta.B.adj)
    }
        
    # organize results
    out.threefold.overall <- c(E, C, I)
    out.threefold.vars <- cbind(E.vars, C.vars, I.vars)
    
    if (!is.null(formula.dummy)) {
      
      E.base <- E - sum(E.vars, na.rm = TRUE)
      C.base <- C - sum(C.vars, na.rm = TRUE)
      I.base <- I - sum(I.vars, na.rm = TRUE)
      
      base.line.threefold <- c(E.base, C.base, I.base)
      out.threefold.vars <- rbind(out.threefold.vars, base.line.threefold)
      rownames(out.threefold.vars)[length(rownames(out.threefold.vars))] <- "(Base)"
    }
    
    names(out.threefold.overall) <- colnames(out.threefold.vars) <- c("coef(endowments)", "coef(coefficients)", "coef(interaction)")
        
    ## two-fold decomposition
        
    # weights: -1 = pooled regression betas; -2 = pooled regression betas (including the group variable)
    w <- c(0, 1, 0.5, proportion.B, -1, -2) 
    w <- c(w, w.extra)
        
    # various values of weights (w)
        
    out.twofold.overall <- c()
    out.beta.star <- c()
    out.twofold.vars <- list()
        
    for (i in 1:length(w)) {
          
      w.this <- w[i]
          
      # calculate coefficients
      
      # overall
      if (w.this == -1) { beta.star <- beta.pooled }        # pooled regression w/o including variable
      else if (w.this == -2) {  beta.star <- beta.pooled.collapse[names(beta.pooled)] }  # pooled regression including group variable
      else { beta.star <- w.this * beta.A + (1 - w.this) * beta.B }  
        
      delta.A <- beta.A - beta.star    
      delta.B <- beta.B - beta.star          
      
      # overall
      Q <- as.numeric( t(x.mean.A - x.mean.B) %*% beta.star)  # explained
      U <- as.numeric( t(x.mean.A) %*% delta.A - t(x.mean.B) %*% delta.B) # unexplained
      U.A <- as.numeric( t(x.mean.A) %*% delta.A )
      U.B <- as.numeric( -t(x.mean.B) %*% delta.B )
        
      if (!is.null(formula.dummy)) {
        if (w.this == -1) { beta.star <- beta.pooled.adj }        # pooled regression w/o including variable
        else if (w.this == -2) {  beta.star <- beta.pooled.collapse[names(beta.pooled)] }  # pooled regression including group variable
        else { beta.star <- w.this * beta.A.adj + (1 - w.this) * beta.B.adj }  
        
        delta.A <- beta.A.adj - beta.star
        delta.B <- beta.B.adj - beta.star
      }
      
      # by variable
      Q.vars <- (x.mean.A - x.mean.B) * beta.star
      U.vars <- x.mean.A * delta.A - x.mean.B * delta.B
      U.A.vars <- x.mean.A * delta.A
      U.B.vars <- -x.mean.B * delta.B
          
      out.twofold.overall <- rbind(out.twofold.overall, cbind(w.this, Q, U, U.A, U.B))
      out.twofold.vars[[i]] <- cbind(w.this, Q.vars, U.vars, U.A.vars, U.B.vars)
      out.beta.star <- rbind(out.beta.star, c(w.this, beta.star))
      
      if (!is.null(formula.dummy)) {
        Q.base <- Q - sum(Q.vars, na.rm = TRUE)
        U.base <- U - sum(U.vars, na.rm = TRUE)
        U.A.base <- U.A - sum(U.A.vars, na.rm = TRUE)
        U.B.base <- U.B - sum(U.B.vars, na.rm = TRUE)
        
        base.line <- c(w.this, Q.base, U.base, U.A.base, U.B.base)
        out.twofold.vars[[i]] <- rbind(out.twofold.vars[[i]], base.line)
        rownames(out.twofold.vars[[i]])[length(rownames(out.twofold.vars[[i]]))] <- "(Base)"
      }

      colnames(out.twofold.overall) <- colnames(out.twofold.vars[[i]]) <- c("weight", "coef(explained)", "coef(unexplained)", "coef(unexplained A)", "coef(unexplained B)")      
    }
          
    out.threefold <- list(out.threefold.overall, out.threefold.vars)
    out.twofold <- list(out.twofold.overall, out.twofold.vars)
    
    names(out.threefold) <- names(out.twofold) <- c("overall","variables")
    
    ## numbers of observations
    out.n <- list()
    out.n$n.A <- n.A; out.n$n.B <- n.B; out.n$n.pooled <- n.pooled
    
    ## values of y
    out.y <- list()
    out.y$y.A <- y.mean.A
    out.y$y.B <- y.mean.B
    out.y$y.diff <- y.diff
    
    ### differences in means
    out.x <- list()
    out.x$x.mean.A <- x.mean.A
    out.x$x.mean.B <- x.mean.B
    out.x$x.mean.diff <- x.mean.A - x.mean.B
    
    ### betas
    out.beta <- list()
    if (!is.null(formula.dummy)) { 
      out.beta$beta.A <- beta.A.adj 
      out.beta$beta.B <- beta.B.adj 
    } 
    else { 
      out.beta$beta.A <- beta.A 
      out.beta$beta.B <- beta.B 
    }
    out.beta$beta.diff <- out.beta$beta.A - out.beta$beta.B
    out.beta$beta.R <- out.beta.star
    
    ### regression objects
    out.reg <- list()
    out.reg$reg.A <- reg.object.A
    out.reg$reg.B <- reg.object.B
    out.reg$reg.pooled.1 <- reg.object.pooled
    out.reg$reg.pooled.2 <- reg.object.pooled.collapse
    
    # get rid of the call components (ugly due to do.call)
    out.reg$reg.A$call <- out.reg$reg.B$call <- out.reg$reg.pooled.1$call <- out.reg$reg.pooled.2$call <- NULL
    
    ### OUTPUT
    out <- list(out.beta, out.n, R, out.reg, out.threefold, out.twofold, out.x, out.y)
    names(out) <- c("beta","n", "R", "reg", "threefold", "twofold","x","y")
      
    return(out)
  } 
  
  
  run.oaxaca.bootstrap <- 
    function(data, indices, formula.reg, formula.reg.collapse, formula.dummy, group.var, y, w.extra, reg.fun, R, ...) {
      
      sdx <- function(x) cbind(apply(x, MARGIN = 1, FUN = sd, na.rm=T))
      
      indices <- list()
      
      # use these data
      d <- data
      
      # subset data set based on group variable; 
      data.pooled <- d[(d[, group.var] == 0) || (d[, group.var] == 1), ]   # only retain observations where group.var is 0 or 1
      data.A <- data.pooled[d[, group.var] == 0,]    
      data.B <- data.pooled[d[, group.var] == 1,] 
      
      # regress on both groups
      reg.object.A <- suppressMessages(do.call(reg.fun, list(formula=formula.reg, data=data.A, ...)))
      reg.object.B <- suppressMessages(do.call(reg.fun, list(formula=formula.reg, data=data.B, ...)))
      reg.object.pooled <- suppressMessages(do.call(reg.fun, list(formula=formula.reg, data=data.pooled, ...)))
      reg.object.pooled.collapse <- suppressMessages(do.call(reg.fun, list(formula=formula.reg.collapse, data=data.pooled, ...)))
      
      # model frames  
      model.matrix.A <- model.matrix(reg.object.A)
      model.matrix.B <- model.matrix(reg.object.B)
      
      # difference in Y means
      used.y.A <- as.vector(model.frame(reg.object.A)[,y])
      used.y.B <- as.vector(model.frame(reg.object.B)[,y])
      used.y.pooled <- as.vector(model.frame(reg.object.pooled)[,y])
      
      y.mean.A <- mean(used.y.A, na.rm = TRUE)
      y.mean.B <- mean(used.y.B, na.rm = TRUE)
      y.diff <- y.mean.A - y.mean.B
      
      # number of observations used
      n.A <- length(used.y.A); n.B <- length(used.y.B); n.pooled <- length(used.y.pooled)
      proportion.B <- n.B / n.pooled  # proportion of observations that have group.var = 1
      
      # weights: -1 = pooled regression betas; -2 = pooled regression betas (including the group variable)
      w <- c(0, 1, 0.5, proportion.B, -1, -2) 
      w <- c(w, w.extra)
      
      E <- C <- I <- rep(NA, times = R)
      
      E.vars <- C.vars <- I.vars <- NULL
      E.base <- C.base <- I.base <- NULL
      
      Q <- U <- U.A <- U.B <- list()
      Q.vars <- U.vars <- U.A.vars <- U.B.vars <- list()
      Q.base <- U.base <- U.A.base <- U.B.base <- list()
      
      length(Q) <- length(U) <- length(U.A) <- length(U.B) <- length(w)
      length(Q.base) <- length(U.base) <- length(U.A.base) <- length(U.B.base) <- length(w)
      length(Q.vars) <- length(U.vars) <- length(U.A.vars) <- length(U.B.vars) <- length(w)
      
      n <- nrow(data)
      
      message("\nBootstrapping standard errors:")
      indicate.quantiles <- quantile(1:R, 
                                     c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), 
                                     type=1) # indicate progress at deciles
      
      for (r in 1:R) {
        
        if (r %in% indicate.quantiles) {
          message(paste(as.character(r)," / ",as.character(R)," (",
                        as.character(round((r/R)*100,2)),"%)", sep=""))
        }
        
        indices[[r]] <- sample(1:n, size = n, replace=T)
            
        # use these data
        d <- data[indices[[r]], ]
          
        # subset data set based on group variable; 
        data.pooled <- d[(d[, group.var] == 0) || (d[, group.var] == 1), ]   # only retain observations where group.var is 0 or 1
        data.A <- data.pooled[d[, group.var] == 0,]    
        data.B <- data.pooled[d[, group.var] == 1,] 
          
        # regress on both groups
        reg.object.A <- suppressMessages(do.call(reg.fun, list(formula=formula.reg, data=data.A, ...)))
        reg.object.B <- suppressMessages(do.call(reg.fun, list(formula=formula.reg, data=data.B, ...)))
        reg.object.pooled <- suppressMessages(do.call(reg.fun, list(formula=formula.reg, data=data.pooled, ...)))
        reg.object.pooled.collapse <- suppressMessages(do.call(reg.fun, list(formula=formula.reg.collapse, data=data.pooled, ...)))
          
        # model frames  
        model.matrix.A <- model.matrix(reg.object.A)
        model.matrix.B <- model.matrix(reg.object.B)
        
        # difference in coefficients
        beta.A <- reg.object.A$coefficients
        beta.B <- reg.object.B$coefficients
        beta.pooled <- reg.object.pooled$coefficients
        beta.pooled.collapse <- reg.object.pooled.collapse$coefficients
        
        # adjust coefficients for reference categories
        if (!is.null(formula.dummy)) {
          for (i in 1:length(formula.dummy)) {
            
            reg.object.dummy <- suppressMessages(do.call(reg.fun, list(formula=formula.dummy[[i]], data=data.pooled, ...))) 
            included.dummy.names <- names(reg.object.dummy$coefficients[-1])
            
            beta.A.adj <- adjust.dummies(beta.A, included.dummy.names)
            beta.B.adj <- adjust.dummies(beta.B, included.dummy.names)
            beta.pooled.adj <- adjust.dummies(beta.pooled, included.dummy.names)
            beta.pooled.collapse.adj <- adjust.dummies(beta.pooled.collapse, included.dummy.names)
          }
        }
        
        # difference in mean endowments
        x.mean.A <- apply(model.matrix.A, MARGIN = 2, FUN = mean)
        x.mean.B <- apply(model.matrix.B, MARGIN = 2, FUN = mean)
        
        ###  three-fold decompositon
        out.vars <- list()
        
        # overall
        E[r] <- as.numeric( t(x.mean.A - x.mean.B) %*% beta.B )              # endowment
        C[r] <- as.numeric( t(x.mean.B) %*% (beta.A - beta.B) )              # coefficients
        I[r] <- as.numeric( t(x.mean.A - x.mean.B) %*% (beta.A - beta.B) )   # interaction
        
        # by variable
        if (is.null(formula.dummy)) {
          E.vars <- cbind(E.vars, cbind((x.mean.A - x.mean.B) * beta.B)) 
          C.vars <- cbind(C.vars, x.mean.B * (beta.A - beta.B))
          I.vars <- cbind(I.vars, (x.mean.A - x.mean.B) * (beta.A - beta.B))
        }
        else {
          E.vars <- cbind(E.vars, cbind((x.mean.A - x.mean.B) * beta.B.adj)) 
          C.vars <- cbind(C.vars, x.mean.B * (beta.A.adj - beta.B.adj))
          I.vars <- cbind(I.vars, (x.mean.A - x.mean.B) * (beta.A.adj - beta.B.adj))
        }  
        
        # organize results
        out.threefold.overall <- c(sd(E, na.rm=T), sd(C, na.rm=T), sd(I, na.rm=T))
        out.threefold.vars <- cbind(sdx(E.vars), sdx(C.vars), sdx(I.vars))
        
        if (!is.null(formula.dummy)) {
          E.base <- c(E.base, E[r] - sum(E.vars[,r], na.rm = TRUE))
          C.base <- c(C.base, C[r] - sum(C.vars[,r], na.rm = TRUE))
          I.base <- c(I.base, I[r] - sum(I.vars[,r], na.rm = TRUE))
          
          base.line.3fold.boot <- c(sd(E.base, na.rm=T), sd(C.base, na.rm=T), sd(I.base, na.rm=T))
          out.threefold.vars <- rbind(out.threefold.vars, base.line.3fold.boot)
          rownames(out.threefold.vars)[length(rownames(out.threefold.vars))] <- "(Base)"
        }
        names(out.threefold.overall) <- colnames(out.threefold.vars) <- c("se(endowments)", "se(coefficients)", "se(interaction)")
        
        ## two-fold decomposition
        
        out.twofold.overall <- c()
        out.twofold.vars <- list()
        
        # various values of weights (w)
        
        for (i in 1:length(w)) {
          
          w.this <- w[i]
          
          # overall
          if (w.this == -1) { beta.star <- beta.pooled }        # pooled regression w/o including variable
          else if (w.this == -2) {  beta.star <- beta.pooled.collapse[names(beta.pooled)] }  # pooled regression including group variable
          else { beta.star <- w.this * beta.A + (1 - w.this) * beta.B }  
          
          delta.A <- beta.A - beta.star    
          delta.B <- beta.B - beta.star          
          
          # overall
          Q[[i]] <- c(Q[[i]], as.numeric( t(x.mean.A - x.mean.B) %*% beta.star))  # explained
          U[[i]] <- c(U[[i]], as.numeric( t(x.mean.A) %*% delta.A - t(x.mean.B) %*% delta.B)) # unexplained
          U.A[[i]] <- c(U.A[[i]], as.numeric( t(x.mean.A) %*% delta.A ))
          U.B[[i]] <- c(U.B[[i]], as.numeric( -t(x.mean.B) %*% delta.B ))
          
          if (!is.null(formula.dummy)) {
            if (w.this == -1) { beta.star <- beta.pooled.adj }        # pooled regression w/o including variable
            else if (w.this == -2) {  beta.star <- beta.pooled.collapse[names(beta.pooled)] }  # pooled regression including group variable
            else { beta.star <- w.this * beta.A.adj + (1 - w.this) * beta.B.adj }  
            
            delta.A <- beta.A.adj - beta.star
            delta.B <- beta.B.adj - beta.star
          }
          
          Q.vars[[i]] <- cbind(Q.vars[[i]], (x.mean.A - x.mean.B) * beta.star)
          U.vars[[i]] <- cbind(U.vars[[i]], x.mean.A * delta.A - x.mean.B * delta.B)
          U.A.vars[[i]] <- cbind(U.A.vars[[i]], x.mean.A * delta.A)
          U.B.vars[[i]] <- cbind(U.B.vars[[i]], -x.mean.B * delta.B)
          
          if (!is.null(formula.dummy)) {
            Q.base[[i]] <- c(Q.base[[i]], Q[[i]][r] - sum(Q.vars[[i]][,r], na.rm = TRUE))
            U.base[[i]] <- c(U.base[[i]], U[[i]][r] - sum(U.vars[[i]][,r], na.rm = TRUE))
            U.A.base[[i]] <- c(U.A.base[[i]], U.A[[i]][r] - sum(U.A.vars[[i]][,r], na.rm = TRUE))
            U.B.base[[i]] <- c(U.B.base[[i]], U.B[[i]][r] - sum(U.B.vars[[i]][,r], na.rm = TRUE))
          }
                    
        }
        
      }
      
      for (i in 1:length(w)) {
        out.twofold.overall <- rbind(out.twofold.overall, cbind(w[i], sd(Q[[i]], na.rm=T), sd(U[[i]], na.rm=T), sd(U.A[[i]], na.rm=T), sd(U.B[[i]], na.rm=T)))
        
        out.twofold.vars[[i]] <- cbind(w[i], sdx(Q.vars[[i]]), sdx(U.vars[[i]]), sdx(U.A.vars[[i]]), sdx(U.B.vars[[i]]))
        if (!is.null(formula.dummy)) {
          base.line.bootstrap <- c(w[i], 
                                   sd(Q.base[[i]], na.rm=TRUE), 
                                   sd(U.base[[i]], na.rm=TRUE), 
                                   sd(U.A.base[[i]], na.rm=TRUE), 
                                   sd(U.B.base[[i]], na.rm=TRUE))   # add base
          out.twofold.vars[[i]] <- rbind(out.twofold.vars[[i]], base.line.bootstrap)
        }
        
        colnames(out.twofold.vars[[i]]) <- c("weight", "se(explained)", "se(unexplained)", "se(unexplained A)", "se(unexplained B)")      
      }
      
      colnames(out.twofold.overall) <- c("weight", "se(explained)", "se(unexplained)", "se(unexplained A)", "se(unexplained B)")      
      
      out.threefold <- list(out.threefold.overall, out.threefold.vars)
      out.twofold <- list(out.twofold.overall, out.twofold.vars)
      names(out.threefold) <- names(out.twofold) <- c("overall","variables")
      
      ### OUTPUT
      out <- list(out.threefold, out.twofold)
      names(out) <- c("threefold", "twofold")
      return(out)
    }


  
  # ------------------------------------------------------------------------------------------------------------------
  
  ## check that inputs are okay
  
  # set defaults
  error.msg <- NULL
  
  if (!is(formula, "formula")) { error.msg <- c(error.msg,"Argument 'formula' must be a formula of the form: y ~ x1 + x2 + x3 + ... | z  or  y ~ x1 + x2 + x3 + ... | z | d1 + d2 + d3 + ... .\n\n") }
  else {
    formula <- as.Formula(formula)   
    
    y <- as.character(formula)[[2]]   
    
    if (length(formula)[2] ==  2) {
      formula.reg <- formula(formula, lhs = 1, rhs = 1)
      formula.reg.collapse <- formula(formula, lhs = 1, collapse = TRUE)
      group.var <- as.character(formula(formula, lhs = 0, rhs = 2))[[2]]
      formula.dummy <- NULL
    } 
    else if (length(formula)[2] ==  3) {
      formula.reg <- formula(formula, lhs = 1, rhs = 1)
      formula.reg.collapse <- formula(formula, lhs = 1, collapse = TRUE)
      group.var <- as.character(formula(formula, lhs = 0, rhs = 2))[[2]]
      
      formula.dummy <- list(); counter <- 0
      for (i in 3:(length(formula)[2])) {
        counter <- counter + 1
        formula.dummy[[counter]] <- formula(formula, lhs = 1, rhs = i)
      }
    } 
    else { error.msg <- c(error.msg,"Argument 'formula' must be a formula of the form: y ~ x1 + x2 + x3 + ... | z  or  y ~ x1 + x2 + x3 + ... | z | d1 + d2 + d3 + ... .\n") }    
    
    # check if group.var is dummy
    if (is.data.frame(data) && (is.null(error.msg))) {
      if (!((is.logical(data[,group.var])) || (all(data[,group.var] %in% c(1,0)))) ) { error.msg <- c(error.msg,"Variable z in argument 'formula' must be logical (TRUE/FALSE) or numeric (with values of only 0 and 1).\n") }    
    
      # check formula dummies
      if (!is.null(formula.dummy)) {
        collect.variables <- NULL   # collect dummy variables names, to check for duplicates
        for (i in 1:length(formula.dummy)) {
          
          # only indicator variables?
          reg.object.dummy <- suppressMessages(do.call(reg.fun, list(formula=formula.dummy[[i]], data=data, ...))) 
          dummy.frame <- model.frame(reg.object.dummy)[,-1, drop = FALSE]
          collect.variables <- c(collect.variables, colnames(dummy.frame))
            
          if (!((is.logical(dummy.frame)) || (all(unlist(dummy.frame) %in% c(1,0)))) ) { error.msg <- c(error.msg,"Variables d1 + d2 + d3 + ... in argument 'formula' must be logical (TRUE/FALSE) or numeric (with values of only 0 and 1).\n") }
          else {
            # do they all add up to one or zero
            if (!all((apply(dummy.frame, MARGIN = 1, sum, na.rm=TRUE)) %in% c(0,1))) { error.msg <- c(error.msg,"Variables d1 + d2 + d3 + ... in argument 'formula' must indicate membership in mutually exclusive categories.\n") }    
          }
              
          # can the variables be found in the main regression
          reg.object.main <- suppressMessages(do.call(reg.fun, list(formula=formula.reg, data=data, ...))) 
          main.frame <- model.frame(reg.object.main)[,-1, drop = FALSE]
          if (!(all(names(dummy.frame) %in% names(main.frame)))) { error.msg <- c(error.msg,"Variables d1 + d2 + d3 + ... in argument 'formula' must be present in the main regression model.\n") }    
          
        }
        if (length(collect.variables) != length(unique(collect.variables))) { error.msg <- c(error.msg,"Sets of variables d1 + d2 + d3 + ... in argument 'formula' must not overlap.\n") }    
      }
    }
  }
  
  if (!is.data.frame(data)) { error.msg <- c(error.msg,"Argument 'data' must contain a data frame.\n")}
  if ((!is.null(R)) && (!.is.all.integers.at.least(R, 1))) { error.msg <- c(error.msg,"Argument 'R' must be NULL or a positive integer.\n")}
  if ((!is.null(R)) && (length(R) !=1 )) { error.msg <- c(error.msg,"Argument 'R' must be of length 1.\n")}
  if ((!is.null(weights)) && (!.is.all.between(weights, 0, 1))) { error.msg <- c(error.msg,"Argument 'weights' must be NULL or a vector of positive numbers between 0 and 1.\n")}
  if (!is.function(reg.fun)) { error.msg <- c(error.msg,"Argument 'reg.fun' must be a function.\n")}
  if (length(reg.fun) != 1) { error.msg <- c(error.msg,"Argument 'reg.fun' must be of length 1.\n")}

  # stop execution if errors found
  if (!is.null(error.msg)) { 
    message(error.msg) 
    return(NULL)
  }
  

  ##################################
  
  message("oaxaca: oaxaca() performing analysis. Please wait.")
  
  # remove repetitive weights
  if (!is.null(weights)) { weights <- unique(weights) }
  
  # call bunch of functions
  coef <- run.oaxaca.coef(formula.reg = formula.reg, formula.reg.collapse = formula.reg.collapse, formula.dummy = formula.dummy, data = data, group.var = group.var, y = y, w.extra = weights, reg.fun = reg.fun, R=R, ...)
  if (!is.null(R)) { se <- run.oaxaca.bootstrap(formula.reg = formula.reg, formula.reg.collapse = formula.reg.collapse, formula.dummy = formula.dummy, data = data, group.var = group.var, y = y, w.extra = weights, reg.fun = reg.fun, R = R, ...) }
  else { se <- NULL }
  
  # create 'oaxaca' object
  out <- create.object(coef, se, cl)
  return(out)
}

################################## PLOT.OAXACA ##################################
# Blinder-Oaxaca plots]
.plot.oaxaca <- function(x, decomposition, type,
                         w, unexplained.split,
                         variables, components,
                         component.left, 
                         component.labels, variable.labels,
                         ci, ci.level, 
                         title, xlab, ylab, bar.color, cl, ...) {
  
  
  # check that arguments are ok
  error.msg <- NULL
  if (!is.object(x)) { error.msg <- c(error.msg, "Argument 'x' must be an object of class \"oaxaca\".\n")}
  else {
    if (class(x)!="oaxaca") { error.msg <- c(error.msg, "Argument 'x' must be an object of class \"oaxaca\".\n")}
  }
  if (!is.null(error.msg)) { message(error.msg); return(NULL)} # early exit if wrong kind of object
  
  if (length(decomposition) != 1) { error.msg <- c(error.msg,"Argument 'decomposition' must be of length 1.\n")}  
  if ((length(decomposition)==1) && (!(decomposition %in% c("twofold", "threefold")))) { error.msg <- c(error.msg,"Argument 'decomposition' must be either 'threefold' (default) or 'twofold'.\n")}  
  
  if (length(type) != 1) { error.msg <- c(error.msg,"Argument 'type' must be of length 1.\n")}  
  if ((length(type)==1) && (!(type %in% c("overall", "variables")))) { error.msg <- c(error.msg,"Argument 'decomposition' must be either 'overall' or 'variables' (default).\n")}
  
  if ((!is.null(w)) && (!is.numeric(w))) { error.msg <- c(error.msg,"Argument 'weight' must be NULL or numeric.\n")}  
  if ((!is.null(w)) && (length(w) != 1)) { error.msg <- c(error.msg,"Argument 'weight' must be of length 1.\n")}  
  
  # find the appropriate twofold decomposition given w
  w.where <- NULL
  if ((is.numeric(w)) && (length(w)==1)) {
    for (i in 1:length(x$twofold$variables)) {
      if (x$twofold$variables[[i]][1,1] == w) {w.where <- i}
    }
    if (is.null(w.where)) { error.msg <- c(error.msg,"Unable to find the weight specified in argument 'w'.\n")}
  }
  
  if ((decomposition == "twofold") && (is.null(w))) { error.msg <- c(error.msg,"Argument 'weight' must be specified when decomposition = 'twofold'.\n")}
                             
  if (!is.logical(unexplained.split)) { error.msg <- c(error.msg,"Argument 'unexplained.split' must be logical (TRUE / FALSE).\n")}  
  if (length(unexplained.split) != 1) { error.msg <- c(error.msg,"Argument 'unexplained.split' must be of length 1.\n")}                           
  
  if ((!is.character(variables)) && (!is.null(variables))) { error.msg <- c(error.msg,"Argument 'variables' must be NULL or a vector of character strings that contains variable names.\n")}
  if (is.character(variables)) {
    if ((!is.null(variables)) && (!all(variables %in% rownames(x$threefold$variables)))) { error.msg <- c(error.msg,"All variables in argument 'variables' must be included in the decomposition.\n")}
  }
  
  if ((!is.character(components)) && (!is.null(components))) { error.msg <- c(error.msg,"Argument 'components' must be NULL or a vector of character strings that contains component names.\n")}
  if (is.character(components)) {
    all.components.available <- c("endowments","coefficients","interaction","explained","unexplained","unexplained A","unexplained B")
    if ((!is.null(components)) && (!all(components %in% all.components.available))) { error.msg <- c(error.msg,"Unrecognized components in argument 'components'.\n")}
    
    # check if you got the right components
    threefold.components <- c("endowments","coefficients","interaction")
    twofold.components <- c("explained","unexplained")
    twofold.components.split <- c("explained","unexplained A","unexplained B")
    if (decomposition == "threefold") {
      if (!(all(components %in% threefold.components))) { error.msg <- c(error.msg,"Argument 'components' contains components that are not appropriate for the requested decomposition.\n")}      
    }
    else if (decomposition == "twofold") {
      if (unexplained.split == FALSE) {
        if (!(all(components %in% twofold.components))) { error.msg <- c(error.msg,"Argument 'components' contains components that are not appropriate for the requested decomposition.\n")}      
      }
      else {
        if (!(all(components %in% twofold.components.split))) { error.msg <- c(error.msg,"Argument 'components' contains components that are not appropriate for the requested decomposition.\n")}      
      }
    }
    
  }
  
  if (!is.logical(component.left)) { error.msg <- c(error.msg,"Argument 'component.left' must be logical (TRUE / FALSE).\n")}  
  if (length(component.left) != 1) { error.msg <- c(error.msg,"Argument 'component.left' must be of length 1.\n")}             
                             
  if ((!is.character(component.labels)) && (!is.null(component.labels))) { error.msg <- c(error.msg,"Argument 'component.labels' must be NULL or a vector of character strings.\n")}
  if ((!is.character(variable.labels)) && (!is.null(variable.labels))) { error.msg <- c(error.msg,"Argument 'variable.labels' must be NULL or a vector of character strings.\n")}
  
  if (!is.logical(ci)) { error.msg <- c(error.msg,"Argument 'ci' must be logical (TRUE / FALSE).\n")}  
  if (length(ci) != 1) { error.msg <- c(error.msg,"Argument 'ci' must be of length 1.\n")}             
  
  if (!is.numeric(ci.level)) { error.msg <- c(error.msg,"Argument 'ci.level' must be numeric.\n")}  
  if (length(ci.level) != 1) { error.msg <- c(error.msg,"Argument 'ci.level' must be of length 1.\n")}             
  if (is.numeric(ci.level) && (length(ci.level) == 1)) { 
    if (!((ci.level > 0) && (ci.level <= 1))) { error.msg <- c(error.msg,"Argument 'ci.level' must be a number between 0 and 1.\n")} 
  }
  
  if (!is.character(title)) { error.msg <- c(error.msg,"Argument 'title' must be a character string that contains the plot title.\n")}  
  if (length(title) != 1) { error.msg <- c(error.msg,"Argument 'title' must be a single character string.\n")}  
  
  if (!is.character(xlab)) { error.msg <- c(error.msg,"Argument 'xlab' must be a character string that contains the x-axis label.\n")}  
  if (length(xlab) != 1) { error.msg <- c(error.msg,"Argument 'xlab' must be a single character string.\n")}  
  
  if (!is.character(ylab)) { error.msg <- c(error.msg,"Argument 'ylab' must be a character string that contains the y-axis label.\n")}  
  if (length(ylab) != 1) { error.msg <- c(error.msg,"Argument 'ylab' must be a single character string.\n")}  
  
  if ((!is.null(bar.color)) && (!is.vector(bar.color))) { error.msg <- c(error.msg,"Argument 'bar.color' must be NULL or a vector.\n")}  
  
  if (!is.null(error.msg)) { message(error.msg); return(NULL) }
  
  ##################################
  ## defaults
  
  # component labels
  component.labels.default <- c("endowments"="Endowments", "coefficients"="Coefficients", "interaction"="Interaction",
                                "explained"="Explained", "unexplained" = "Unexplained", "unexplained A" = "Unexplained (A)", "unexplained B" = "Unexplained (B)")
  

  # fill in from user-specified values
  component.labels <- c(component.labels, component.labels.default[!(names(component.labels.default) %in% names(component.labels))])
  
  # which variables to get
  if (!is.null(variables)) { variables.use <- variables }
  else (variables.use <- rownames(x$threefold$variables) )
  
  
  ##################################
  
  if (type == "overall") { component.left <- TRUE }   # for overall decomposition, components always on the left
  
  if (decomposition == "threefold") {
    if (type == "variables") {
      coef <- x$threefold$variables[variables.use, c(1,3,5)]
      se <- x$threefold$variables[variables.use, c(2,4,6)]
    }
    else if (type == "overall") {
      coef <- rbind(x$threefold$overall[c(1,3,5)])
      se <- rbind(x$threefold$overall[c(2,4,6)])
    }
  } 
  else if (decomposition == "twofold") {
    if (type == "variables") {
      if (unexplained.split == FALSE) {
        coef <- x$twofold$variables[[w.where]][variables.use, c(2,4)]   
        se <- x$twofold$variables[[w.where]][variables.use, c(3,5)]
      } else {
        coef <- x$twofold$variables[[w.where]][variables.use, c(2,6,8)] 
        se <- x$twofold$variables[[w.where]][variables.use, c(3,7,9)]
      }
    } 
    else if (type == "overall") {
      if (unexplained.split == FALSE) {
        coef <- rbind(x$twofold$overall[w.where, c(2,4)])
        se <- rbind(x$twofold$overall[w.where, c(3,5)])
      } else {
        coef <- rbind(x$twofold$overall[w.where, c(2,6,8)])
        se <- rbind(x$twofold$overall[w.where, c(3,7,9)])
      }
    }
  }
  
  # check if all standard errors are NA: if so, omit confidence intervals
  if (all(is.na(se))) { ci <- FALSE }
  
  # keep only components we want, if so specified
  if (!is.null(components)) {
    
    keep.components.coef <- keep.components.se <- c()
    for (i in 1:length(components)) { 
      keep.components.coef[i] <- paste("coef(",components[i],")",sep="") 
      keep.components.se[i] <- paste("se(",components[i],")",sep="") 
    }
    
    coef <- coef[,keep.components.coef, drop=FALSE]
    se <- se[,keep.components.se, drop=FALSE]
  } 
  
  # diagram
  coef.melt <- melt(coef)
  se.melt <- melt(se)
  se.melt$se <- se.melt$value
  se.only <- as.data.frame(se.melt$se)
  
  # from now on, work with this data set
  d <- cbind(coef.melt, se.only, .confidence.interval(coef.melt$value, se.only$se, ci.level)[,c("ci.lower","ci.upper")])
  
  if (is.null(components)) { d$components <- d$Var2 } else {d$components <- factor(d$Var2, levels = keep.components.coef)}
  d$vars <- factor(d$Var1, levels = variables.use)
  
  # apply components labels
  for (i in 1:length(component.labels)) {
    component.which <- names(component.labels)[i]
    levels(d$components)[levels(d$components) == paste("coef(", component.which, ")", sep="")] <- component.labels[component.which]  
  }

  # apply variable labels
  for (i in 1:length(variable.labels)) {
    variable.which <- names(variable.labels)[i]
    levels(d$vars)[levels(d$vars) == variable.which] <- variable.labels[variable.which]  
    rownames(coef)[rownames(coef) == variable.which] <- variable.labels[variable.which] 
    var.labels.use <- rownames(coef)
  }
  
  # orientation of the diagram
  if (component.left == FALSE) {
    d$left.side <- d$vars
    d$group <- d$components
    left.side.labels <- var.labels.use
    
    if (is.null(components)) {
      # order factors
      if (decomposition == "threefold") {
        d$group <- factor(d$group, levels=c(component.labels["endowments"], component.labels["coefficients"], component.labels["interaction"]))
      }
      else if (decomposition == "twofold") {
        if (unexplained.split == FALSE) { d$group <- factor(d$group, levels=c(component.labels["explained"], component.labels["unexplained"])) }
        else { d$group <- factor(d$group, levels=c(component.labels["explained"], component.labels["unexplained A"], component.labels["unexplained B"])) }
      }
    }
    
  }
  else { 
    d$left.side <- d$components
    d$group <- d$vars
    
    left.side.labels <- factor(d$left.side, levels = component.labels[components] )
    
    if (is.null(components)) {
      # order factors
      if (decomposition == "threefold") {
        left.side.labels <- factor(d$left.side, levels=c(component.labels["endowments"], component.labels["coefficients"], component.labels["interaction"]))
      }
      else if (decomposition == "twofold") {
        if (unexplained.split == FALSE) { left.side.labels <- factor(d$left.side, levels=c(component.labels["explained"], component.labels["unexplained"])) }
        else { left.side.labels <- factor(d$left.side, levels=c(component.labels["explained"], component.labels["unexplained A"], component.labels["unexplained B"])) }
      }
    }
    
    left.side.labels <- levels(left.side.labels)
    
  }
  
  left.side.labels <- left.side.labels[length(left.side.labels):1]  # reverse order for left-size labels
  how.many.groups <- length(unique(d$group))
  
  if (type == "variables") { 
    p.facet_wrap <- facet_wrap(~group,nrow=how.many.groups) 
  }
  else if (type == "overall") { 
    p.facet_wrap <- geom_blank()
  }
  
  if (!is.null(bar.color)) { 
    
    color.vector <- bar.color
    length(color.vector) <- length(levels(d$left.side))
    color.vector[is.na(color.vector)] <- bar.color[length(bar.color)]
    
    if (is.null(names(bar.color))) {
      names(color.vector) <- left.side.labels[length(left.side.labels):1]
    }
  
    p.scale_fill_manual <- scale_fill_manual(values = color.vector)
    
  } 
  else{ # no color scheme specified
    p.scale_fill_manual <- geom_blank()
  }
  
  ### the following assignments are only to appease R CMD CHECK
  ci.lower <- ci.upper <- value <- left.side <- NULL
  ### 
  
  # toggle confidence interval
  if (ci == TRUE) {
    p.geom_point <- geom_point()
    p.geom_errorbar <- geom_errorbar(width=.3, aes(ymin=ci.lower, ymax=ci.upper),position = position_dodge(0.9)) 
  }
  else {
    p.geom_point <- geom_blank()
    p.geom_errorbar <- geom_blank()
  }
  
  # plot the result
  ggplot(d,
    aes(x = left.side, y = value, fill = left.side, ...)) +
    geom_bar(stat = "identity",position = position_dodge(0.9)) +
    p.scale_fill_manual +
    geom_hline(yintercept = 0, linetype = "dotted") + 
    p.geom_point +
    p.geom_errorbar +
    xlab(ylab) +
    ylab(xlab) +
    labs(title = title) +
    coord_flip() +
    p.facet_wrap +
    scale_x_discrete(limits=left.side.labels) + 
    theme_bw() +
    guides(fill=FALSE)

}
