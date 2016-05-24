
### convert format of parameter from numeric vector to list

paramVect2List <- function(paramVect, cl, iso)
{
  if (!(is.numeric(paramVect)))
    stop("'parameter' is not a numeric vector")
  names(paramVect) <- NULL
  ncl <- length(cl)
  lcl <- sapply(cl, length)
  lcl[which(iso)] <- 1
  if (!(length(paramVect) == ncl+sum(lcl)))
    stop("'parameter' does not fit to 'cl'")
  posTheta <- ncl+c(0, cumsum(lcl))
  paramList <- list()
  for (i in 1:ncl)
      paramList[[i]] <- list("alpha"=paramVect[i], 
        "theta"=paramVect[(posTheta[i]+1):posTheta[i+1]])
  paramList
}

### convert format of parameter from list to numeric vector

paramList2Vect <- function(paramList, cl, iso)
{
  if (!(is.list(paramList) & all(sapply(paramList, length)==2)))
    stop("'parameter' is not a list of pairs of lists")
  ncl <- length(paramList)
  if (!(ncl == length(cl)))
    stop("'parameter' must be list of same length as 'cl'")
  for (i in 1:ncl)
    names(paramList[[i]]) <- c("alpha", "theta")
  if (!(all(sapply(paramList, function(x) length(x$alpha))==1)))
    stop("in 'parameter': 'alpha' must be of length 1 for each clique")
  if (!(sum(sapply(paramList, function(x) x$alpha))==1))
    stop("sum of alpha in 'parameter' must be equal to 1")
  lcl <- sapply(cl, length)
  lcl[which(iso)] <- 1
  if (!(all(sapply(paramList, function(x) length(x$theta)) == lcl)))
    stop("in 'parameter': length of 'theta' must be corresponding to 'cl' and 'iso' in each clique")
  paramVect <- unlist(paramList)
  posAlpha <- c(cumsum(c(0,lcl[-ncl]) + 1))
  paramVect <- c(paramVect[posAlpha], paramVect[-posAlpha])
  names(paramVect)[1:ncl] <- paste("alpha",1:ncl,sep="")
  names(paramVect)[-(1:ncl)] <- paste("theta",1:sum(lcl), sep="")
  paramVect
}

### Likelihood for constrained optimization

LogLConstrained <- function(parameter, y, x, DM, n, covtype, 
    eps.R = 1e-12, n.cl, cl, LogLik = "ML", iso) {
    parameter <- c(parameter[1:(n.cl - 1)], 1 - sum(parameter[1:(n.cl - 
        1)]), parameter[(n.cl):length(parameter)])
    alpha <- parameter[1:n.cl]
    theta <- parameter[-(1:n.cl)]
    R <- Rfunc(as.matrix(x), theta, alpha, covtype, n, n.cl, cl, iso) + 
        diag(eps.R, ncol = n, nrow = n)
    F <- matrix(rep(1, n), ncol = 1)
    T <- chol(R)
    M <- backsolve(t(T), F, upper.tri = FALSE)
    Tinv_y <- backsolve(t(T), matrix(y, ncol = 1), upper.tri = FALSE)
    Q <- qr.Q(qr(M))
    H <- Q %*% t(Q)
    z <- Tinv_y - H %*% Tinv_y 
    v <- t(z) %*% z/n
    loglik <- -0.5 * (n * log(2 * pi * v) + 2 * sum(log(diag(T))) + 
        n)
    return(loglik)
}

### Correlation function between x and x

Rfunc <- function(x, theta, alpha, covtype, n, n.cl, cl, iso) {
    R <- matrix(0, ncol = n, nrow = n)
    theta.list <- cl
    ntemp <- 1
    # building up theta.list (theta in clique structure)
    n.cl.ani <- n.cl - sum(iso)
    if (n.cl.ani>0) { 
      # for (j in (1:n.cl.ani)) {  # WRONG: DOES NOT RESPECT THE ORDER
      for (j in (1:n.cl)[!iso]) {
          theta.list[[j]] <- theta[(ntemp):(ntemp -1 + length(cl[[j]]))]
          ntemp <- ntemp + length(cl[[j]])
      }
    }
    if (n.cl.ani < n.cl) {
        # for (j in ((n.cl.ani + 1):n.cl)) {   #WRONG AGAIN
        for (j in (1:n.cl)[iso]) {
            theta.list[[j]] <- theta[ntemp]
            ntemp <- ntemp + 1
        }
    }
    for (j in 1:n.cl) {
        cor.str <- covStruct.create(covtype = covtype, d = length(cl[[j]]), 
            var.names = NULL, known.covparam = "None", coef.cov = theta.list[[j]], 
            coef.var = alpha[j], iso = iso[j])
        R <- R + covMatrix(object = cor.str, X = x[, cl[[j]],drop=FALSE])[[1]]
    }
    return(R)
}

### Correlation function between newdata and x

rfunc <- function(newdata, x, theta, alpha, covtype, n.cl, cl, iso) {
    theta.list <- cl
    ntemp <- 1
    n.cl.ani <- n.cl - sum(iso)
    if (n.cl.ani>0){
      # for (j in (1:n.cl.ani)) {  WRONG: DOES NOT RESPECT THE ORDER
      for (j in (1:n.cl)[!iso]) {
          theta.list[[j]] <- theta[(ntemp):(ntemp -1 + length(cl[[j]]))]
          ntemp <- ntemp + length(cl[[j]])
      }
    }
    if (n.cl.ani < n.cl) {
        #for (j in ((n.cl.ani + 1):n.cl)) {   WRONG AGAIN
        for (j in (1:n.cl)[iso]) {
            theta.list[[j]] <- theta[ntemp]
            ntemp <- ntemp + 1
        }
    }
    r <- matrix(0, nrow(x), nrow(newdata))
    for (j in 1:n.cl) {
        cor.str <- covStruct.create(covtype = covtype, d = length(cl[[j]]), 
            var.names = NULL, known.covparam = "None", coef.cov = theta.list[[j]], 
            coef.var = alpha[j], iso = iso[j])
        r <- r + covMat1Mat2(object = cor.str, X1 = as.matrix(x[, 
            cl[[j]]]), X2 = as.matrix(newdata[,cl[[j]]]))
    }
    return(r)
}

### estimation of kriging parameters with block additive kernel

kmAdditive <- function(x, y, n.initial.tries = 50, 
    limits = NULL, eps.R = 1e-08, cl, covtype = "gauss", eps.Var = 1e-06, 
    max.it = 1000, iso = FALSE) {
    ### based on  LogLcl4
    n.cl <- length(cl)
    if (identical(iso, FALSE)) {
        iso <- rep(FALSE, n.cl)
    }
    if (n.cl == 1) {
        parameter <- km(~1, design = data.frame(x), response = y, covtype = covtype, 
            iso = iso)
        message("full clique: used DiceKriging::km, return object of class km")
        return(parameter)
    } else {
        n.cl.ani <- n.cl - sum(iso)
        n.cl.iso <- sum(iso)
        d <- ncol(x)
        n <- length(y)
        theta.n <- sum(as.numeric(lapply(cl[which(iso == FALSE)], length)))
        theta.n <- theta.n + sum(iso)
        if (is.null(limits))#get limits for theta
        {
          lower <- rep(1e-10, theta.n)
          upper <- numeric()
          for (i in 1:n.cl)
          {
            if (iso[i]==FALSE)
            upper <- c(upper, 2 * diff(apply(x[, cl[[i]],drop=FALSE], 2, range)))
            else upper <- c(upper, max(2*diff(apply(x[, cl[[i]],drop=FALSE], 2, range))))
          }
          limits <- data.frame(lower=lower,upper=upper)
        }
        DM <- matrix(1, ncol = 1, nrow = n)
        Ui <- c(rep(-1, n.cl - 1), rep(0, theta.n))
        Ui <- rbind(Ui, cbind(diag(rep(1, n.cl - 1)), matrix(0, ncol = theta.n, 
            nrow = n.cl - 1)))
        Ui <- rbind(Ui, cbind(matrix(0, ncol = n.cl - 1, nrow = theta.n), 
            diag(rep(1, theta.n))))
        Ui <- rbind(Ui, cbind(matrix(0, ncol = n.cl - 1, nrow = theta.n), 
            diag(rep(-1, theta.n))))
        Ci <- c(-(1 - eps.Var), rep(eps.Var, n.cl - 1), limits$lower, -limits$upper)
        ## choosing initial points
        alphastart <- matrix(runif(n.initial.tries * n.cl), ncol = n.cl)
        alphastart <- (alphastart/apply(alphastart, 1, sum))[, -n.cl]
        thetastart01 <- matrix(runif(n.initial.tries * theta.n), nrow = n.initial.tries) #  choosing a random start vector
        thetastart <- t(apply(thetastart01,1, function(a) a*(limits$upper-limits$lower)+limits$lower))
        parameterstart <- cbind(alphastart, thetastart)
    
        LLinitial <- apply(parameterstart, 1, LogLConstrained, y = y, 
            x = x, DM = DM, n = n, covtype = covtype, eps.R = eps.R, n.cl = n.cl, 
            cl = cl, iso = iso)
        parameter <- parameterstart[which.max(LLinitial), ]
        
        ### optimization the constrained LL
        test <- try(constrOptim(theta = parameter, f = LogLConstrained, 
            grad = NULL, ui = Ui, ci = Ci, y = y, x = x, DM = DM, n = n, 
            covtype = covtype, eps.R = eps.R, n.cl = n.cl, cl = cl, iso = iso, 
            outer.iterations = 10, outer.eps = 1e-05, control = list(fnscale = -1, 
                maxit = max.it, trace = 0)))
        parameter <- test$par
        parameter <- c(parameter[1:(n.cl - 1)], 1 - sum(parameter[1:(n.cl - 
            1)]), parameter[n.cl:length(parameter)])
        paramList <- paramVect2List(parameter, cl, iso)
        return(paramList)
    }
}

### kriging model prediction with block additive kernel

predictAdditive <- function(newdata, x, y, parameter, covtype = "gauss", eps.R = 1e-08, 
    cl, iso = FALSE, se.compute=FALSE) {
    n.cl <- length(cl)
    if (identical(iso, FALSE)) {
        iso <- rep(FALSE, n.cl)
    }
    if (n.cl == 1 & class(parameter) == "km") {
        message("full clique: used DiceKriging::predict.km") 
        pred <- predict(parameter, newdata = newdata, type = "UK", 
                           se.compute=se.compute, checkNames = FALSE)
        if (se.compute)
          result <- data.frame(mean=pred$mean, sd=pred$sd)
        else 
          result <- data.frame(mean=pred$mean)
    } else {
        if (n.cl == 1 & class(parameter) != "km") 
            stop("for only one clique a kriging model of class km is required for 'parameter'")
        parameter <- paramList2Vect(parameter, cl, iso)
        x <- as.matrix(x)
        alpha <- parameter[(1:n.cl)]
        theta <- parameter[-(1:n.cl)]
        n <- length(y)
        DM <- matrix(1, ncol = 1, nrow = n)
        p <- 1  #  number of trend parameters
        R <- Rfunc(as.matrix(x), theta, alpha, covtype, n, n.cl, cl, 
            iso) + diag(eps.R, ncol = n, nrow = n)
        Rinvs <- solve(R)
        beta <- solve(t(DM) %*% Rinvs %*% DM) %*% t(DM) %*% Rinvs %*% 
            y
        factor2 <- Rinvs %*% (y - DM %*% beta)
         yhat <- numeric(length(newdata[, 1]))
        r.newdata <- rfunc(newdata, x, theta, alpha, covtype, n.cl, cl, iso)
          yhat <- as.numeric(beta) + t(r.newdata)%*%factor2
    if (se.compute) {
      sigma2 <- (1/(n - p)) * t((y - DM %*% beta)) %*% 
            Rinvs %*% (y - DM %*% beta)
      sigma2 <- drop(sigma2)
  # below code from DiceKriging, one line changed because of Cor/Cov-differences (JF)   
            newdata <- as.matrix(newdata)
            C.mat <- R
            T <- chol(C.mat)
            M <- solve(t(T))%*%DM
            F.newdata <- matrix(1,nrow=nrow(newdata),1)
            c.newdata <- r.newdata
            Tinv.c.newdata <- backsolve(t(T), c.newdata, upper.tri = FALSE)
            s2.predict.1 <- apply(Tinv.c.newdata, 2, crossprod)
            T.M <- chol(t(M) %*% M)
            s2.predict.mat <- backsolve(t(T.M), t(F.newdata - 
              t(Tinv.c.newdata) %*% M), upper.tri = FALSE)
            s2.predict.2 <- apply(s2.predict.mat, 2, crossprod)
  #    line below adjusted because R and r.newdata are Cor, not Cov (JF):
            s2.predict <- sigma2*pmax(1 - s2.predict.1 + s2.predict.2,0)
            s2.predict <- as.numeric(s2.predict)
            s2.predict <- s2.predict * n/(n - p)
            sd.predict <- sqrt(s2.predict)
            result <- data.frame(mean=yhat, sd=sd.predict)
    } else result <- data.frame(mean=yhat)
    
    }
    return(result)
}
