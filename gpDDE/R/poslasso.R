###Positive-lasso###
###Based on lars.pos provided by David Reiss, further modified by Kun Chen###
###Use positive=T, type="lasso"###
###Not fully modified for positive-lar and positive-stagewise.


#Example:
#source("poslasso.R")
#X <- matrix(ncol=20,nrow=100,rnorm(2000))
#beta <- as.vector(c(rep(1,5),rep(-1,5),rep(0,10)))
#E <- as.vector(rnorm(100,0,1))
#Y <- X%*%beta + E
#fit <- lars.pos(X,Y,positive=TRUE,type="lasso",intercept=FALSE,normalize=FALSE)
#coef(fit)
#plot(fit$Cp)
#plot(fit)
###The last solution is the same as the LS-fit of the selected predictors
#lastsolution <- coef(fit)[nrow(coef(fit)),]
#as.vector(coef(lm(Y~X[,lastsolution!=0]-1)))


lars.pos <-
function (x, y, type = c("lasso", "lar", "forward.stagewise"),
    trace = FALSE, Gram, eps = .Machine$double.eps, max.steps,
    use.Gram = TRUE, positive = FALSE, normalize=TRUE,intercept=FALSE)
{
    call <- match.call()
    type <- match.arg(type)
    TYPE <- switch(type, lasso = "LASSO", lar = "LAR", forward.stagewise = "Forward Stagewise")
    if (trace)
        cat(paste(TYPE, "sequence\n"))
    nm <- dim(x)
    n <- nm[1]
    m <- nm[2]
    im <- inactive <- seq(m)
    one <- rep(1, n)
    vn <- dimnames(x)[[2]]



    if(intercept){
      meanx <- drop(one %*% x)/n
      x <- scale(x, meanx, FALSE)	# centers x
      mu <- mean(y)
      y <- drop(y - mu)
    }else {
      meanx <- rep(0,m)
      mu <- 0
      y <- drop(y)
    }


    #meanx <- (one %*% x)[1, ]/n
    #x <- scale(x, meanx, FALSE)
    #mu <- mean(y)
    #y <- my.drop(y - mu)




  if(normalize){
    normx <- sqrt(drop(one %*% (x^2)))
    nosignal<-normx/sqrt(n) < eps
    if(any(nosignal))# ignore variables with too small a variance
      {
        ignores<-im[nosignal]
        inactive<-im[-ignores]
        normx[nosignal]<-eps*sqrt(n)
        if(trace)
          cat("LARS Step 0 :\t", sum(nosignal), "Variables with Variance < eps; dropped for good\n")
          #
      }else ignores <- NULL #singularities; augmented later as well
    names(normx) <- NULL
    x <- scale(x, FALSE, normx)	# scales x
  }else {
    normx <- rep(1,m)
    ignores <- NULL
  }



    if (use.Gram & missing(Gram)) {
        if (m > 500 && n < m)
            cat("There are more than 500 variables and n<m;\nYou may wish to restart and set
            use.Gram=FALSE\n")
        if (trace)
            cat("Computing X'X .....\n")
        Gram <- t(x) %*% x
    }





    Cvec <- my.drop(t(y) %*% x)
    ssy <- sum(y^2)
    residuals <- y
    if (missing(max.steps)) max.steps <- 8 * min(m, n - intercept)
    beta <- matrix(0, max.steps + 1, m)
    Gamrat <- NULL
    arc.length <- NULL
    R2 <- 1
    RSS <- ssy
    first.in <- integer(m)
    active <- NULL
    actions <- as.list(seq(max.steps))
    drops <- FALSE
    Sign <- NULL
    R <- NULL


    k <- 0
    pos_stop <- FALSE
    while((k < max.steps) & (length(active) < min(m - length(ignores),n-intercept)) &
    pos_stop==FALSE) {
        action <- NULL
        k <- k + 1

        C <- Cvec[inactive]

        if (!positive){
          Cmax <- max(abs(C))
        }else
        {
          Cmax <- max(C)
          ##if(Cmax < -eps) pos_stop <- TRUE
        }

        if (!any(drops)) {
            if (!positive)
                new <- abs(C) >= Cmax - eps
            else new <- C >= Cmax - eps
            C <- C[!new]
            new <- inactive[new]
            for (inew in new) {
                if (use.Gram) {
                  R <- updateR.faster(Gram[inew, inew], R, my.drop(Gram[inew,
                    active]), Gram = TRUE, eps = eps)
                }else {
                  R <- updateR.faster(x[, inew], R, x[, active],
                    Gram = FALSE, eps = eps)
                }
                if (attr(R, "rank") == length(active)) {
                  nR <- seq(length(active))
                  R <- R[nR, nR, drop = FALSE]
                  attr(R, "rank") <- length(active)
                  ignores <- c(ignores, inew)
                  action <- c(action, -inew)
                  if (trace)
                    cat("LARS Step", k, ":\t Variable", inew,
                      "\tcollinear; dropped for good\n")
                }else {
                  if (first.in[inew] == 0)
                    first.in[inew] <- k
                  active <- c(active, inew)
                  if (!positive)
                    Sign <- c(Sign, sign(Cvec[inew]))
                  else Sign <- c(Sign, abs(sign(Cvec[inew])))
                  action <- c(action, inew)
                  if (trace)
                    cat("LARS Step", k, ":\t Variable", inew,
                      "\tadded\n")
                }
            }
        }
        else action <- -dropid
        Gi1 <- backsolve(R, lars::backsolvet(R, Sign))
        dropouts <- NULL
        if (type == "forward.stagewise") {
            directions <- Gi1 * Sign
            if (!all(directions > 0)) {
                if (use.Gram) {
                  nnls.object <- nnls.lars.faster(active, Sign,
                    R, directions, Gram[active, active], trace = trace,
                    use.Gram = TRUE, eps = eps)
                }
                else {
                  nnls.object <- nnls.lars.faster(active, Sign,
                    R, directions, x[, active], trace = trace,
                    use.Gram = FALSE, eps = eps)
                }
                positive <- nnls.object$positive
                dropouts <- active[-positive]
                action <- c(action, -dropouts)
                active <- nnls.object$active
                Sign <- Sign[positive]
                Gi1 <- nnls.object$beta[positive] * Sign
                R <- nnls.object$R
                C <- Cvec[-c(active, ignores)]
            }
        }
        A <- 1/sqrt(sum(Gi1 * Sign))
        w <- A * Gi1
        if (!use.Gram)
            u <- (x[, active, drop = FALSE] %*% w)[, 1]
        if (length(active) >= min(n - intercept, m - length(ignores))) {
            gamhat <- Cmax/A
        }
        else {
            if (use.Gram) {
                a <- my.drop(w %*% Gram[active, -c(active, ignores),
                  drop = FALSE])
            }
            else {
                a <- (u %*% x[, -c(active, ignores), drop = FALSE])[1,
                  ]
            }
            if (!positive)
                gam <- c((Cmax - C)/(A - a), (Cmax + C)/(A +
                  a))
            else gam <- (Cmax - C)/(A - a)
            gamhat <- min(gam[gam > eps], Cmax/A)
        }
        if (type == "lasso") {
            dropid <- NULL
            b1 <- beta[k, active]
            z1 <- -b1/w
            zmin <- min(z1[z1 > eps], gamhat)
            if (zmin < gamhat) {
                gamhat <- zmin
                drops <- z1 == zmin
            }
            else drops <- FALSE
        }
        beta[k + 1, ] <- beta[k, ]
        beta[k + 1, active] <- beta[k + 1, active] + gamhat * w
        if (use.Gram) {
            Cvec <- Cvec - gamhat * Gram[, active, drop = FALSE] %*% w
        }
        else {
            residuals <- residuals - gamhat * u
            Cvec <- (t(residuals) %*% x)[1, ]
        }
        Gamrat <- c(Gamrat, gamhat/(Cmax/A))
        arc.length <- c(arc.length, gamhat)
        if (type == "lasso" && any(drops)) {
            dropid <- seq(drops)[drops]
            for (id in rev(dropid)) {
                if (trace)
                  cat("Lasso Step", k + 1, ":\t Variable", active[id],
                    "\tdropped\n")
                R <- lars::downdateR(R, id)
            }
            dropid <- active[drops]
            beta[k + 1, dropid] <- 0
            active <- active[!drops]
            Sign <- Sign[!drops]
        }
        if (!is.null(vn))
            names(action) <- vn[abs(action)]
        actions[[k]] <- action
        inactive <- im[-c(active, ignores)]


        ###Note that for positive-lasso, some predictors may never enter the model
        ###We shall stop when all the predictors left have negative correlations with y.
        ###The last step corresponds to LS fit of the selected predictors.
        if (positive)
        {
          C <- Cvec[inactive]
          Cmax <- max(C)
          if(Cmax < -eps) pos_stop <- TRUE
        }



    }
    beta <- beta[seq(k + 1), ]
    dimnames(beta) <- list(paste(0:k), vn)
    if (trace)
        cat("Computing residuals, RSS etc .....\n")
    residuals <- y - x %*% t(beta)
    beta <- scale(beta, FALSE, normx)
    RSS <- apply(residuals^2, 2, sum)
    R2 <- 1 - RSS/RSS[1]


    ###Cp selection
    actions=actions[seq(k)]
    netdf=sapply(actions,function(x)sum(sign(x)))
    df=cumsum(netdf)### This takes into account drops
    if(intercept)df=c(Intercept=1,df+1)
    else df=c(Null=0,df)
    rss.big=rev(RSS)[1]
    df.big=n-rev(df)[1]
    if(rss.big<eps|df.big<eps)sigma2=NaN
    else
    sigma2=rss.big/df.big
    Cp <- RSS/sigma2 - n + 2 * df
    attr(Cp,"sigma2")=sigma2
    attr(Cp,"n")=n
    ##Cp <- ((n - k - 1) * RSS)/rev(RSS)[1] - n + 2 * seq(k + 1)


    object <- list(call = call, type = TYPE, R2 = R2, RSS = RSS,
        Cp = Cp, actions = actions[seq(k)], entry = first.in,
        Gamrat = Gamrat, arc.length = arc.length, Gram = if (use.Gram) Gram else NULL,
        beta = beta, mu = mu, normx = normx, meanx = meanx)
    class(object) <- c("lars", "lars.pos")
    object
}








my.drop <-
function (x)
{
    q = dim(x)
    if (length(q) == 0)
        x
    else if (q[1] == 1)
        x[1, ]
    else if (q[2] == 1)
        x[, 1]
    else x
}





nnls.lars.faster <-
function (active, Sign, R, beta, Gram, eps = 1e-10, trace = FALSE,
    use.Gram = TRUE)
{
    if (!use.Gram)
        x <- Gram
    M <- m <- length(active)
    im <- seq(m)
    positive <- im
    zero <- NULL
    while (m > 1) {
        zero.old <- c(m, zero)
        R.old <- lars::downdateR(R, m)
        beta0 <- backsolve(R.old, lars::backsolvet(R.old, Sign[-zero.old])) *
            Sign[-zero.old]
        beta.old <- c(beta0, rep(0, length(zero.old)))
        if (all(beta0 > 0))
            break
        m <- m - 1
        zero <- zero.old
        positive <- im[-zero]
        R <- R.old
        beta <- beta.old
    }
    while (TRUE) {
        while (!all(beta[positive] > 0)) {
            alpha0 <- beta.old/(beta.old - beta)
            alpha <- min(alpha0[positive][(beta <= 0)[positive]])
            beta.old <- beta.old + alpha * (beta - beta.old)
            dropouts <- match(alpha, alpha0[positive], 0)
            for (i in rev(dropouts)) R <- lars::downdateR(R, i)
            positive <- positive[-dropouts]
            zero <- im[-positive]
            beta0 <- backsolve(R, lars::backsolvet(R, Sign[positive])) *
                Sign[positive]
            beta <- beta.old * 0
            beta[positive] <- beta0
        }
        if (use.Gram) {
            w <- 1 - Sign * my.drop(Gram %*% (Sign * beta))
        }
        else {
            jw <- x %*% (Sign * beta)
            w <- 1 - Sign * my.drop(t(jw) %*% x)
        }
        if ((length(zero) == 0) || all(w[zero] <= 0))
            break
        add <- order(w)[M]
        if (use.Gram) {
            R <- updateR.faster(Gram[add, add], R, my.drop(Gram[add,
                positive]), Gram = TRUE, eps = eps)
        }
        else {
            R <- updateR.faster(x[, add], R, x[, positive], Gram = FALSE,
                eps = eps)
        }
        positive <- c(positive, add)
        zero <-  setdiff(zero, add)
        beta0 <- backsolve(R, lars::backsolvet(R, Sign[positive])) *
            Sign[positive]
        beta[positive] <- beta0
    }
    if (trace) {
        dropouts <- active[-positive]
        for (i in dropouts) {
            cat("NNLS Step:\t Variable", i, "\tdropped\n")
        }
    }
    list(active = active[positive], R = R, beta = beta, positive = positive)
}







updateR.faster <-
function (xnew, R = NULL, xold, eps = .Machine$double.eps, Gram = FALSE)
{
    xtx <- if (Gram)
        xnew
    else sum(xnew^2)
    norm.xnew <- sqrt(xtx)
    if (is.null(R)) {
        R <- matrix(norm.xnew, 1, 1)
        attr(R, "rank") <- 1
        return(R)
    }
    Xtx <- if (Gram)
        xold
    else (t(xnew) %*% xold)[1, ]
    r <- lars::backsolvet(R, Xtx)
    rpp <- norm.xnew^2 - sum(r^2)
    rank <- attr(R, "rank")
    if (rpp <= eps)
        rpp <- eps
    else {
        rpp <- sqrt(rpp)
        rank <- rank + 1
    }
    R <- cbind(rbind(R, 0), c(r, rpp))
    attr(R, "rank") <- rank
    R
}


#
# lasso.adapt.pos <-function(x,y, pars, beta, pars.c = 1){
# # adaptive lasso from lars with BIC stopping rule
# # this one uses the "known variance" version of BIC with RSS/(full model mse)
# # must use a recent version of R so that normalize=FALSE can be used in lars
# #  standardize variables like lars does
#     w.pars <- abs(pars) * pars.c
#     w.beta <- abs(beta)                      # weights for adaptive lasso
#     beta.ind <- which(w.beta > 0)
#     w.beta <- w.beta[beta.ind]
#     x <- x[, c(1:length(pars), beta.ind + length(par))]
#
#     one <- rep(1, n)
#     meanx <- drop(one %*% x)/n
#     xc <- scale(x, meanx, FALSE)         # first subtracts mean
#     normx <- sqrt(drop(one %*% (xc^2)))
#     names(normx) <- NULL
#     xs <- scale(xc, FALSE, normx)        # now rescales with norm (not sd)
#
# xs=scale(xs,center=FALSE,scale=1/w)  # xs times the weights
# object=lars(xs,y,type="lasso",normalize=FALSE)
#
# # get min BIC
# # bic=log(n)*object$df+n*log(as.vector(object$RSS)/n)   # rss/n version
# sig2f=summary(out.ls)$sigma^2        # full model mse
# bic2=log(n)*object$df+as.vector(object$RSS)/sig2f       # Cp version
# step.bic2=which.min(bic2)            # step with min BIC
#
# fit=predict.lars(object,xs,s=step.bic2,type="fit",mode="step")$fit
# coeff=predict.lars(object,xs,s=step.bic2,type="coef",mode="step")$coefficients
# coeff=coeff*w/normx                  # get back in right scale
# st=sum(coeff !=0)                    # number nonzero
# mse=sum((y-fit)^2)/(n-st-1)          # 1 for the intercept
#
# # this next line just finds the variable id of coeff. not equal 0
# if(st>0) x.ind<-as.vector(which(coeff !=0)) else x.ind<-0
# intercept=as.numeric(mean(y)-meanx%*%coeff)
# return(list(fit=fit,st=st,mse=mse,x.ind=x.ind,coeff=coeff,intercept=intercept,object=object,
#             bic2=bic2,step.bic2=step.bic2))
# }
