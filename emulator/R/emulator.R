"betahat.fun" <-
function (xold, Ainv, d, give.variance = FALSE, func = regressor.basis) 
{
    H <- regressor.multi(xold, func = func)
    H <- as.matrix(H)
    out <- solve(quad.form(Ainv, H), cprod(cprod(Ainv, 
        H), d))
    out <- as.vector(out)
    names(out) <- colnames(H)
    if (give.variance) {
        jj.sigma <- drop(sigmahatsquared(H, Ainv, d))
        return(list(betahat = out, sigmahatsquared = jj.sigma, 
            variance = jj.sigma * solve(quad.form(Ainv, H))))
    }
    else {
        return(out)
    }
}
"betahat.fun.A" <-
function (xold, A, d, give.variance = FALSE, func = regressor.basis) 
{
    H <- regressor.multi(xold, func = func)
    H <- as.matrix(H)
    colnames(H)[1] <- "constant"
    out <- solve(quad.form.inv(A, H), cprod(H, solve(A, d)))
    out <- as.vector(out)
    names(out) <- rownames(H)
    if (give.variance) {
        jj.sigma <- drop(sigmahatsquared.A(H, A, d))
        return(list(betahat = out, sigmahatsquared = jj.sigma, 
            variance = jj.sigma * solve(quad.form.inv(A, H))))
    }
    else {
        return(out)
    }
}
"corr" <-
function (x1, x2, scales = NULL, pos.def.matrix = NULL, 
    coords = "cartesian", spherical.distance.function = NULL) 
{
    if (is.null(scales) & is.null(pos.def.matrix)) {
        stop("need either scales or a pos.definite.matrix")
    }
    if (!is.null(scales) & !is.null(pos.def.matrix)) {
        stop("scales *and* pos.def.matrix supplied.  corr() needs one only.")
    }
    if (is.null(pos.def.matrix)) {
        pos.def.matrix <- diag(scales, nrow = length(scales))
    }
    x1 <- as.vector(unlist(x1))
    x2 <- as.vector(unlist(x2))
        corrscale <- function(x) {
            exp(-abs(x))
        }
    m <- switch(coords,
                cartesian = x1 - x2,
                spherical = spherical.distance.function(x1, x2),
                stop("coords must be either Cartesian or Spherical")
                )
    return(corrscale(quad.form(pos.def.matrix, m)))
}
"corr.matrix" <-
function (xold, yold = NULL, method = 1, distance.function = corr, 
    ...) 
{
    if (is.null(yold)) {
      nully <- TRUE
        yold <- xold
    } else {
      nully <- FALSE
    }

    if(!identical(distance.function, corr) & identical(method,1)){
      method <- 2
    }

    if(identical(method,1)) {
      a <- list(...)
      scales <- a$scales
      if(length(scales)==1){
        scales <- rep(scales,ncol(xold))
      }
      
      pos.def.matrix <- a$pos.def.matrix
      if (is.null(scales) & is.null(pos.def.matrix)) {
        stop("need either scales or a pos.definite.matrix")
      }
      if (!is.null(scales) & !is.null(pos.def.matrix)) {
        stop("scales *and* pos.def.matrix supplied.  corr() needs one only.")
      }
      if (is.null(pos.def.matrix)) {
        pos.def.matrix <- diag(scales, nrow = length(scales))
      }
      jj <- function(x){as.matrix(diag(as.matrix(x)))}
      if(nully){
        R <- quad.tform(pos.def.matrix, xold)
        S <- kronecker(jj(R),t(rep(1,nrow(xold))))
        return(exp(Re(R+t(R)-S-t(S))))  #sic
      }
      txold <- ht(xold)
      tyold <- ht(yold)
      R1 <- quad.form(pos.def.matrix,txold)
      R2 <- quad.form(pos.def.matrix,tyold)
      

      S1 <- kronecker(t(jj(R1)),rep(1,nrow(yold)))
      S2 <- kronecker(jj(R2),t(rep(1,nrow(xold))))

      a1 <- cprod(txold, pos.def.matrix) %*% tyold
      a2 <- cprod(tyold, pos.def.matrix) %*% txold
      return(exp(Re(t(a1)+a2-S1-S2)))
    } else if (identical(method,2)){
      out <- apply(xold, 1, function(y) {
        apply(yold, 1, function(x) {
          distance.function(x, y, ...)
        })
      })
      return(as.matrix(out))
    }
    else if (identical(method,3)){
      n <- nrow(xold)
      m <- nrow(yold)
      A <- matrix(NA, m, n)
      for (i in 1:n) {
        for (j in 1:m) {
          A[j, i] <- distance.function(xold[i, ], yold[j, 
                                                       ], ...)
        }
      }
      colnames(A) <- rownames(xold)
      rownames(A) <- rownames(yold)
      return(as.matrix(A))
    } else {
      stop("method must be 1, 2, or 3")
    }
}

"estimator" <-
function (val, A, d, scales = NULL, pos.def.matrix = NULL, func=regressor.basis) 
{
    d.missing.estimated <- d + NA
    for (i in 1:nrow(val)) {
        val.oneshort <- val[-i, ,drop=FALSE]
        val.missing <- val[i, ]
        d.oneshort <- d[-i]
        d.missing <- d[i]
        A.oneshort <- A[-i, -i]
        Ainv.oneshort <- solve(A.oneshort)
        d.missing.estimated[i] <- interpolant(val.missing, d.oneshort, 
            val.oneshort, Ainv = Ainv.oneshort, scales = scales, 
            pos.def.matrix = pos.def.matrix, func=func, give.full.list = TRUE)$mstar.star
    }
    return(d.missing.estimated)
}

"interpolant" <-
function (x, d, xold, Ainv = NULL, A = NULL, use.Ainv = TRUE, 
          scales = NULL, pos.def.matrix = NULL, func = regressor.basis,
          give.full.list = FALSE, distance.function=corr, ...) 
{
    if (is.null(scales) & is.null(pos.def.matrix)) {
        stop("need either scales or a pos.definite.matrix (used to calculate tx)")
    }
    if (!is.null(scales) & !is.null(pos.def.matrix)) {
        stop("scales *and* pos.def.matrix supplied.  corr() needs one only.")
    }
    if (is.null(pos.def.matrix)) {
        pos.def.matrix <- diag(scales,nrow=length(scales))
    }
    if (is.null(A)) {
        A <- corr.matrix(xold=xold, pos.def.matrix = pos.def.matrix,
                         distance.function=distance.function, ...)
    }
    if (is.null(Ainv) & use.Ainv) {
        Ainv <- solve(A)
    }
#    tx <- apply(xold, 1, distance.function, x2 = x, pos.def.matrix = pos.def.matrix, ...)
    tx <- Conj(drop(corr.matrix(xold=xold, yold=t(x),
                      distance.function=distance.function,
                      pos.def.matrix = pos.def.matrix, ...)))
    hx <- unlist(func(x))
    H <- regressor.multi(xold, func = func)
    if (use.Ainv) {
        n <- nrow(Ainv)
        jj <- betahat.fun(xold, Ainv, d, give.variance = TRUE, 
            func = func)
        betahat <- jj$betahat
        beta.var <- jj$variance
        sigmahat.square <- jj$sigmahatsquared
        beta.marginal.sd <- sqrt(diag(beta.var))
        prior <- cprod(hx,betahat)
        mstar.star <- prior + cprod(cprod(Ainv, 
            tx), d - H %*% betahat)
        cstar.x.x <- corr(x, x, pos.def.matrix = pos.def.matrix, ...) - quad.form(Ainv, tx)
        cstar.star <- cstar.x.x + quad.form.inv(quad.form(Ainv, 
            H), hx - Conj(cprod(H, cprod(Ainv, tx))))
    }
    else {
        n <- nrow(A)
        jj <- betahat.fun.A(xold, A, d, give.variance = TRUE, 
            func = func)
        betahat <- jj$betahat
        sigmahat.square <- jj$sigmahatsquared
        beta.var <- jj$variance
        beta.marginal.sd <- sqrt(diag(beta.var))
        prior <- drop(cprod(hx,betahat))
        mstar.star <- drop(prior + cprod(solve(A, 
            tx), d - H %*% betahat))
        cstar.x.x <- corr(x, x, pos.def.matrix = pos.def.matrix, 
            ...) - quad.form.inv(A, tx)
        cstar.star <- cstar.x.x + quad.form.inv(quad.form.inv(A, 
            H), hx - Conj(cprod(H, solve(A, tx))))
    }
    if (give.full.list) {
        return(list(betahat = betahat, prior = prior,
                    beta.var = beta.var, beta.marginal.sd = beta.marginal.sd, 
            sigmahat.square = sigmahat.square, mstar.star = mstar.star, 
            cstar = cstar.x.x, cstar.star = cstar.star, Z = sqrt(abs(sigmahat.square * 
                cstar.star))))
    }
    else {
        return(as.vector(mstar.star))
    }
}

"int.qq" <- function(x, d, xold, Ainv, pos.def.matrix, func=regressor.basis){
  bhat <- betahat.fun(xold,Ainv,d,func=func)
  out <- 
    cprod(apply(x,1,func),bhat) + 
      cprod(
                cprod(Ainv,corr.matrix(x,xold,pos.def.matrix=pos.def.matrix)),
                d-cprod(apply(xold,1,func), bhat))
return(as.vector(out))
}

"interpolant.quick" <-
function (x, d, xold, Ainv=NULL, scales = NULL, pos.def.matrix = NULL, 
    func = regressor.basis, give.Z = FALSE, distance.function=corr, ...) 
{
    if (is.null(scales) & is.null(pos.def.matrix)) {
        stop("need either scales or a pos.definite.matrix")
    }
    if (!is.null(scales) & !is.null(pos.def.matrix)) {
        stop("scales *and* pos.def.matrix supplied.  corr() needs one only.")
    }
    if (is.null(pos.def.matrix)) {
        pos.def.matrix <- diag(scales,nrow=length(scales))
    }
    if (is.null(Ainv)){
      Ainv <- solve(corr.matrix(xold=xold,pos.def.matrix=pos.def.matrix))
    }
      
    
    betahat <- betahat.fun(xold, Ainv, d, func = func)
    H <- regressor.multi(xold, func = func)
    mstar.star <- rep(NA, nrow(x))
    prior <- rep(NA,nrow(x))

    if (give.Z) {
        Z <- rep(NA, nrow(x))
        sigmahat.square <- sigmahatsquared(H, Ainv, d)
    }
    for (i in 1:nrow(x)) {
        hx <- func(x[i, ])
        tx <- as.vector(corr.matrix(yold=xold,
                                    xold = x[i,,drop=FALSE],
                                    method=1,
                                    pos.def.matrix = pos.def.matrix, 
                                    distance.function=distance.function,
                                    ...
                                    )
                        )
        prior[i] <- cprod(hx,betahat)

        mstar.star[i] <- prior[i] + cprod(cprod(Ainv, 
            tx), (d - H %*% betahat))
        if (give.Z) {
            cstar.x.x <- 1 - quad.form(Ainv, tx)
            cstar.star <- cstar.x.x + quad.form.inv(quad.form(Ainv, 
                H), hx - Conj(cprod(H, cprod(Ainv, tx))))
            Z[i] <- sqrt(abs(sigmahat.square * cstar.star))
        }
    }
    if (give.Z) {
        return(list(mstar.star = mstar.star, Z = Z, prior = prior))
    }
    else {
        return(mstar.star)
    }
}

"var.conditional" <-
function (x, xold, d, A, Ainv, scales = NULL, pos.def.matrix = NULL, 
    func = regressor.basis, distance.function=corr, ...) {
    if (is.null(scales) & is.null(pos.def.matrix)) {
      stop("need either scales or a pos.definite.matrix")
    }
    if (!is.null(scales) & !is.null(pos.def.matrix)) {
      stop("scales *and* pos.def.matrix supplied.  corr() needs one only.")
    }
    if (is.null(pos.def.matrix)) {
      pos.def.matrix <- diag(scales,nrow=length(scales))
    }
    
    H <- regressor.multi(xold, func = func)
    hx <- regressor.multi(x,func=func)
#    txvec <- 
#      apply(x,1,function(y){apply(xold,1,function(x){distance.function(x,y,pos.def.matrix=pos.def.matrix, ...)})})

    txvec <-
    corr.matrix(xold=x,yold=xold,distance.function=distance.function,
    pos.def.matrix=pos.def.matrix, ...)
    bit1 <- quad.form(Ainv,txvec)
    jj <- hx - cprod(txvec,Ainv) %*% H
    bit2 <- quad.form.inv(quad.form(Ainv,H),t(jj))
    
    cstar <- corr.matrix(xold=x,pos.def.matrix=pos.def.matrix,
    distance.function=distance.function, ...) - bit1 + bit2
    cstar <- as.matrix(cstar) 
    
    rownames(cstar) <- rownames(x)
    colnames(cstar) <- rownames(x)
    
    return(
           cstar*
           sigmahatsquared(H=regressor.multi(xold, func = func),
                           Ainv=Ainv, d=d)
           )
}

"cond.sample" <- function (n=1, x, xold, d, A, Ainv, scales = NULL, pos.def.matrix = NULL, 
    func = regressor.basis, ...) {
  mstar <- 
    interpolant.quick(x=x, d=d, xold=xold, Ainv=Ainv, scales=scales,
                      pos.def.matrix=pos.def.matrix, func=func,
                      give.Z=FALSE)

jj.sigma <- var.conditional(x, xold, d, A, Ainv, scales = scales, pos.def.matrix = pos.def.matrix, 
               func = func, ...)

  random.bit <- 
    rmvt(n=n,sigma=jj.sigma, df=length(d))

  out <- sweep(random.bit,2,mstar,"+")
  colnames(out) <- rownames(x)
  return(out)
}

"latin.hypercube" <-
function (n, d, names=NULL, normalize = FALSE, complex=FALSE) 
{
    if(complex){
        out <- Recall(n,2*d, normalize=normalize, complex=FALSE)
        out <- out[,seq_len(d)] + 1i*out[,-seq_len(d)]
        colnames(out) <- names
        return(out)
    }

        
    if (normalize) {
        f <- function(...) {
            sample(0:(n - 1))/(n - 1)
        }
    }
    else {
        f <- function(...) {
            (sample(1:n) - 0.5)/n
        }
    }
    out <- sapply(1:d, f)
    colnames(out) <- names
    return(out)
}

"makeinputfiles" <-
function (number.of.runs = 100, gaussian = TRUE, directoryname = "~/goldstein/genie-cgoldstein/", 
    filename = "QWERTYgoin", expert.estimates, area.outside = 0.05) 
{
    pad <- function(y, len, padchar = "0") {
        paste(paste(rep(padchar, len - nchar(y)), collapse = ""), 
            y, collapse = "", sep = "")
    }
    a <- expert.estimates
    minimum.values <- a$low
    maximum.values <- a$high
    lh.normalized.uniformdist <- latin.hypercube(number.of.runs, 
        nrow(a))
    lh.real <- lh.normalized.uniformdist
    lh.real[] <- NA
    if (gaussian) {
        for (i in 1:ncol(lh.real)) {
            a <- minimum.values[i]
            b <- maximum.values[i]
            if (a > 0) {
                lh.real[, i] <- exp(qnorm(p = lh.normalized.uniformdist[, 
                  i], mean = (log(a) + log(b))/2, sd = (log(b) - 
                  log(a))/(2 * qnorm(1 - area.outside/2))))
            }
            else {
                lh.real[, i] <- qunif(p = lh.normalized.uniformdist[, 
                  i], min = 0, max = b)
            }
        }
        lh.real[is.nan(lh.real)] <- 0
    }
    else {
        lh.normalized <- lh.normalized.uniformdist
        unnormalize <- function(normalized.vector) {
            minimum.values + normalized.vector * (maximum.values - 
                minimum.values)
        }
        lh.real <- t(apply(lh.normalized, 1, unnormalize))
    }
    write.table(x = lh.real, file = paste(directoryname, "list_of_inputfiles.txt", 
        sep = ""), quote = FALSE, row.names = TRUE, col.names = FALSE)
    for (i in 1:number.of.runs) {
        fullfilename <- paste(directoryname, filename, ".", i - 
            1, sep = "")
        f <- function(string, append = TRUE) {
            write(string, file = fullfilename, append = append)
        }
        f("400001 400000 400000  4000  400", append = FALSE)
        f("n")
        f("100   5")
        f("20.0")
        f("20.0")
        f("0.90")
        for (j in 1:8) {
            f(lh.real[i, j])
        }
        f("0.0")
        f("0.0")
        for (j in 9:11) {
            f(lh.real[i, j])
        }
        f("0.0")
        for (j in 12:12) {
            f(lh.real[i, j])
        }
        f("1")
        for (j in 13:14) {
            f(lh.real[i, j])
        }
        f("0.0")
        f("0.0")
        f("0.0")
        for (j in 16:17) {
            f(lh.real[i, j])
        }
        f(paste("tmp/", pad(i - 1, 3), sep = ""))
        f("tmp/tmp.avg")
    }
    return(0)
}

"model" <-
function (x) 
{
    if (1 == 1) {
        centrepoint <- x
        x[] <- 4
        return(sum((x - centrepoint)^2))
    }
    else {
        return(as.double(x[1] < x[2]))
    }
}

"optimal.scales" <-
function (val, scales.start, d, use.like = TRUE, give.answers = FALSE, func=regressor.basis,
    ...) 
{
    if (use.like) {
        objective.fun <- function(scales, val, d) {
            -scales.likelihood(scales = exp(scales), xold = val, 
                d = d, give_log=TRUE, func=func)
        }
    }
    else {
        jj <- function(scales, val, d) {
            A <- corr.matrix(xold=val, scales = exp(scales))
            error <- abs(d - estimator(val, A, d, scales = exp(scales)))
            return(sum(error^2))
        }
        initial.value <- jj(scales=scales.start, val=val, d=d)
        objective.fun <- function(scales,val,d){jj(scales,val,d)/initial.value}
        
    }
    jj <- optim(par = log(scales.start), objective.fun, val = val, 
        d = d, ...)
    if (give.answers) {
        return(jj)
    }
    else {
        return(exp(jj$par))
    }
}

"optimal.scale" <-
function (val, d, use.like = TRUE, give.answers = FALSE,  func=regressor.basis,
    ...) 
{
    n <- ncol(val)
    if (use.like) {
        objective.fun <- function(scale, val, d) {
            out <- -scales.likelihood(scales = rep(exp(scale),n), xold = val, 
                d = d, give_log=TRUE, func=func)
        }
    }
    else {
        objective.fun <- function(scale, val, d) {
            A <- corr.matrix(xold=val, scales = rep(exp(scale),n))
            error <- abs(d - estimator(val, A, d, scales = rep(exp(scale),n), func=func))
            return(sum(error^2))
        }
    }
    jj <- optimize(f=objective.fun, interval=c(-4,10), maximum=FALSE, val = val, d = d, ...)
    if (give.answers) {
        return(jj)
    }
    else {
      return(exp(jj$minimum))
    }
}

"pad" <-
function (x, len, padchar = "0", strict = TRUE) 
{
    n <- nchar(x)
    if (nchar(padchar) > 1) {
        stop("padchar must be a single character")
    }
    if (n > len) {
        if (strict) {
            stop("input arg too long")
        }
        else {
            return(substr(as.character(x), n - len + 1, n))
        }
    }
    return(paste(paste(rep(padchar, len - n), collapse = ""), 
        x, sep = "", collapse = ""))
}

"prior.B" <-
function (H, Ainv, B0 = NULL) 
{
    if (is.null(B0)) {
        return(quad.form(Ainv, H))
    }
    else {
        return(B0 + quad.form(Ainv, H))
    }
}

"prior.b" <-
function (H, Ainv, d, b0 = NULL, B0 = NULL) 
{
    B <- prior.B(H, Ainv, B0)
    if (is.null(b0)) {
        b <- cprod(solve(B), cprod(H, Ainv)) %*% d
    }
    else {
        b <- solve(B) %*% (B0 %*% b0 + cprod(H, Ainv) %*% 
            d)
    }
    return(b)
}

"tr" <-
function (a) 
{
    i <- 1:nrow(a)
    return(sum(a[cbind(i, i)]))
}  

"regressor.basis" <-
function (x) 
{
    x <- c(1, x)
    names(x)[1] <- "const"
    return(x)
}

"regressor.multi" <-
function (x.df, func = regressor.basis) 
{
  out <- t(apply(x.df, 1, func))
  if(nrow(out) == 1){
    return(t(out))
  } else {
    return(out)
  }
}

"s.chi" <-
function (H, Ainv, d, s0 = 0, fast.but.opaque = TRUE) 
{
    if (fast.but.opaque) {
        out <- s0 + quad.form(Ainv - quad.form(quad.form(solve(quad.form(Ainv, 
            H)), t(H)), Ainv), d)
    }
    else {
        out <- s0 + t(d) %*% (Ainv - Ainv %*% H %*% solve(t(H) %*% 
            Ainv %*% H) %*% t(H) %*% Ainv) %*% d
    }
    return(out)
}

"sample.from.exp.est" <-
function (number.of.runs, expert.estimates, gaussian = TRUE, 
    area.outside = 0.05) 
{
    a <- expert.estimates
    minimum.values <- a$low
    maximum.values <- a$high
    lh.normalized.uniformdist <- latin.hypercube(number.of.runs, 
        16)
    lh.real <- lh.normalized.uniformdist
    lh.real[] <- NA
    if (gaussian) {
        for (i in 1:ncol(lh.real)) {
            a <- minimum.values[i]
            b <- maximum.values[i]
            lh.real[, i] <- exp(qnorm(p = lh.normalized.uniformdist[, 
                i], mean = (log(a) + log(b))/2, sd = (log(b) - 
                log(a))/(2 * qnorm(1 - area.outside/2))))
        }
        lh.real[is.nan(lh.real)] <- 0
    }
    else {
        lh.normalized <- lh.normalized.uniformdist
        unnormalize <- function(normalized.vector) {
            minimum.values + normalized.vector * (maximum.values - 
                minimum.values)
        }
        lh.real <- t(apply(lh.normalized, 1, unnormalize))
    }
    return(lh.real)
}

"sample.n.fit" <-
  function(n=10 , scales.generate=100 , scales.fit =
           100, func=regressor.basis, ...
           )
{
  
  if(is.null(func)){
    func <-
      function(x){out <- c(1,x)
                  names(out)[1] <- "const"
                  return(out)
                }
  }   
  toy <- as.matrix(seq(from=0,to=1,len=n))
  colnames(toy) <- "a"
  rownames(toy) <- paste("obs",1:nrow(toy),sep=".")
  
  x <- seq(from=-1,to=2,len=200)
  A <- corr.matrix(xold=toy,scales=scales.fit)
  Ainv <- solve(A)
  
  d.noisy <- as.vector(rmvnorm(n=1 , mean=toy*0 ,
                               sigma=corr.matrix(xold=toy,scales=scales.generate)
                               ))
  
  jj <- interpolant.quick(as.matrix(x), d.noisy, toy, scales=scales.fit,
                          func=func, 
                          Ainv=Ainv, give.Z=TRUE)
  
  plot(x,jj$mstar.star,xlim=range(x),type="l",col="black",lwd=3, ...)
  lines(x,jj$prior,col="green",type="l")
  lines(x,jj$mstar.star+jj$Z,type="l",col="red",lty=2)
  lines(x,jj$mstar.star-jj$Z,type="l",col="red",lty=2)
  points(toy,d.noisy,pch=16,cex=2)
  legend("topright",lty=c(1:2,0),col=c("black","red","green","black"),pch=c(NA,NA,NA,16),legend=c("best estimate","+/- 1sd","prior","training set"))
}

"scales.likelihood" <-
function (pos.def.matrix = NULL, scales = NULL,  xold,
          use.Ainv = TRUE, d, give_log=TRUE, func = regressor.basis) 
{
    if (is.null(scales) & is.null(pos.def.matrix)) {
        stop("need either scales or a pos.definite.matrix")
    }
    if (!is.null(scales) & !is.null(pos.def.matrix)) {
        stop("scales *and* pos.def.matrix supplied.  corr() needs one only.")
    }
    if (is.null(pos.def.matrix)) {
        pos.def.matrix <- diag(scales,nrow=length(scales))
    }
    H <- regressor.multi(xold, func = func)
    q <- ncol(H) 
    n <- nrow(H)
    A <- corr.matrix(xold=xold, pos.def.matrix = pos.def.matrix)

    ## Define a function that returns log(1/sqrt(det(x)))
    f <- function(M){(-0.5) * sum(log(eigen(M,TRUE,TRUE)$values))}
    ## Also note that f() traps any nonpositive eigenvalue [ie non-positive definite M]
    
    bit2 <- f(A)

    if(use.Ainv){
      Ainv <- solve(A)
      bit1 <- log(sigmahatsquared(H, Ainv, d)) * (-(n - q)/2)
      bit3 <- f(quad.form(Ainv, H))
    } else {
      bit1 <- log(sigmahatsquared.A(H, A, d))  * (-(n - q)/2)
      bit3 <- f(quad.form.inv(A, H))
    }
    out <- drop(bit1 + bit2 + bit3)
    if(give_log){
      return(out)
    } else {
      return(exp(out))
    }
}

"sigmahatsquared" <-
function (H, Ainv, d) 
{
    n <- nrow(Ainv)
    q <- ncol(H) 
    H <- as.matrix(H)
    out <- quad.form(Ainv - quad.form.inv(quad.form(Ainv, H), 
        cprod(H, Ainv)), d)/(n - q - 2)
    return(Re(out))
}
"sigmahatsquared.A" <-
function (H, A, d) 
{
    n <- nrow(A)
    q <- ncol(H)
    (quad.form.inv(A, d) - quad.form(quad.form.inv(quad.form.inv(A, 
        H), t(solve(A, H))), d))/(n - q - 2)
}
