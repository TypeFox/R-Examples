require(emulator)

"Afun" <-
function (level, Di, Dj, hpa) 
{
    corr.matrix(Di, Dj, pos.def.matrix = hpa$B[[level]])
}

"H.fun.app" <-
function (D1, subsets, basis, hpa) 
{
   p.matrix <-
     function (rho) 
       {
         n <- length(rho)
         out <- diag(1, nrow = n + 1)
         out[cbind((1:n) + 1, 1:n)] <- -rho
         solve(out)
       }
    L <- as.sublist(D1, subsets)
    X <- p.matrix(hpa$rhos)
    f <- function(j) {
        out <- kronecker(X[j, , drop = FALSE], basis(L[[j]]))
        rownames(out) <- rownames(L[[j]])
        return(out)
    }
    tmp <- sapply(1:length(L), f)
    out <- do.call("rbind", tmp)
    jj <- paste("level", rep(1:length(L), each = ncol(basis(L[[1]]))), 
        sep = "")
    jj <- paste(jj, colnames(basis(L[[1]])), sep = ".")
    colnames(out) <- jj
    return(out)
}
"Pi" <-
function (hpa, i, j) 
{
    if (i > j + 1) {
        return(0)
    }
    else if (i == j + 1) {
        return(1)
    }
    else {
        return(prod(hpa$rhos[seq(from = i, to = j)]))
    }
}

"V.fun.app" <-
function (D1, subsets, hpa) 
{
    n <- sum(as.numeric(lapply(subsets, length)))
    out <- matrix(NA, n, n)
    f <- function(k, l) {
        Dk <- D1[subsets[[k]], , drop = FALSE]
        Dl <- D1[subsets[[l]], , drop = FALSE]
        out <- matrix(0, nrow(Dl), nrow(Dk))
        for (i in 1:min(k, l)) {
            out <- out + Pi(hpa, i, k - 1) * Pi(hpa, i, l - 1) * 
                hpa$sigma_squareds[i] * Afun(level = i, Dk, Dl, hpa)
        }
        return(out)
    }
    jj <- 1:length(subsets)
    out <- do.call("cbind", lapply(jj, function(x) {
        do.call("rbind", lapply(jj, f, k = x))
    }))
    rownames(out) <- c(lapply(subsets,names),recursive=TRUE)
    colnames(out) <- rownames(out)
    return(out)
}

"as.sublist" <-
function (D1, subsets) 
{
    out <- NULL
    for (i in 1:length(subsets)) {
        out[[i]] <- D1[subsets[[i]], , drop = FALSE]
        rownames(out[[i]]) <- names(subsets[[i]])
    }
    return(out)
}

"basis.toy" <-
function (x) 
{
    if (is.vector(x)) {
        stopifnot(length(x) == 3)
        out <- c(1, x, x[1]^2)
        names(out) <- c("const", LETTERS[1:3], "quad")
        return(out)
    }
    else {
        return(t(apply(x, 1, match.fun(sys.call()[[1]]))))
    }
}



"betahat.app" <-
function(D1, subsets, basis, hpa, z, use.Vinv=TRUE)
{
  H <- H.fun.app(D1=D1,subsets=subsets,basis=basis,hpa=hpa)
  V <- V.fun.app(D1=D1,subsets=subsets,hpa=hpa)
  if(use.Vinv){
    return(betahat.app.H(H=H,Vinv=solve(V),z=z))
  } else {
    return(betahat.app.H(H=H,V=V,z=z))
  }
}

"betahat.app.H" <-
function (H, V = NULL, Vinv = NULL, z) 
{
    z <- unlist(z)
    if (!is.null(Vinv)) {
      return(drop(crossprod(solve(quad.form(Vinv,H)),crossprod(H, crossprod(Vinv, z)))))
    } else {
      return(drop(solve(quad.form.inv(V, H), crossprod(H, solve(V, z)))))
    }
}

"generate.toy.observations" <-
function (D1, subsets, basis.fun, hpa = NULL, betas = NULL, export.truth=FALSE)
{

    if(is.null(hpa)){
      hpa <-
        hpa.fun.toy(
                    c(
          c(sigma_squared1 = 0.4, sigma_squared2 = 0.5, sigma_squared3 = 0.2, sigma_squared4 = 0.1), 
          scales.level1 = c(A = 1.1, B = 1.5, C = 1.9),
          scales.level2 = c(A = 1.6, B = 2.2, C = 2.9),
          scales.level3 = c(A = 0.4, B = 0.9, C = 1.2),
          scales.level4 = c(A = 0.7, B = 0.5, C = 1.6), 
          rhos = c(level1 = 1.1, level2 = 1.2, level3 = 1.3)
                      )
                    )
    }
      
      if(is.null(betas)){
        betas <-
          rbind(c(1, 2, 3, 4, 5),
                c(1, 1, 3, 4, 5),
                c(1, 1, 1, 4, 5),
                c(1, 1, 1, 1, 5))
        colnames(betas) <- c("const", LETTERS[1:3], "quad")
        rownames(betas) <- paste("level", 1:4, sep = "")
      }
      
      if(export.truth){
        return(list(
                    hpa=hpa,
                    betas=betas
                    )
               )
      }


    sigma_squareds <- hpa$sigma_squareds
    B <- hpa$B
    rhos <- hpa$rhos
    delta <- function(i) {
      out <- rmvnorm(n = 1,
                     mean = basis.fun(D1[subsets[[i]], , drop =
                       FALSE]) %*% betas[i, ],
                     sigma = sigma_squareds[i] * corr.matrix(xold = D1[subsets[[i]], , drop = FALSE], pos.def.matrix = B[[i]])
                     )
      out <- drop(out)
      names(out) <- rownames(D1[subsets[[i]], , drop = FALSE])
      return(out)
    }
    
    use.clever.but.untested.method <- FALSE
    
    if(use.clever.but.untested.method){
      z1 <- delta(1)
      z2 <- delta(2) + rhos[1] * z1[match(subsets[[2]], subsets[[1]])]
      z3 <- delta(3) + rhos[2] * z2[match(subsets[[3]], subsets[[2]])]
      z4 <- delta(4) + rhos[3] * z3[match(subsets[[4]], subsets[[3]])]
      return(list(z1 = z1, z2 = z2, z3 = z3, z4 = z4))
    } else {
      out <- NULL
      out[[1]] <- delta(1)
      for(i in 2:length(subsets)){
        out[[i]] <- delta(i) + rhos[i-1] *
    out[[i-1]][match(subsets[[i]], subsets[[i-1]])]
      }
      return(out)
    }
  }

"hdash.fun" <-
function (x, hpa, basis) 
{
  rhos <- hpa$rhos
  jj <- prod(rhos)/c(1, cumprod(rhos))
  names(jj) <- names(hpa$sigma_squareds)
  if (is.vector(x)) {
    x <- t(x)
  }
  out <- kronecker(jj, t(basis(x)), make.dimnames = TRUE)
  colnames(out) <- rownames(x)
  return(drop(out))
}

"hpa.fun.toy" <-
function (x) 
{
    if (length(x) != 19) {
        stop("x must have 19 elements")
    }
    "pdm.maker" <- function(x) {
        jj <- diag(x[1:3])
        rownames(jj) <- LETTERS[1:3]
        colnames(jj) <- LETTERS[1:3]
        return(jj)
    }
    sigma_squareds <- x[1:4]
    names(sigma_squareds) <- paste("level", 1:4, sep = "")
    B <- list()
    B[[1]] <- pdm.maker(x[5:7])
    B[[2]] <- pdm.maker(x[8:10])
    B[[3]] <- pdm.maker(x[11:13])
    B[[4]] <- pdm.maker(x[14:16])
    rhos <- x[17:19]
    names(rhos) <- paste("level", 1:3, sep = "")
    return(list(sigma_squareds = sigma_squareds, B = B, rhos = rhos))
}

"is.consistent" <-
function (subsets, z) 
{
    all(unlist(lapply(subsets, length)) == unlist(lapply(z, length)))
}

"is.nested" <-
function (subsets) 
{
    !any(sapply(mapply(setdiff, subsets[-1], subsets[-length(subsets)]), 
        length))
}

"is.strict" <-
function (subsets) 
{
    all(diff(unlist(lapply(subsets, length))) < 0)
}

"mdash.fun" <-
function (x, D1, subsets, hpa, Vinv = NULL, use.Vinv = TRUE, 
    z, basis) 
{
    H <- H.fun.app(D1 = D1, subsets = subsets, basis = basis, hpa = hpa)
    hdash <- hdash.fun(x, hpa, basis)
    tx <- tee.fun(x, D1, subsets, hpa)
    if (use.Vinv) {
        if (is.null(Vinv)) {
            Vinv <- solve(V.fun.app(D1, subsets, hpa))
        }
        bhat <- betahat.app.H(H = H, Vinv = Vinv, z = z)
        out <- drop(crossprod(hdash, bhat) + crossprod(tx , Vinv) %*% (unlist(z) - 
            H %*% bhat))
        names(out) <- rownames(x)
        return(out)
    }
    else {
        V <- V.fun.app(D1, subsets, hpa)
        bhat <- betahat.app.H(H = H, V = V, z = z)
        out <- drop(crossprod(hdash,bhat) + crossprod(tx , solve(V, unlist(z) - H %*% 
            bhat)))
        names(out) <- rownames(x)
        return(out)
    }
}

"object" <-
function(level,D,z,basis,subsets,hpa)
{
  rhos <- hpa$rhos
  sigma_squared <- hpa$sigma_squareds
  if(FALSE){
    bit1 <- log(det(Afun(level=level,Di=D[subsets[[level]],,drop=FALSE],Dj=NULL,hpa=hpa)))
  }
  jj <- Afun(level=level,Di=D[subsets[[level]],,drop=FALSE],Dj=NULL,hpa=hpa)
  
  bit1 <- sum(log(abs(eigen(jj,TRUE,TRUE)$values)))

  n.i <- length(z[[level]])
  bit2 <- n.i*log(sigma_squared[[level]])


  if(level == 1){
    d <- z[[1]]
  } else {
    d <- z[[level]] -
      rhos[level-1]*z[[level-1]][match(subsets[[level]],subsets[[level-1]])]
  }
  jj <- betahat.app(D1=D,subsets=subsets,basis=basis,hpa=hpa,z=z)
  u <- length(jj)/length(z)
  betahat <- jj[ (1+(level-1)*u):(level*u)]
  mismatch <- d -  basis(D[subsets[[level]],]) %*% betahat
  mat <- solve(sigma_squared[[level]]*Afun(level=level,Di=D[subsets[[level]],,drop=FALSE],Dj=NULL,hpa=hpa))
  bit3 <-  quad.form(mat,mismatch)
  return(bit1 + bit2 + bit3)
}
  
"subsets.fun" <-
function (n, levels = 4, prob = 0.7) 
{
    "select" <-
    function (a, prob) 
    {
        a[c(FALSE, TRUE)[1 + rbinom(n = length(a), size = 1, prob = prob)]]
    }
    a <- 1:n
    out <- NULL
    for (i in 1:levels) {
        out[[i]] <- a
        jj <- numeric(0)
        while (length(jj <- select(a, prob = prob)) == 0) {
        }
        a <- jj
    }
    return(out)
}

"tee.fun" <-
function (x, D1, subsets, hpa) 
{
    if (!is.vector(x)) {
        out <- apply(x, 1, match.fun(sys.call()[[1]]), D1 = D1, 
            subsets = subsets, hpa = hpa)
        return(out)
    }
    rhos <- hpa$rhos
    sigma_squareds <- hpa$sigma_squareds
    x <- t(x)
    t.old <- prod(rhos) * sigma_squareds[1] * as.vector(Afun(level = 1, 
        Di = x, Dj = D1[subsets[[1]],,drop=FALSE ], hpa = hpa))
    out <- t.old
    s <- length(subsets)
    stopifnot(s > 1)
    for (i in 2:s) {
        t.dash <- t.old[match(subsets[[i]], subsets[[i - 1]])]
        t.new <- rhos[i - 1] * t.dash + sigma_squareds[i] * Pi(hpa, i, 
            s - 1) * as.vector(Afun(level = i, Di = x, Dj = D1[subsets[[i]], 
            , drop = FALSE], hpa = hpa))
        out <- c(out, t.new)
        t.old <- t.new
    }
    names(out) <- c(lapply(subsets,names),recursive=TRUE)
    return(out)
}

"opt.1" <-
function(D, z, basis, subsets, hpa.start, give.answers=FALSE, ...)
  {
    f <- function(candidate){
      hpa.temp <- hpa.start
      hpa.temp$B[[1]] <- diag(exp(candidate[1]),nrow=ncol(D))
      hpa.temp$sigma_squareds[1] <- exp(candidate[2])

      jj <- object(level=1,
             D=D,z=z, basis=basis, subsets=subsets, hpa=hpa.temp)
      return(jj)
    }

    start.point <- log(c(hpa.start$B[[1]][1,1], hpa.start$sigma_squareds[1]))
    optim.output <- optim(start.point, f, ...)
    u <- optim.output$par
    hpa.out <- hpa.start
    diag(hpa.out$B[[1]]) <- exp(u[1])
    hpa.out$sigma_squareds[1] <- exp(u[2])
    if(give.answers){
      return(list(optim.output=optim.output,hpa=hpa.out))
    } else {
      return(hpa.out)
    }
  }
    

#opt1(D=D1.toy,z=z.toy,basis=basis.toy,subsets=subsets.toy, hpa.start=hpa.toy,control=list(maxit=1))

"opt.gt.1" <-
function(level, D, z, basis, subsets, hpa.start, give.answers=FALSE,
...)
{
  if(level<2){stop("level should be 2 or greater")}
  
  f <- function(candidate){
    hpa.temp <- hpa.start
    hpa.temp$B[[level]] <- diag(exp(candidate[1]), nrow=ncol(D))
    hpa.temp$sigma_squareds[level] <- exp(candidate[2])
    hpa.temp$rhos[level-1] <- exp(candidate[3])
    object(level=level,
           D=D,z=z, basis=basis, subsets=subsets, hpa=hpa.temp)
  }
  start.point <- log(c(hpa.start$B[[level]][1,1], hpa.start$sigma_squareds[level], hpa.start$rhos[level-1]))
  optim.output <- optim(start.point, f, ...)
  u <- optim.output$par
  hpa.out <- hpa.start
  diag(hpa.out$B[[level]]) <- exp(u[1])
  hpa.out$sigma_squareds[level] <- exp(u[2])
  hpa.out$rhos[level-1] <- exp(u[3])
  if(give.answers){
    return(list(optim.output=optim.output,hpa.out=hpa.out))
  } else {
    return(hpa.out)
  }
}

"c_fun" <- function(x,xdash=x,subsets,hpa)
{
  s <- length(subsets)
  out <- 0
  for(i in 1:s){
    out <- out +
      Pi(hpa,i,s-1)^2*hpa$sigma_squareds[i]*Afun(level=i,x,xdash,hpa=hpa)
  }
  return(out)
}

"cdash.fun" <- function(x, xdash=x, V=NULL, Vinv=NULL, D1, subsets,
  basis,hpa, method=2){
  if(is.vector(x)){x <- t(x)}
  if(is.vector(xdash)){xdash <- t(xdash)}
  if(is.null(V)){
    V <- V.fun.app(D1=D1,subsets=subsets,hpa=hpa)
  }

  tx <- tee.fun(x,D1=D1,subsets=subsets,hpa=hpa)
  txdash <- tee.fun(xdash,D1=D1,subsets=subsets,hpa=hpa)
  hx <- hdash.fun(x=x,hpa=hpa,basis=basis)
  hxdash <- hdash.fun(x=xdash,hpa=hpa,basis=basis)
  H <- H.fun.app(D1=D1,subsets=subsets,basis=basis,hpa=hpa)

  cxx <- c_fun(x=x,xdash=xdash,subsets=subsets,hpa=hpa) 

  if(method == 1){
    if(is.null(Vinv)){
      Vinv <- solve(V)
    }
    cxxdash <-
      cxx - 
        quad.3form(Vinv,txdash,tx) + 
          quad.3form(
                     solve(quad.form(Vinv,H)),
                     hxdash-quad.3form(Vinv,H,txdash),
                     hx    -quad.3form(Vinv,H,tx    )
                     )
    return(cxxdash)
  } else if (method == 2){
    if(is.null(Vinv)){
      Vinv <- solve(V)
    }
    U <- crossprod(Vinv,H)
    cxxdash <- 
      cxx - 
        quad.3form(Vinv,txdash,tx) + 
          quad.3form(
                     solve(crossprod(H,U)),
                     hxdash - crossprod(U,txdash),
                     hx - crossprod(U,tx)
                     )
  } else if (method == 3){  # no matrix inversion
    U <- crossprod(V,H)
    Y <- solve(V,H)
    cxxdash <- 
      cxx -
        crossprod(txdash,solve(V,tx)) +
          crossprod(
          solve(crossprod(H,Y), hxdash - crossprod(Y,txdash)),
                        hx - crossprod(Y,tx)
                        )
  } else if (method ==4) {  # bog-standard method, a direct
                            # transliteration of KOH2000:
    if(is.null(Vinv)){
      Vinv <- solve(V)
    }
    cxxdash <- 
      cxx -
        t(txdash) %*% Vinv %*% tx + 
          t(hxdash-t(t(txdash) %*% Vinv %*% H)) %*%
            solve(t(H) %*% Vinv  %*% H) %*%
          (hx-t(t(tx) %*% Vinv %*% H))
    } else {
      stop("method must be an integer in the range 1-4")
    }
  return(cxxdash)
}

"basis.genie" <- 
function (x) 
{
    if (is.vector(x)) {
        stopifnot(length(x) == 4)
        out <- c(1, x)
        names(out) <- c("const", LETTERS[1:4])
        return(out)
    }
    else {
        return(t(apply(x, 1, match.fun(sys.call()[[1]]))))
    }
}

"hpa.fun.genie" <- 
function (x) 
{
    if (length(x) != 17) {
        stop("x must have 17 elements")
    }
    pdm.maker <- function(x) {
        jj <- diag(x[1:4])
        rownames(jj) <- c("InvDrag","FWF.scl","powercloud","fluxfactor")
        colnames(jj) <- rownames(jj)
        return(jj)
    }
    sigma_squareds <- x[1:3]
    names(sigma_squareds) <- paste("level", 1:3, sep = "")
    B <- NULL
    B[[1]] <- pdm.maker(x[4:7])
    B[[2]] <- pdm.maker(x[8:11])
    B[[3]] <- pdm.maker(x[12:15])
    names(B) <- paste("level", 1:3, sep = ".")
    rhos <- x[16:17]
    names(rhos) <- paste("level", 1:2, sep = ".")
    return(list(sigma_squareds = sigma_squareds, B = B, rhos = rhos))
}

"subset_maker" <-
function(x)
{
  out <- lapply(x,function(i){1:i})
  names(out) <- paste("level",1:length(out),sep=".")
  return(out)
}
  
