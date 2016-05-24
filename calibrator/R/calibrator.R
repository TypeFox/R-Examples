"C1" <-
function (D1, D2, theta, phi) 
{
    omega_x <- phi$omega_x
    omega_t <- phi$omega_t
    sigma1squared <- phi$sigma1squared
    pos.def.matrix <- blockdiag(omega_x, omega_t)
    corr.matrix(xold = D1, yold = D2.fun(D2 = D2, theta = theta), 
        pos.def.matrix = pos.def.matrix)*
        sigma1squared
}
"D1.fun" <-
function (x.star, t.vec) 
{
    if (is.vector(x.star) & is.vector(t.vec)) {
        x.star <- as.matrix(t(x.star))
        t.vec <- as.matrix(t(t.vec))
    }
    if (is.vector(t.vec)) {
        t.vec <- do.call("rbind", lapply(1:nrow(x.star), function(...) {
            t.vec
        }))
    }
    if (is.vector(x.star)) {
        x.star <- do.call("rbind", lapply(1:nrow(t.vec), function(...) {
            x.star
        }))
    }
    D1 <- cbind(x.star, t.vec)
    if (is.null(rownames(x.star)) & is.null(rownames(x.star))) {
        rownames(D1) <- paste("code", 1:nrow(t.vec), sep = ".")
    }
    else {
        if (is.null(x.star)) {
            rownames(D1) <- rownames(t.vec)
        }
        else {
            rownames(D1) <- rownames(x.star)
        }
    }
    return(D1)
}
"D2.fun" <-
function (D2, theta) 
{
    jj <- t(array(theta, c(length(theta), nrow(D2))))
    colnames(jj) <- names(theta)
    jj <- cbind(D2, jj)
    code.ind <- ncol(D2) + 1:length(theta)
    colnames(jj)[code.ind] <- LETTERS[1:length(theta)]
    if (nrow(jj) > 0) {
        rownames(jj) <- paste("obs", 1:nrow(jj), sep = ".")
    }
    return(jj)
}
"E.theta.toy" <-
function (D2 = NULL, H1 = NULL, x1 = NULL, x2 = NULL, phi, give.mean = TRUE) 
{
    if (give.mean) {
        m_theta <- phi$theta.apriori$mean
        return(H1(D1.fun(D2, t.vec = m_theta)))
    }
    else {
        out <- matrix(0, 6, 6)
        out[4:6, 4:6] <- phi$theta.apriori$sigma
        return(out)
    }
}
"Edash.theta.toy" <-
function (x, t.vec, k, H1, fast.but.opaque = FALSE, a = NULL, 
    b = NULL, phi = NULL) 
{
    if (fast.but.opaque) {
        edash.mean <- a + crossprod(b, t.vec[k, ])
    }
    else {
        V_theta <- phi$theta.apriori$sigma
        m_theta <- phi$theta.apriori$mean
        omega_t <- phi$omega_t
        edash.mean <- solve(solve(V_theta) + 2 * omega_t, solve(V_theta, 
            m_theta) + 2 * crossprod(omega_t, t.vec[k, ]))
    }
    jj <- as.vector(edash.mean)
    names(jj) <- rownames(edash.mean)
    edash.mean <- jj
    return(H1(D1.fun(x, edash.mean)))
}

"hbar.fun.toy" <- function(theta, X.dist, phi){
  if(is.vector(theta)){theta <- t(theta)}
  first.bit <- phi$rho*H1.toy(D1.fun(X.dist$mean,theta))
  second.bit <- H2.toy(X.dist$mean)
  jj.names <- colnames(second.bit)
  second.bit <- kronecker(second.bit,rep(1,nrow(first.bit)))
  colnames(second.bit) <- jj.names
  return(t(cbind(first.bit, second.bit)))
}
 
"Ez.eqn7.supp" <-
function (z, D1, H1, D2, H2, extractor, beta2, y, E.theta, phi) 
{
    rho <- phi$rho
    m_theta <- phi$theta.apriori$mean
    V_theta <- phi$theta.apriori$sigma
    expectation <- E.theta(D2 = D2, H1 = H1, phi = phi)
    beta1hat <- beta1hat.fun(D1 = D1, H1 = H1, y = y, phi = phi)
    f <- function(i) {
        t.fun(x = D2[i, ], D1 = D1, extractor = extractor, phi = phi)
    }
    bit1 <- H2(D2) %*% beta2
    bit2 <- rho * expectation %*% beta1hat
    bit3 <- rho * crossprod(sapply(1:nrow(D2), f), solve(V1(D1, 
        phi = phi), y - H1(D1) %*% beta1hat))
    return(bit1 + bit2 + bit3)
}
"Ez.eqn9.supp" <-
function (x, theta, d, D1, D2, H1, H2, phi) 
{
    if (is.vector(theta)) {
        return(Ez.eqn9.supp.vector(x = x, theta = theta, d = d, 
            D1 = D1, D2 = D2, H1 = H1, H2 = H2, phi = phi))
    }
    f <- function(jj.theta) {
        Ez.eqn9.supp.vector(x = x, theta = jj.theta, d = d, D1 = D1, 
            D2 = D2, H1 = H1, H2 = H2, phi = phi)
    }
    out <- apply(theta, 1, f)
    rownames(out) <- rownames(x)
    return(out)
}

"h.fun" <-
function(x, theta, H1, H2, phi) {
    rho <- phi$rho
    t(cbind(rho * H1(D1.fun(x.star = x, t.vec = theta)), 
            H2(x)))
  }

"Ez.eqn9.supp.vector" <-
function (x, theta, d, D1, D2, H1, H2, phi) 
{
    bhat <- betahat.fun.koh(theta = theta, d = d, D1 = D1, D2 = D2, 
        H1 = H1, H2 = H2, phi = phi)
    h.x.th <- h.fun(x, theta = theta, H1 = H1, H2 = H2,phi=phi)
    H.th <- H.fun(theta = theta, D1 = D1, D2 = D2, H1 = H1, H2 = H2, 
        phi = phi)
    t.x.th <- tee(x, theta = theta, D1 = D1, D2 = D2, phi = phi)
    Vd.th <- Vd(theta = theta, D1 = D1, D2 = D2, phi = phi)
    out <- crossprod(h.x.th, bhat) + crossprod(t.x.th, solve(Vd.th, 
        d - H.th %*% bhat))
      return(out)
}

"Cov.eqn9.supp" <-
function (x, xdash=NULL, theta, d, D1, D2, H1, H2, phi) 
{
  if(is.null(xdash)){xdash <- x}
  bit1 <- V1(D1.fun(x,theta),D1.fun(xdash,theta), phi=phi)* phi$rho^2
  bit2 <- V2(x,xdash,phi=phi)
  
  t.x.th <-     tee(x=x    , theta=theta,D1=D1,D2=D2,phi=phi)
  t.xdash.th <- tee(x=xdash, theta=theta,D1=D1,D2=D2,phi=phi)
    
  Vd.inv <- solve(Vd(D1=D1,D2=D2,theta=theta,phi=phi))
  
  jj.fun <- function(x){
    h.fun(x=x,theta=theta,H1=H1,H2=H2,phi=phi) -
      quad.3form(
                 M     = Vd.inv,
                 left  = H.fun(theta=theta, D1=D1, D2=D2, H1=H1, H2=H2, phi=phi),
                 right = tee(x=x,theta=theta, D1=D1, D2=D2, phi=phi)
                 )
  }
 
  bit3 <-
    quad.3form(
               M     = Vd.inv,
               left  = t.x.th,
               right = t.xdash.th
               )

  bit4 <-
    quad.3form(
               M = W(D1=D1, D2=D2, H1=H1, H2=H2, theta=theta, det=FALSE, phi=phi),
               left  = jj.fun(x),
               right = jj.fun(xdash)
               )
  return(bit1 + bit2 - bit3 + bit4) 
}


"H.fun" <-
function (theta, D1, D2, H1, H2, phi) 
{
    rho <- phi$rho
    top.left <- H1(D1)
    low.left <- H1(D1.fun(D2, theta)) * rho
    low.right <- H2(D2)
    top.right <- matrix(0, nrow = nrow(top.left), ncol = ncol(low.right))
    out <- rbind(cbind(top.left, top.right), cbind(low.left, 
        low.right))
    colnames(out)[(ncol(top.left) + 1):ncol(out)] <- colnames(low.right)
    return(out)
}
"H1.toy" <-
function (D1) 
{
    if (is.vector(D1)) {
        D1 <- t(D1)
    }
    out <- t(apply(D1, 1, h1.toy))
    colnames(out)[1] <- "h1.const"
    return(out)
}
"H2.toy" <-
function (D2) 
{
    if (is.vector(D2)) {
        D2 <- t(D2)
    }
    out <- t(apply(D2, 1, h2.toy))
    colnames(out) <- names(h2.toy(D2[1,,drop=TRUE]))
    return(out)
}
"V.fun" <-
function (D1, D2, H1, H2, extractor, E.theta, Edash.theta, give.answers = FALSE, 
    test.for.symmetry = FALSE, phi) 
{

    lambda <- phi$lambda
    rho <- phi$rho
    x.star <- extractor(D1)$x.star
    t.vec <- extractor(D1)$t.vec
    v1.d1 <- V1(D1 = D1, other = NULL, phi = phi)
    v1.d1.inv <- solve(v1.d1)
    w1.d1 <- W1(D1 = D1, H1 = H1, phi = phi)
    h1.d1 <- H1(D1)

    jj <- crossprod(h1.d1,v1.d1.inv)
    w1.blah <- crossprod(w1.d1, jj)
#    v1.blah <- crossprod(v1.d1.inv, h1.d1) %*% w1.d1
    v1.blah <- crossprod(jj, w1.d1)
    v1.blah.di.blah <- v1.blah %*% jj

    line1 <- line2 <- line3 <- line4 <- line5 <- line6 <- matrix(0, 
        nrow(D2), nrow(D2))
    tvec.arbitrary <- phi$theta.apriori$mean
    for (i in 1:nrow(D2)) {
        if (test.for.symmetry == FALSE) {
            j.vector <- i:nrow(D2)
        }
        else {
            j.vector <- 1:nrow(D2)
        }
        for (j in j.vector) {
            line1[i, j] <- V1(D1 = D1.fun(x.star = D2[i, ], t.vec = tvec.arbitrary), 
                other = D1.fun(x.star = D2[j,,drop=TRUE], t.vec = tvec.arbitrary), 
                phi = phi)
            jj.tt <- tt.fun(D1 = D1, extractor = extractor, x.i =
                            D2[i,,drop=TRUE], x.j = D2[j,,drop=TRUE], phi = phi)
            line2[i, j] <- sum(rowSums(v1.d1.inv * jj.tt))
            jj.hh <- hh.fun(x.i = D2[i,,drop=TRUE], x.j = D2[j,,drop=TRUE], E.theta = E.theta, 
                H1 = H1, phi = phi)
            line3[i, j] <- sum(rowSums(w1.d1 * jj.hh))
            jj.th <- ht.fun(x.i = D2[j,,drop=TRUE], x.j = D2[i,,drop=TRUE], D1 = D1, 
                extractor = extractor, Edash.theta = Edash.theta, 
                H1 = H1, fast.but.opaque = TRUE, x.star = x.star, 
                t.vec = t.vec, phi = phi)
            line4[i, j] <- sum(rowSums(w1.blah * t(jj.th)))
            jj.ht <- ht.fun(x.i = D2[i,,drop=TRUE], x.j = D2[j,,drop=TRUE], D1 = D1, 
                extractor = extractor, Edash.theta = Edash.theta, 
                H1 = H1, fast.but.opaque = TRUE, x.star = x.star, 
                t.vec = t.vec, phi = phi)
            line5[i, j] <- sum(rowSums(v1.blah * jj.ht))
            line6[i, j] <- sum(rowSums(jj.tt * v1.blah.di.blah))
        } 
    }
    C.m <- +line1 - line2 + line3 - line4 - line5 + line6
    if (test.for.symmetry == FALSE) {
        C.m <- symmetrize(C.m)
    }
    C.m <- rho^2 * C.m
    rownames(C.m) <- rownames(D2)
    colnames(C.m) <- rownames(D2)
    C.lambda <- lambda * diag(nrow = nrow(D2))
    C.v2 <- V2(D2, other = NULL, phi = phi)
    out <- C.lambda + C.v2 + C.m
    if (give.answers) {
        attributes(line1) <- attributes(out)
        attributes(line2) <- attributes(out)
        attributes(line3) <- attributes(out)
        attributes(line4) <- attributes(out)
        attributes(line5) <- attributes(out)
        attributes(line6) <- attributes(out)
        return(list(line1 = line1, line2 = line2, line3 = line3, 
            line4 = line4, line5 = line5, line6 = line6, v1.d1 = v1.d1, 
            v1.d1.inv = v1.d1.inv, w1.d1 = w1.d1, h1.d1 = h1.d1, 
            w1.blah = w1.blah, v1.blah = v1.blah, v1.blah.di.blah = v1.blah.di.blah, 
            C.m = C.m, C.lambda = C.lambda, C.v2 = C.v2, out = out))
    }
    else {
        return(out)
    }
}
"V1" <-
function (D1, other = NULL, phi) 
{
    return(phi$sigma1squared * t(corr.matrix(xold = D1, yold = other, 
        pos.def.matrix = blockdiag(phi$omega_x, phi$omega_t))))

}
"V2" <-
function (x, other = NULL, phi) 
{
    if (is.vector(x)) {
        x <- as.data.frame(t(x))
    }
    return(t(phi$sigma2squared * corr.matrix(xold = x, yold = other, 
        pos.def.matrix = phi$omegastar_x)))
}
"Vd" <-
function (D1, D2, theta, phi) 
{
    rho <- phi$rho
    lambda <- phi$lambda
    top.left <- V1(D1 = D1, other = NULL, phi = phi)
    top.right <- rho * t(C1(D1 = D1, D2 = D2, theta = theta, 
        phi = phi))
    low.left <- t(top.right)
    low.right <- lambda * diag(nrow(low.left)) + rho^2 * V1(D1 = D1.fun(D2, 
        theta), other = NULL, phi = phi) + V2(x = D2, other = NULL, 
        phi = phi)
    return(rbind(cbind(top.left, top.right), cbind(low.left, 
        low.right)))
}
"W" <-
function (D1, D2, H1, H2, theta, det = FALSE, phi) 
{
    jj.H <- H.fun(theta = theta, D1 = D1, D2 = D2, H1 = H1, H2 = H2, 
        phi = phi)
    jj.V <- Vd(theta = theta, D1 = D1, D2 = D2, phi = phi)
    out <- quad.form.inv(jj.V, jj.H)
    if (det) {
        return(1/det(out))
    }
    else {
        return(solve(out))
    }
}
"W1" <-
function (D1, H1, det = FALSE, phi) 
{
    out <- quad.form.inv(V1(D1, phi = phi), H1(D1))
    if (det) {
        return(1/det(out))
    }
    else {
        return(solve(out))
    }
}
"W2" <-
function (D2, H2, V, det = FALSE) 
{
    out <- as.matrix(quad.form(solve(V), H2(D2)))
    if (det) {
        return(1/det(out))
    }
    else {
        return(solve(out))
    }
}
"beta1hat.fun" <-
function (D1, H1, y, phi) 
{
    out <- crossprod(W1(D1 = D1, H1 = H1, phi = phi), crossprod(H1(D1), 
        solve(V1(D1 = D1, phi = phi), y)))
    return(out)
}
"beta2hat.fun" <-
function (D1, D2, H1, H2, V = NULL, z, etahat.d2, extractor, 
    E.theta, Edash.theta, phi) 
{
    rho <- phi$rho
    if (is.null(V)) {
        V <- V.fun(D1 = D1, D2 = D2, H1 = H1, H2 = H2, extractor = extractor, 
            E.theta = E.theta, Edash.theta = Edash.theta, phi = phi)
    }
    out <- crossprod(W2(D2, H2, V), crossprod(H2(D2), solve(V, 
        z - rho * etahat.d2)))
    return(out)
}
"betahat.fun.koh" <-
function (D1, D2, H1, H2, theta, d, phi) 
{
    if (is.vector(theta)) {
        return(betahat.fun.koh.vector(D1 = D1, D2 = D2, H1 = H1, 
            H2 = H2, theta = theta, d = d, phi = phi))
    }
    f <- function(jj.theta) {
        betahat.fun.koh.vector(D1 = D1, D2 = D2, H1 = H1, H2 = H2, 
            theta = jj.theta, d = d, phi = phi)
    }
    out <- apply(theta, 1, f)
    rownames(out) <- colnames(H.fun(theta = theta[1, ], D1 = D1, 
        D2 = D2, H1 = H1, H2 = H2, phi = phi))
    return(out)
}
"betahat.fun.koh.vector" <-
function (D1, D2, H1, H2, theta, d, phi) 
{
    jj.W <- W(theta, D1 = D1, D2 = D2, H1 = H1, H2 = H2, phi = phi)
    jj.H <- H.fun(theta, D1 = D1, D2 = D2, H1 = H1, H2 = H2, 
        phi = phi)
    jj.V <- Vd(theta, D1 = D1, D2 = D2, phi = phi)
    out <- crossprod(jj.W, crossprod(jj.H, solve(jj.V, d)))
    return(out)
}
"blockdiag" <-
function (m1, m2, p.tr = 0, p.ll = 0) 
{
    topleft <- m1
    topright <- matrix(p.tr, nrow(m1), ncol(m2))
    colnames(topright) <- colnames(m2)
    lowleft <- matrix(p.ll, nrow(m2), ncol(m1))
    lowright <- m2
    rbind(cbind(topleft, topright), cbind(lowleft, lowright))
}
"computer.model" <-
function (X, params = NULL, set.seed.to.zero = TRUE,
          draw.from.prior=FALSE, export.true.hyperparameters = FALSE, phi=NULL) 
{
    if(draw.from.prior){
      jj.psi1 <-
      rmvnorm(n=1,mean=phi$psi1.apriori$mean,sigma=phi$psi1.apriori$sigma) 
      REAL.SIGMA1SQUARED <- c(sigma1squared = jj.psi1[4])
      REAL.SCALES <- c(x=jj.psi1[1], y=jj.psi1[2],
      A=jj.psi1[3],B=jj.psi1[4], C=jj.psi1[5])
    } else {
      REAL.SIGMA1SQUARED <- c(sigma1squared = 0.3)
      REAL.SCALES <- c(x = 10, y = 0.1, A = 10, B = 1, C = 0.1)
    }
    REAL.COEFFS <- c(const.co.coeff = 4, x.co.coeff = 5, y.co.coeff = 6, 
                     A.co.coeff = 7, B.co.coeff = 8, C.co.coeff = 9)

    if (export.true.hyperparameters) {
        warning("Only omniscient programmers should set export.true.hyperparameters to TRUE")
        return(list(REAL.SIGMA1SQUARED = REAL.SIGMA1SQUARED, 
            REAL.SCALES = REAL.SCALES, REAL.COEFFS = REAL.COEFFS))
    }
    if (set.seed.to.zero) {
        set.seed(0)
    }
    
    if (is.vector(X) & is.vector(params)) {
        X <- t(X)
        params <- t(params)
    }
    if (is.vector(X)) {
        X <- kronecker(t(X), rep(1, nrow(params)))
        rownames(X) <- rownames(params)
    }
    if (ncol(X) == 5) {
        if (!is.null(params)) {
            stop("X has 5 columns; params argument superfluous")
        }
        params <- X[, 3:5]
        X <- X[, 1:2]
    }
    if (is.vector(params)) {
        params <- kronecker(t(params), rep(1, nrow(X)))
        rownames(params) <- rownames(X)
    }
    x <- X[, 1]
    y <- X[, 2]
    A <- params[, 1]
    B <- params[, 2]
    C <- params[, 3]
    jj.mean <- as.vector(H1.toy(cbind(X, params)) %*% REAL.COEFFS)
    names(jj.mean) <- rownames(params)
    D <- cbind(X, params)
    var.matrix <- REAL.SIGMA1SQUARED * corr.matrix(D, scales = REAL.SCALES)
    out <- rmvnorm(n = 1, mean = jj.mean, sigma = as.matrix(var.matrix))
    out <- as.vector(out)
    names(out) <- rownames(X)
    return(out)
}
"cov.p5.supp" <-
function (x, xdash=NULL, theta, d, D1, D2, H1, H2, phi) 
{
  if(!is.vector(x)){stop("x should be a vector; consider Cov.eqn9.supp() for vector x and xdash with a single value for theta")}
  x <- t(x)
  if(is.null(xdash)){xdash <- x}
  f <- function(theta){Cov.eqn9.supp(
                                     x=x,
                                     xdash=xdash,
                                     theta=theta,
                                     d=d,
                                     D1=D1,
                                     D2=D2,
                                     H1=H1,
                                     H2=H2,
                                     phi=phi
                                     )}
  apply(theta,1,f)
}
"create.new.toy.datasets" <-
function (D1,D2,export=FALSE)
{

    REAL.PARAMS <- c(A = 0.9, B = 0.1, C = 0.3)
    REAL.RHO <- 1
    REAL.LAMBDA <- 0.5

        if(export){
          return(list(
                      REAL.PARAMS = REAL.PARAMS,
                      REAH.RHO = REAL.RHO,
                      REAL.LAMBDA = REAL.LAMBDA
                      ))
        }
      
    N <- nrow(D1)
    n <- nrow(D2)

    
    D1.total <- rbind(D1,
                      cbind(D2,kronecker(rep(1,n),t(REAL.PARAMS)) )
    )

    jj <- computer.model(D1.total)
    
    y.toy <- jj[1:N]

    i <- (N+1):(N+n)

    error.term <- rnorm(n=n, mean=0,sd=REAL.LAMBDA)
    
    z.toy <- jj[i]*REAL.RHO + model.inadequacy(X=D1.total[i,1:2]) + error.term
    d.toy <- c(y.toy, z.toy)
    return(list(y.toy = y.toy, z.toy = z.toy, d.toy = d.toy))
  }
"dists.2frames" <-
function (a, b = NULL, A = NULL, A.lower = NULL, test.for.symmetry = TRUE) 
{
    if (is.null(b)) {
        b <- a
    }
    if (!is.null(A)) {
        distance <- function(x, y) {
            m <- x - y
            return(as.vector(crossprod(t(crossprod(m, A)), m)))
        }
    }
    else {
        distance <- function(x, y) {
            m <- x - y
            jj <- crossprod(A.lower, m)
            return(crossprod(jj, jj))
        }
    }
    if (test.for.symmetry == TRUE) {
        return(apply(a, 1, function(y) {
            apply(b, 1, function(x) {
                distance(x, y)
            })
        }))
    }
    else {
        out <- matrix(0, nrow(a), nrow(b))
        for (i in 1:nrow(a)) {
            for (j in i:nrow(b)) {
                out[i, j] <- distance(a[i, ], b[j, ])
            }
        }
        return(symmetrize(out))
    }
}

"etahat" <-
function (D1, D2, H1, y, E.theta, extractor, phi) 
{
    if(is.function(E.theta)){
      jj.hfun <- E.theta(D2=D2, H1=H1, phi=phi, give.mean=TRUE)
    } else {
      jj.hfun <- H1(D1.fun(x.star=D2, t.vec=E.theta))
    }
    jj.beta <- beta1hat.fun(D1=D1,H1=H1,y=y,phi=phi)
    bit1 <- drop(jj.hfun %*% jj.beta)

    b <- beta1hat.fun(D1=D1, H1=H1, y=y, phi=phi)
    f <- function(x){ t.fun(x, D1=D1, extractor=extractor, phi=phi) }

    jj.tfun <- apply(D2,1,f)
    bit2 <- crossprod(jj.tfun, solve(V1(D1, phi=phi),y - H1(D1) %*% b))
    return(drop(bit1 + bit2))
}
"extractor.toy" <-
function (D1) 
{
    return(list(x.star = D1[, 1:2, drop = FALSE], t.vec = D1[, 
        3:5, drop = FALSE]))
}
"h1.toy" <-
function (x) 
{
    x <- c(1, x)
    names(x)[1] <- "const"
    return(x)
}
"h2.toy" <-
function (x) 
{
#    x <- c(x[1] > 0.5, x[2] > 0.5)
    x <- c(1,x[1],x[2]^2)
    names(x) <- c("constant.h2","alpha","beta")
    return(x)
}
"hh.fun" <-
function (x.i, x.j, H1, E.theta, phi) 
{
    h.i <- E.theta(D2 = x.i, H1 = H1, phi = phi)
    h.j <- E.theta(D2 = x.j, H1 = H1, phi = phi)
    out <- crossprod(h.i, h.j)
    out <- out + E.theta(x1 = x.i, x2 = x.j, give.mean = FALSE, 
        phi = phi)
    return(out)
}
"ht.fun" <-
function (x.i, x.j, D1, extractor, Edash.theta, H1, fast.but.opaque = TRUE, 
    x.star = NULL, t.vec = NULL, phi) 
{
    V_theta <- phi$theta.apriori$sigma
    m_theta <- phi$theta.apriori$mean
    omega_x <- phi$omega_x
    omega_t <- phi$omega_t
    if (fast.but.opaque == TRUE) {
        if (is.null(x.star) | is.null(t.vec)) {
            stop("argument 'fast.but.opaque' being TRUE means caller has to supply x.star and t.vec (use extractor.toy() for this)")
        }
    }
    else {
        x.star <- extractor(D1)$x.star
        t.vec <- extractor(D1)$t.vec
    }
    f1 <- function(k) {
        as.vector(quad.form(omega_x, x.i - x.star[k, ]))
    }
    matrix.in.f2 <- 2 * V_theta + solve(omega_t)
    f2 <- function(k) {
        as.vector(quad.form.inv(matrix.in.f2, m_theta - t.vec[k, 
            ]))
    }
    jj.a <- phi$a
    jj.b <- phi$b
    jj.c <- phi$c
    f3 <- function(k) {
        if (fast.but.opaque == TRUE) {
            return(Edash.theta(x = x.j, t.vec = t.vec, k = k, 
                H1 = H1, fast.but.opaque = TRUE, a = jj.a, b = jj.b, 
                phi = phi))
        }
        else {
            return(Edash.theta(x = x.j, t.vec = t.vec, k = k, 
                H1 = H1, phi = phi))
        }
    }
    out <- sapply(1:nrow(t.vec), function(i) {
        jj.c * exp(-f1(i) - f2(i)) * f3(i)
    })
    rownames(out) <- colnames(f3(1))
    colnames(out) <- rownames(t.vec)
    return(t(out))
}
"is.positive.definite" <-
function (a, ...) 
{
    all(eigen(a, only.values = TRUE, ...)$values > 0)
}
"p.eqn4.supp" <-
function (D1, y, H1, include.prior = TRUE, lognormally.distributed, 
    return.log = FALSE, phi) 
{
    v1d1 <- V1(D1 = D1, phi = phi)
    w1d1 <- W1(D1 = D1, H1 = H1, phi = phi)
    vec <- y - H1(D1) %*% beta1hat.fun(D1, H1, y, phi)
    if (return.log) {
        jj <- as.vector(-quad.form.inv(v1d1, vec)/2) + log(det(w1d1))/2 - 
            log(det(v1d1))/2
        if (include.prior) {
            return(jj + log(prob.psi1(phi, lognormally.distributed = lognormally.distributed)))
        }
        else {
            return(jj)
        }
    }
    else {
        jj <- as.vector(exp(-1/2 * quad.form.inv(v1d1, vec))) * 
            sqrt(det(w1d1)/det(v1d1))
        if (include.prior) {
            return(jj * prob.psi1(phi, lognormally.distributed = lognormally.distributed))
        }
        else {
            return(jj)
        }
    }
}
"p.eqn8.supp" <-
function (theta, D1, D2, H1, H2, d, include.prior = FALSE, lognormally.distributed = FALSE, 
    return.log = FALSE, phi) 
{
    if (is.vector(theta)) {
        return(p.eqn8.supp.vector(theta = theta, D1 = D1, D2 = D2, 
            H1 = H1, H2 = H2, d = d, include.prior = include.prior, 
            lognormally.distributed = lognormally.distributed, 
            return.log = return.log, phi = phi))
    }
    f <- function(jj.theta) {
        return(p.eqn8.supp.vector(theta = jj.theta, D1 = D1, 
            D2 = D2, H1 = H1, H2 = H2, d = d, include.prior = include.prior, 
            lognormally.distributed = lognormally.distributed, 
            return.log = return.log, phi = phi))
    }
    out <- apply(theta, 1, f)
    return(out)
}
"p.eqn8.supp.vector" <-
function (theta, D1, D2, H1, H2, d, include.prior = FALSE, lognormally.distributed = FALSE, 
    return.log = FALSE, phi) 
{
    betahat <- betahat.fun.koh(theta = theta, D1 = D1, D2 = D2, 
        H1 = H1, H2 = H2, phi = phi, d = d)
    Vd.theta <- Vd(theta = theta, D1 = D1, D2 = D2, phi = phi)
    det.W.theta <- W(theta = theta, D1 = D1, D2 = D2, H1 = H1, 
        H2 = H2, det = TRUE, phi = phi)
    H.theta <- H.fun(theta = theta, D1 = D1, D2 = D2, H1 = H1, 
        H2 = H2, phi = phi)
    mismatch <- d - H.theta %*% betahat
    if (return.log) {
        out <- -log(det(Vd.theta))/2 + log(det.W.theta)/2 - quad.form.inv(Vd.theta, 
            mismatch)/2
        if (include.prior) {
            out <- out + log(prob.theta(theta, phi, lognormally.distributed = lognormally.distributed))
        }
    }
    else {
        out <- det(Vd.theta)^(-1/2) * det.W.theta^(1/2) * exp(-0.5 * 
            quad.form.inv(Vd.theta, mismatch))
        if (include.prior) {
            out <- out * prob.theta(theta, phi, lognormally.distributed = lognormally.distributed)
        }
    }
    return(as.vector(out))
}
"p.page4" <-
function (D1, D2, H1, H2, V, y, z, E.theta, Edash.theta, extractor, include.prior = FALSE, 
    lognormally.distributed = FALSE, return.log = FALSE, phi) 
{
    if(is.null(V)){
      V <-
        V.fun(D1=D1,D2=D2,H1=H1,H2=H2,extractor=extractor,E.theta=E.theta,
              Edash.theta=Edash.theta, give.answers=FALSE,phi=phi)
    } 
    etahat.d2 <- etahat(D1 = D1, D2=D2, H1 = H1, y=y,
        E.theta = E.theta, extractor = extractor, phi = phi)
    beta2hat <- beta2hat.fun(D1 = D1, D2 = D2, H1=H1, H2 = H2, V = V, 
        z = z, etahat.d2 = etahat.d2, extractor=extractor,
    E.theta=E.theta, Edash.theta=Edash.theta, phi = phi)
    h2.d2 <- H2(D2)
    rho <- phi$rho
    vec <- z - h2.d2 %*% beta2hat - rho * etahat.d2
    if (return.log) {
        bit1 <- log(W2(D2, H2, V, det = TRUE))/2 - log(det(V))/2
        bit2 <- -0.5 * quad.form.inv(V, vec)
        if (include.prior) {
            return(bit1 + bit2 + log(prob.psi2(phi, lognormally.distributed = lognormally.distributed)))
        }
        else {
            return(bit1 + bit2)
        }
    }
    else {
        bit1 <- sqrt(W2(D2, H2, V, det = TRUE)/det(V))
        bit2 <- exp(-0.5 * quad.form.inv(V, vec))
        if (include.prior) {
            return(bit1 * bit2 * prob.psi2(phi, lognormally.distributed = lognormally.distributed))
        }
        else {
            return(bit1 * bit2)
        }
    }
}
"phi.change" <-
function (phi.fun, old.phi = NULL, rho = NULL, lambda = NULL, 
    psi1 = NULL, psi1.apriori = NULL, psi1.apriori.mean = NULL, 
    psi1.apriori.sigma = NULL, psi2 = NULL, psi2.apriori = NULL, 
    psi2.apriori.mean = NULL, psi2.apriori.sigma = NULL, theta.apriori = NULL, 
    theta.apriori.mean = NULL, theta.apriori.sigma = NULL) 
{
    if (is.null(rho)) {
        rho <- old.phi$rho
    }
    if (is.null(lambda)) {
        lambda <- old.phi$lambda
    }
    if (is.null(psi1)) {
        psi1 <- old.phi$psi1
    }
    else {
        if (is.null(names(psi1))) {
            names(psi1) <- names(old.phi$psi1)
        }
    }
    if (is.null(psi1.apriori)) {
        psi1.apriori <- old.phi$psi1.apriori
    }
    if (!is.null(psi1.apriori.mean)) {
        psi1.apriori$mean <- psi1.apriori.mean
        if (is.null(names(psi1.apriori.mean))) {
            names(psi1.apriori$mean) <- names(old.phi$psi1.apriori$mean)
        }
    }
    if (!is.null(psi1.apriori.sigma)) {
        psi1.apriori$sigma <- psi1.apriori.sigma
        if (is.null(rownames(psi1.apriori.sigma))) {
            rownames(psi1.apriori$sigma) <- rownames(old.phi$psi1.apriori$sigma)
        }
        if (is.null(colnames(psi1.apriori.sigma))) {
            colnames(psi1.apriori$sigma) <- colnames(old.phi$psi1.apriori$sigma)
        }
    }
    if (is.null(psi2)) {
        psi2 <- old.phi$psi2
    }
    else {
        if (is.null(names(psi2))) {
            names(psi2) <- names(old.phi$psi2)
        }
    }
    if (is.null(psi2.apriori)) {
        psi2.apriori <- old.phi$psi2.apriori
        if (!is.null(psi2.apriori.mean)) {
            psi2.apriori$mean <- psi2.apriori.mean
            if (is.null(names(psi2.apriori.mean))) {
                names(psi2.apriori$mean) <- names(old.phi$psi2.apriori$mean)
            }
        }
        if (!is.null(psi2.apriori.sigma)) {
            psi2.apriori$sigma <- psi2.apriori.sigma
            if (is.null(rownames(psi2.apriori.sigma))) {
                rownames(psi2.apriori$sigma) <- rownames(old.phi$psi2.apriori$sigma)
            }
            if (is.null(colnames(psi2.apriori.sigma))) {
                colnames(psi2.apriori$sigma) <- colnames(old.phi$psi2.apriori$sigma)
            }
        }
    }
    if (is.null(theta.apriori)) {
        theta.apriori <- old.phi$theta.apriori
        if (!is.null(theta.apriori.mean)) {
            theta.apriori$mean <- theta.apriori.mean
            if (is.null(names(theta.apriori.mean))) {
                names(theta.apriori$mean) <- names(old.phi$theta.apriori$mean)
            }
        }
        if (!is.null(theta.apriori.sigma)) {
            theta.apriori$sigma <- theta.apriori.sigma
            if (is.null(rownames(theta.apriori.sigma))) {
                rownames(theta.apriori$sigma) <- rownames(old.phi$theta.apriori$sigma)
            }
            if (is.null(colnames(theta.apriori.sigma))) {
                colnames(theta.apriori$sigma) <- colnames(old.phi$theta.apriori$sigma)
            }
        }
    }
    return(phi.fun(rho = rho, lambda = lambda, psi1 = psi1, psi1.apriori = psi1.apriori, 
        psi2 = psi2, psi2.apriori = psi2.apriori, theta.apriori = theta.apriori))
}
"phi.fun.toy" <-
function (rho, lambda, psi1, psi1.apriori, psi2, psi2.apriori, 
    theta.apriori) 
{
    "pdm.maker.psi1" <- function(psi1) {
        jj.omega_x <- diag(psi1[1:2])
        rownames(jj.omega_x) <- names(psi1[1:2])
        colnames(jj.omega_x) <- names(psi1[1:2])
        jj.omega_t <- diag(psi1[3:5])
        rownames(jj.omega_t) <- names(psi1[3:5])
        colnames(jj.omega_t) <- names(psi1[3:5])
        sigma1squared <- psi1[6]
        return(list(omega_x = jj.omega_x, omega_t = jj.omega_t, 
            sigma1squared = sigma1squared))
    }
    "pdm.maker.psi2" <- function(psi2) {
        jj.omegastar_x <- diag(psi2[1:2])
        sigma2squared <- psi2[3]
        return(list(omegastar_x = jj.omegastar_x, sigma2squared = sigma2squared))
    }
    jj.mean <- theta.apriori$mean
    jj.V_theta <- theta.apriori$sigma
    jj.discard.psi1 <- pdm.maker.psi1(psi1)
    jj.omega_t <- jj.discard.psi1$omega_t
    jj.omega_x <- jj.discard.psi1$omega_x
    jj.sigma1squared <- jj.discard.psi1$sigma1squared
    jj.discard.psi2 <- pdm.maker.psi2(psi2)
    jj.omegastar_x <- jj.discard.psi2$omegastar_x
    jj.sigma2squared <- jj.discard.psi2$sigma2squared
    jj.omega_t.upper <- chol(jj.omega_t)
    jj.omega_t.lower <- t(jj.omega_t.upper)
    jj.omega_x.upper <- chol(jj.omega_x)
    jj.omega_x.lower <- t(jj.omega_x.upper)
    jj.a <- solve(solve(jj.V_theta) + 2 * jj.omega_t, solve(jj.V_theta, 
        jj.mean))
    jj.b <- t(2 * solve(solve(jj.V_theta) + 2 * jj.omega_t) %*% 
        jj.omega_t)
    jj.c <- jj.sigma1squared/sqrt(det(diag(nrow = nrow(jj.V_theta)) + 
        2 * jj.V_theta %*% jj.omega_t))
    names(jj.c) <- "ht.fun.precalc"
    jj.A <- solve(jj.V_theta + solve(jj.omega_t)/4)
    jj.A.upper <- chol(jj.A)
    jj.A.lower <- t(jj.A.upper)
    list(rho = rho, lambda = lambda, psi1 = psi1, psi1.apriori = psi1.apriori, 
        psi2 = psi2, psi2.apriori = psi2.apriori, theta.apriori = theta.apriori, 
        omega_x = jj.omega_x, omega_t = jj.omega_t, 
        omegastar_x = jj.omegastar_x, sigma1squared = jj.sigma1squared, 
        sigma2squared = jj.sigma2squared, omega_x.upper = jj.omega_x.upper, 
        omega_x.lower = jj.omega_x.lower, omega_t.upper = jj.omega_t.upper, 
        omega_t.lower = jj.omega_t.lower, a = jj.a, b = jj.b, 
        c = jj.c, A = jj.A, A.upper = jj.A.upper, A.lower = jj.A.lower)
}
"phi.true.toy" <-
function (phi) 
{
    re <- model.inadequacy(export.true.hyperparameters = TRUE)
    co <- computer.model(export.true.hyperparameters = TRUE)
    phi.change(old.phi = phi, phi.fun = phi.fun.toy, rho = re$REAL.RHO, 
        lambda = re$REAL.LAMBDA, psi1 = c(co$REAL.SCALES, co$REAL.SIGMA1SQUARED), 
        psi2 = c(re$REAL.ROUGHNESS, re$REAL.SIGMA2SQUARED))
}
"prob.psi1" <-
function (phi, lognormally.distributed = TRUE) 
{
    if (lognormally.distributed) {
        if (any(phi$psi1 <= 0)) {
            return(0)
        }
        return(dmvnorm(x = log(phi$psi1), mean = phi$psi1.apriori$mean, 
            sigma = phi$psi1.apriori$sigma))
    }
    else {
        return(dmvnorm(x = phi$psi1, mean = phi$psi1.apriori$mean, 
            sigma = phi$psi1.apriori$sigma))
    }
}
"prob.psi2" <-
function (phi, lognormally.distributed = TRUE) 
{
    if (lognormally.distributed) {
        if (any(phi$psi2 <= 0)) {
            return(0)
        }
        return(dmvnorm(x = log(c(phi$rho, phi$lambda, phi$psi2)), 
            mean = phi$psi2.apriori$mean, sigma = phi$psi2.apriori$sigma))
    }
    else {
        return(dmvnorm(x = c(phi$rho, phi$lambda, phi$psi2), 
            mean = phi$psi2.apriori$mean, sigma = phi$psi2.apriori$sigma))
    }
}
"prob.theta" <-
function (theta, phi, lognormally.distributed = FALSE) 
{
    if (lognormally.distributed) {
        if (any(theta) < 0) {
            return(0)
        }
        return(dmvnorm(x = log(theta), mean = log(phi$theta.aprior$mean), 
            sigma = phi$theta.aprior$sigma))
    }
    else {
        return(dmvnorm(x = theta, mean = phi$theta.aprior$mean, 
            sigma = phi$theta.aprior$sigma))
    }
}
"model.inadequacy" <-
function (X, set.seed.to.zero = TRUE, draw.from.prior=FALSE,
    export.true.hyperparameters = FALSE, phi=NULL) 
{
      if (set.seed.to.zero) {
        set.seed(0)
      }
      
      if(draw.from.prior){
        jj.params <- rmvnorm(n=1, mean=phi$theta.apriori$mean,
                             sigma=phi$theta.apriori$sigma)
        jj.psi2 <- rmvnorm(n=1, mean=phi$psi2.apriori$mean ,
                           sigma=phi$psi2.apriori$sigma)
        
        REAL.BETA2 <- c(int.re.coeff = 0.4, tt = -0.8,fish=1)
        REAL.SIGMA2SQUARED <- c(sigma2squared=jj.psi2[5])
        REAL.ROUGHNESS <- c(x=jj.psi2[3],y=jj.psi2[4])
      } else {
        REAL.BETA2 <- c(int.re.coeff = 0.4, tt = -0.8,fish=1)
        REAL.SIGMA2SQUARED <- c(sigma2squared = 0.13)
        REAL.ROUGHNESS <- c(x = 2, y = 3)
      }
    
    if (export.true.hyperparameters) {
        warning("Only omniscient statisticians should access the true hyperparameters")
        return(list(
                    REAL.BETA2 = REAL.BETA2, 
                    REAL.SIGMA2SQUARED = REAL.SIGMA2SQUARED,
                    REAL.ROUGHNESS = REAL.ROUGHNESS
                    ))

    }
    if (is.vector(X)) {
        X <- t(X)
    }
      out <- as.vector(rmvnorm(    n = 1,
                                mean = as.vector(H2.toy(X) %*% REAL.BETA2),
                               sigma = REAL.SIGMA2SQUARED * as.matrix(corr.matrix(X, 
                                 scales = REAL.ROUGHNESS)
                               )
                       ))
      return(out)
}
"sample.theta" <-
function (n = 1, phi) 
{
    out <- rmvnorm(n = n, mean = phi$theta.apriori$mean, sigma = phi$theta.apriori$sigma)
    colnames(out) <- names(phi$theta.apriori$mean)
    return(out)
}
"stage1" <-
function (D1, y, H1, maxit, trace=100, method = "Nelder-Mead", directory = ".",
          do.filewrite=FALSE, do.print=TRUE,
    phi.fun, lognormally.distributed = FALSE, include.prior = TRUE, 
    phi) 
{
    if(do.filewrite & !isTRUE(file.info(directory)$isdir)){
      stop("do.filewrite = TRUE; directory name supplied does not exist")
    }
       
    f <- function(candidate) {
        phi.temp <- phi.change(phi.fun = phi.fun, old.phi = phi, 
            psi1 = exp(candidate))
        f.out <- p.eqn4.supp(D1 = D1, y = y, H1 = H1, lognormally.distributed = lognormally.distributed, 
            include.prior = include.prior, return.log = TRUE, 
            phi = phi.temp)
        if (do.print) {
             print(candidate)
             print(f.out)
           }
        if (do.filewrite){
          save(
               phi.temp, f.out, ascii = TRUE, file =
               file.path(directory,  paste("stage1.",gsub("[ :]", "_", date()),sep=""))
               )
          save(
               phi.temp, f.out, ascii = TRUE, file =
               file.path(directory, "stage1.latest")
               )
        }
    return(f.out)
  }
   psi1.start <- log(phi$psi1)
    e <- optim(psi1.start, f, method = method, control = list(fnscale = -1, 
        maxit = maxit, trace = trace))
    phi.out <- phi.change(phi.fun = phi.fun, old.phi = phi, psi1 = exp(e$par))
    return(invisible(phi.out))
}

"stage2" <-
  function (D1, D2, H1, H2, y, z, maxit, trace=100, method = "Nelder-Mead",
            directory = ".", do.filewrite=FALSE, do.print=TRUE,
            extractor, phi.fun, E.theta, Edash.theta, isotropic = FALSE,
            lognormally.distributed = FALSE, include.prior = TRUE,
            use.standin = FALSE, rho.eq.1 = TRUE, phi) 
{
    if(do.filewrite & !isTRUE(file.info(directory)$isdir)){
      stop("do.filewrite = TRUE; directory name supplied does not exist")
    }
    f <- function(candidate) {
      if(isotropic){
        phi.temp <- phi.change(phi.fun = phi.fun, old.phi = phi, 
                               rho = exp(candidate[1]), lambda = exp(candidate[2]), 
                               psi2 = c(rep(exp(candidate[3]),nrow(phi$omegastar_x)),candidate[4])
                               )
      } else {
        phi.temp <- phi.change(phi.fun = phi.fun, old.phi = phi, 
                               rho = exp(candidate[1]), lambda = exp(candidate[2]), 
                               psi2 = exp(candidate[-c(1,2)])
                               )
      }
      if (use.standin) {
        V.temp <- 1 + diag(3, nrow = nrow(D2))
      } else {
        V.temp <- V.fun(D1 = D1, D2 = D2, H1 = H1, H2 = H2, 
                        extractor = extractor, E.theta = E.theta, Edash.theta = Edash.theta, 
                        phi = phi.temp)
      }
      if(rho.eq.1){
        phi.temp <- phi.change(phi.fun = phi.fun, old.phi = phi, rho=1)
      }
      f.out <- drop(p.page4(D1 = D1, D2 = D2, H1 = H1, H2 = H2, V = V.temp,
                            z = z, y = y, E.theta = E.theta, Edash.theta=Edash.theta,
                            extractor=extractor,include.prior = include.prior, 
                            lognormally.distributed = lognormally.distributed, 
                            return.log = TRUE, phi = phi.temp)
                    )
        if (do.print) {
          print(candidate)
          print(f.out)
        }
        if(do.filewrite){
          save(phi.temp, f.out, ascii = TRUE, file =
               file.path(directory,  paste("stage2.",gsub("[ :]", "_", date()),sep=""))
               )
          save(
               phi.temp, f.out, ascii = TRUE, file =
               file.path(directory, "stage2.latest")
               )
        }
        return(f.out)
    }
    if(isotropic){
      rho.lambda.psi2.start <- log(c(phi$rho, phi$lambda,  phi$psi2[1], phi$sigma2squared))
    } else {
      rho.lambda.psi2.start <- log(c(phi$rho, phi$lambda, phi$psi2))
    }
    e <- optim(rho.lambda.psi2.start, f, method = method,
               control = list(fnscale = -1, trace = trace, maxit = maxit))
    jj <- exp(e$par)
    if(isotropic){
      phi.out <- phi.change(phi.fun = phi.fun, old.phi = phi, rho = jj[1], 
                            lambda = jj[2], psi2 = c(rep(jj[3],nrow(phi$omegastar_x)),jj[4])
                            )
    } else {
      phi.out <- phi.change(phi.fun = phi.fun, old.phi = phi, rho = jj[1], 
                            lambda = jj[2], psi2 = jj[-c(1,2)]
                            )
    }
    if(rho.eq.1){
      phi.out <- phi.change(phi.fun = phi.fun, old.phi=phi.out, rho=1)
    }
    return(invisible(phi.out))
}

"stage3" <-
function (D1, D2, H1, H2, d, maxit, trace=100, method = "Nelder-Mead", directory
        = ".", do.filewrite=FALSE, do.print=TRUE,
    include.prior = TRUE, lognormally.distributed = FALSE, theta.start = NULL, 
    phi) 
{

  if(do.filewrite & !isTRUE(file.info(directory)$isdir)){
    stop("do.filewrite = TRUE; directory name supplied does not exist")
  }
 
    f <- function(theta) {
        print(theta)
        f.out <- p.eqn8.supp(theta, D1 = D1, D2 = D2, H1 = H1, 
            H2 = H2, d = d, include.prior = include.prior, lognormally.distributed = lognormally.distributed, 
            return.log = TRUE, phi = phi)
        if (do.print) {
            print(theta)
            print(f.out)
          }
        if(do.filewrite){
            save(theta, f.out, ascii = TRUE, file =
               file.path(directory,  paste("stage3.",gsub("[ :]", "_", date()),sep="")))
            save(theta, f.out, ascii = TRUE, file =
               file.path(directory, "stage3.latest"))
        }
        return(f.out)
    }
    if (is.null(theta.start)) {
        theta.start <- phi$theta.apriori$mean
    }
    e <- optim(theta.start, f, method = method, control = list(fnscale = -1, 
        maxit = maxit, trace = trace))
    optimal.theta <- e$par
    names(optimal.theta) <- names(phi$theta.apriori$mean)
    return(optimal.theta)
}
"symmetrize" <-
function (a) 
{
    a + t(a) - diag(diag(a))
}
"t.fun" <-
function (x, D1, extractor, phi) 
{
    sigma1squared <- phi$sigma1squared
    rho <- phi$rho
    omega_x <- phi$omega_x
    omega_t <- phi$omega_t
    x.star <- extractor(D1)$x.star
    t.vec <- extractor(D1)$t.vec
    m_theta <- phi$theta.apriori$mean
    V_theta <- phi$theta.apriori$sigma
    jj.mat <- solve(2 * V_theta + solve(omega_t))
    prod.1 <- sigma1squared/sqrt(det(diag(nrow = nrow(omega_t)) + 
        2 * crossprod(V_theta, omega_t)))
    out <- rep(NA, nrow(D1))
    for (j in 1:length(out)) {
        jj.x <- as.vector(x.star[j, ] - x)
        jj.m <- as.vector(m_theta - t.vec[j, ])
        prod.2 <- exp(-quad.form(omega_x, jj.x))
        prod.3 <- exp(-quad.form(jj.mat, jj.m))
        out[j] <- prod.1 * prod.2 * prod.3
    }
    names(out) <- rownames(t.vec)
    return(out)
}
"tee" <-
function (x, theta, D1, D2, phi) 
{
    if (is.vector(x)) {
        x <- t(x)
    }
    rho <- phi$rho
    jj.v1 <- D1.fun(x, theta)
    d2.theta <- D2.fun(D2 = D2, theta = theta)
    jj.D1 <- do.call("rbind", lapply(1:nrow(d2.theta), function(i) {
        jj.v1
    }))
    bit1 <- rho * V1(D1 = jj.v1, other = D1, phi = phi)
    bit2a <- rho^2 * V1(D1 = jj.v1, other = d2.theta, phi = phi)
    bit2b <- V2(x, D2, phi = phi)
    out <- t(cbind(bit1, bit2a + bit2b))
    return(out)
}


"tt.fun" <-
function (D1, extractor, x.i, x.j, test.for.symmetry = FALSE, 
    method = 1, phi) 
{
    V_theta <- phi$theta.apriori$sigma
    m_theta <- phi$theta.apriori$mean
    x.star <- extractor(D1)$x.star
    t.vec <- extractor(D1)$t.vec
    omega_x <- phi$omega_x
    omega_t <- phi$omega_t
    sigma1squared <- phi$sigma1squared
    A.lower <- phi$A.lower
    if (method == 0) {
        x.jframe <- kronecker(rep(1, nrow(t.vec)), t(x.j))
        x.iframe <- kronecker(rep(1, nrow(t.vec)), t(x.i))
        attributes(x.jframe) <- attributes(x.star)
        attributes(x.iframe) <- attributes(x.star)
        exp.x.jk <- exp(-dists.2frames(x.jframe, x.star, omega_x))
        exp.x.il <- exp(-dists.2frames(x.star, x.iframe, omega_x))
        exp.t.kl <- exp(-dists.2frames(t.vec, t.vec, omega_t)/2)
    }
    else {
        exp.x.jk <- kronecker(exp(-dists.2frames(t(x.j), x.star, 
            omega_x)), t(rep(1, nrow(D1))))
        exp.x.il <- kronecker(exp(-t(dists.2frames(x.star, t(x.i), 
            omega_x))), rep(1, nrow(D1)))
        exp.t.kl <- exp(-dists.2frames(t.vec, t.vec, omega_t)/2)
        rownames(exp.x.jk) <- rownames(D1)
        colnames(exp.x.jk) <- rownames(D1)
        rownames(exp.x.il) <- rownames(D1)
        colnames(exp.x.il) <- rownames(D1)
    }
    mtheta.minus.tmean <- matrix(0, nrow(x.star), nrow(x.star))
    for (k in 1:nrow(x.star)) {
        if (test.for.symmetry == TRUE) {
            l.vector <- 1:nrow(x.star)
        }
        else {
            l.vector <- k:nrow(x.star)
        }
        for (l in l.vector) {
            jj <- m_theta - (t.vec[k, ] + t.vec[l, ])/2
            A <- phi$A
            mtheta.minus.tmean[k, l] <- exp(-as.vector(crossprod(jj, 
                A) %*% jj)/2)
        }
    }
    if (test.for.symmetry == FALSE) {
        mtheta.minus.tmean <- symmetrize(mtheta.minus.tmean)
    }
    fish <- 1/sqrt(det(diag(nrow = nrow(V_theta)) + 4 * V_theta %*% 
        omega_t))
    out <- sigma1squared^2 * fish * exp.x.jk * exp.x.il * exp.t.kl * 
        mtheta.minus.tmean
    colnames(out) <- rownames(out)
    return(out)
}

"EK.eqn10.supp" <- function(X.dist, D1, D2, H1, H2, d, hbar.fun,
  lower.theta, upper.theta, extractor, give.info=FALSE, include.prior=FALSE, phi, ...){


  if(!exists("adapt")){
    adapt <- function(...){stop("not the original adapt")}
    stop("the adapt package is no longer available on CRAN: so
          the adapt() function, which is needed for calculate_B(),
          is not available either.

  You may be able to install the adapt package notwithstanding its
  availability on CRAN or is license.  If you are happy with this (I
  am), everything should work.

I am working on providing a replacement for adapt(), but this  is
low on my list of priorities.  Sorry about this.")

  }      


  
  scalar1 <-
    phi$sigma1squared/sqrt(det(diag(nrow=nrow(phi$omega_x)) +
                                       2*X.dist$var %*% phi$omega_x))
  scalar2 <-
    phi$sigma2squared/sqrt(det(diag(nrow=nrow(phi$omegastar_x)) +
                                   2*X.dist$var %*% phi$omegastar_x))
  
  matrix.1 <- solve(2*X.dist$var + solve(phi$omega_x))
  matrix.2 <- solve(2*X.dist$var + solve(phi$omegastar_x))
  
  "t1bar" <- function(theta, D1, D2, X.dist){ # Eqn 14 of the supplement
    f1 <- function(i){
      as.vector(quad.form(phi$omega_t, theta - t.vec[i,]))
    }
    f2 <- function(i){
      as.vector(quad.form(matrix.1, X.dist$mean - x.star[i,]))
    }
    out <- sapply(1:nrow(t.vec), function(i) { exp(-f1(i) - f2(i)) })
    return( (phi$rho) * scalar1 * out)
  }
  
  "t2bar" <- function(D1, D2, X.dist){ # This is eqn 15 of the
                                       # supplement.  Note absence of
                                       # theta!
    f1 <- function(i){
      as.vector(quad.form(matrix.1 , X.dist$mean - D2[i,]))
    }
    f2 <- function(i){
      as.vector(quad.form(matrix.2,  X.dist$mean - D2[i,]))
    }
    
    out1 <- sapply(1:nrow(D2), function(i) { exp(-f1(i)) })
    out2 <- sapply(1:nrow(D2), function(i) { exp(-f2(i)) })
    
    return( (phi$rho)^2*scalar1*out1 + scalar2*out2)
  }
  
  "tbar.fun" <- function(theta, D1, D2, X.dist, phi){
    c(t1bar(theta=theta, D1=D1, D2=D2, X.dist=X.dist),
      t2bar(D1=D1, D2=D2, X.dist=X.dist))
  }
  
  "integrand.numerator" <- function(theta){
    bhat <- betahat.fun.koh(D1=D1, D2=D2, H1=H1, H2=H2,
                            theta=theta,d=d,phi=phi)
    hbar <- hbar.fun(theta=theta, X.dist=X.dist, phi=phi)
    tbar <- tbar.fun(theta=theta, D1=D1, D2=D2, X.dist=X.dist, phi=phi)
    
    out <- crossprod(hbar,bhat) +
      crossprod(tbar, solve(Vd(D1=D1, D2=D2, theta=theta, phi=phi),
                            d-H.fun(theta=theta, D1=D1, D2=D2, H1=H1,
      H2=H2, phi=phi) %*% bhat))
    
    out <- out*p.eqn8.supp(theta=theta, D1=D1, D2=D2, H1=H1, H2=H2,
                           d=d, include.prior = include.prior,
                           lognormally.distributed = FALSE,
                           return.log = FALSE, phi=phi)
    return(out)
  }

  "integrand.denominator" <- function(theta){
    p.eqn8.supp(theta=theta, D1=D1, D2=D2, H1=H1, H2=H2, d=d,
                include.prior = include.prior,
                lognormally.distributed = FALSE,
                return.log = FALSE, phi=phi)
  }
  
  t.vec <- extractor(D1)$t.vec
  x.star <- extractor(D1)$x.star

  multi.dimensional <- length(phi$theta.apriori$mean)>1
  if(multi.dimensional){ # multi dimensional case

    
    numerator   <- adaptIntegrate(f=integrand.numerator  , lowerLimit = lower.theta, upperLimit = upper.theta, ...)
    denominator <- adaptIntegrate(f=integrand.denominator, lowerLimit = lower.theta, upperLimit = upper.theta, ...)

    out <- numerator$integral/denominator$integral
    
  } else { # one dimensional case
    integrand.numerator.vectorized <- function(theta.vec){
      sapply(theta.vec,integrand.numerator)
    }
    numerator <-
      integrate(f=integrand.numerator.vectorized,
                lower=lower.theta, upper=upper.theta, ...)
    
    
    integrand.denominator.vectorized <- function(theta.vec){
      sapply(theta.vec,integrand.denominator)
    }
    denominator <-
      integrate(f=integrand.denominator.vectorized,
                lower=lower.theta, upper=upper.theta, ...)
    
    out <- numerator$value/denominator$value
  }
  

  if(give.info){
    return(list(answer=out, numerator=numerator, denominator=denominator))
  } else {
    return(out)
  }    
}

"MH" <- function(n, start, sigma, pi){
  out <- matrix(NA,n,length(start))
  out[1,] <- start
  
  for(i in 2:n){  
    proposed <- drop(out[i-1,] + rmvnorm(n=1,sigma=sigma))
    num <- pi(proposed)
    den <- pi(out[i-1,])

    if( (num==0) & (den==0)){
      print("this cannot happen")
      alpha <- 0
    } else {
      alpha <- min(1, num/den)
    }
    if(runif(1)<alpha){
      out[i,] <- proposed
    } else {
      out[i,] <- out[i-1,]
    }
  }
  return(out)
}
