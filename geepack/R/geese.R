geese <- function(formula = formula(data),
                  sformula = ~ 1,
                  id, waves = NULL,
                  data = parent.frame(), subset = NULL, na.action = na.omit,
                  contrasts = NULL, weights = NULL,
                  ## zcor is design matrix for alpha,
                  ## corp is known paratemers to correlation coef. rho
                  zcor = NULL, corp = NULL,
                  ## zsca is constructed from sformula
                  ## control parameters
                  control = geese.control(...),
                  ## param 
                  b = NULL, alpha = NULL, gm = NULL,
                  ## geestr
                  family = gaussian(),
                  mean.link = NULL,
                  variance = NULL,
                  cor.link = "identity",
                  sca.link = "identity",
                  link.same = TRUE,
                  scale.fix = FALSE, scale.value = 1.0,
                  ## corr
                  corstr = "independence",
                  ...) {
  scall <- match.call()
  mnames <- c("", "formula", "data", "offset", "weights", "subset", "na.action", "id", "waves", "corp")
  cnames <- names(scall)
  cnames <- cnames[match(mnames,cnames,0)]
  mcall <- scall[cnames]
  if (is.null(mcall$id)) mcall$id <- as.name("id")
  mcall[[1]] <- as.name("model.frame")
  m <- eval(mcall, parent.frame())

  y <- model.extract(m, "response")
  if (is.null(dim(y))) N <- length(y) else N <- dim(y)[1]
  mterms <- attr(m, "terms")
  x <- model.matrix(mterms, m, contrasts)
  offset <- model.extract(m, "offset")
  if (is.null(offset)) offset <- rep(0, N)
  w <- model.extract(m, "weights")
  if (is.null(w)) w <- rep(1, N)
  id <- model.extract(m, id)
  waves <- model.extract(m, "waves")
  corp <- model.extract(m, "corp")
  if (is.null(id)) stop("id variable not found.")

  ##print(control)
  
  ## setting up the scale model;
  ## borrowed idea from S+ function dglm by Gordon Smyth
  mcall$formula <- formula
  mcall$formula[3] <- switch(match(length(sformula), c(0,2,3)),
                             1, sformula[2], sformula[3])
  m <- eval(mcall, parent.frame())
  terms <- attr(m, "terms")
  zsca <- model.matrix(terms, m, contrasts)
  soffset <- model.extract(m, "offset")
  if (is.null(soffset)) soffset <- rep(0, N)  
 
  if (is.character(family)) family <- get(family)
  if (is.function(family))  family <- family()
  ans <- geese.fit(x, y, id, offset, soffset, w,
                   waves, zsca, zcor, corp, 
                   control,
                   b, alpha, gm,
                   family, mean.link, variance, cor.link, sca.link,
                   link.same, scale.fix, scale.value, 
                   corstr, ...)
  ans <- c(ans, list(call=scall, formula=formula)) 
  class(ans) <- "geese"
  ans
}

geese.fit <- function(x, y, id,
                      offset=rep(0,N), soffset=rep(0,N), weights=rep(1,N),
                      waves = NULL, zsca = matrix(1,N,1),
                      zcor = NULL, corp = NULL,
                      control = geese.control(...),
                      ## param 
                      b = NULL, alpha = NULL, gm = NULL,
                      ## geestr
                      family = gaussian(),
                      mean.link = NULL,
                      variance = NULL,
                      cor.link = "identity",
                      sca.link = "identity",
                      link.same = TRUE,
                      scale.fix = FALSE, scale.value = 1.0,
                      ## corr
                      corstr = "independence", ...) {
  N <- length(id)
  ##clusz <- unlist(lapply(split(id, id), length))
  clusnew <- c(which(diff(as.numeric(id)) != 0), length(id))
  clusz <- c(clusnew[1], diff(clusnew))
  maxclsz <- max(clusz)
  if (is.null(waves)) waves <- unlist(sapply(clusz, function(x) 1:x))
  waves <- as.integer(waves)

  LINKS <- c("identity", "logit", "probit", "cloglog", "log", "inverse", "fisherz", "lwybc2", "lwylog")
  VARIANCES <- c("gaussian", "binomial", "poisson", "Gamma") ## quasi is not supported yet

  if (is.null(mean.link)) mean.link <- family$link
  if (is.null(variance)) variance <- family$family
  mean.link.v <- pmatch(mean.link, LINKS, -1, TRUE)
  cor.link.v <- pmatch(cor.link, LINKS, -1, TRUE)
  sca.link.v <- pmatch(sca.link, LINKS, -1, TRUE)
  variance.v <- pmatch(variance, VARIANCES, -1, TRUE)
  if (any(mean.link.v == -1)) stop("mean.link invalid.")
  if (any(cor.link.v == -1)) stop("cor.link invalid.")
  if (any(sca.link.v == -1)) stop("sca.link invalid.")
  if (any(variance.v == -1)) stop("variance invalid.")
  if (length(mean.link.v) != length(variance.v))
    stop("mean.link and variance not same length.")
  if (length(mean.link.v) != length(sca.link.v))
    stop("mean.link and sca.link not same lehgnt.")
      
  if (length(id) != length(y)) stop("id and y not same length.")
  if (length(offset) != length(y)) stop("offset and y not same length")
  if (length(soffset) != length(y)) stop("sca.offset and y not same length")
  if (nrow(zsca) != length(y)) stop("nrow(zsca) and length(y) not match")
  
  if (link.same) linkwaves <- rep(1, N)
  else {
    if (max(waves) != maxclsz) stop("maximum waves and maximum cluster size not equal")
    if (length(mean.link.v) != maxclsz) stop("length of mean.link not equal to the maximum cluster size.")
    linkwaves <- waves
  }
  linkwaves <- as.integer(linkwaves)
  geestr <- list(length(mean.link.v), as.integer(mean.link.v),
                 as.integer(variance.v), as.integer(sca.link.v),
                 as.integer(cor.link.v), as.integer(scale.fix))

  CORSTRS <- c("independence", "exchangeable", "ar1", "unstructured", "userdefined", "fixed")
  corstrv <- pmatch(corstr, CORSTRS, -1)
  if (corstrv == -1) stop("invalid corstr.")
  corr <- list(as.integer(corstrv), maxclsz)
  
  if (is.null(zcor)) {
    if (corstrv == 5) stop("need zcor matrix for userdefined corstr.") 
    else zcor <- genZcor(clusz, waves, corstrv)
  }
  else {
    if (!is.matrix(zcor)) zcor <- as.matrix(zcor)
    if (corstrv >= 4 && nrow(zcor) != sum(clusz * (clusz - 1) / 2)) stop("nrow(zcor) need to be equal sum(clusz * (clusz - 1) / 2) for unstructured or userdefined corstr.")
    if (corstrv %in% c(2,3) && nrow(zcor) != length(clusz)) stop("nrow(zcor) need to be equal to the number of clusters for exchangeable or ar1 corstr.")
  }
  if (!is.matrix(zcor)) zcor <- as.matrix(zcor)
  if (is.null(corp)) corp <- as.double(waves)

  p <- ncol(x)
  q <- ncol(zcor)
  r <- ncol(zsca)
  
  ## Initial values setup
  ## This may fail for binomial model with log link (relative risk)
  ## fit0 <- glm.fit(x, y, weights=weights, offset=offset, family=family)
  if (is.null(b)){
    ##b <- rep(1,p)
    fit0 <- glm.fit(x, y, weights=weights, offset=offset, family=family)
    b <- fit0$coef
  }
  if (is.null(alpha)) {
    if (corstrv == 6) alpha <- 1
    else alpha <- rep(0,q)
  }
  if (is.null(gm)) {
    ##gm <- rep(scale.value, r)
    qlf <- quasi(LINKS[sca.link.v])$linkfun
    ## pr2 <- (residuals.glm(fit0, type="pearson")) ^ 2
    mu <- quasi(LINKS[mean.link.v])$linkinv(x %*% b)
    pr2 <- (y - mu) ^ 2 / family$variance(mu)
    gm <- lm.fit(zsca, qlf(pr2), offset = soffset)$coef
  }
  param <- list(b, alpha, gm)

  ans <- .Call("gee_rap", y, x, offset, soffset, weights,
               linkwaves, zsca, zcor, corp,
               clusz, geestr, corr, param, control, PACKAGE = "geepack")
  names(ans) <- c("beta", "alpha", "gamma", "vbeta", "valpha", "vgamma",
                  "vbeta.naiv", "valpha.naiv", "valpha.stab",
                  "vbeta.ajs", "valpha.ajs", "vgamma.ajs",
                  "vbeta.j1s", "valpha.j1s", "vgamma.j1s",
                  "vbeta.fij", "valpha.fij", "vgamma.fij",
                  "error")
  ans$xnames <- dimnames(x)[[2]]
  ans$zsca.names <- dimnames(zsca)[[2]]
  ans$zcor.names <- dimnames(zcor)[[2]]
  if (is.null(ans$zcor.names)) ans$zcor.names = paste("alpha", 1:ncol(zcor), sep=":")
  names(ans$beta) <- ans$xnames
  names(ans$gamma) <- ans$zsca.names
  if (length(ans$alpha) > 0)  names(ans$alpha) <- ans$zcor.names

  param <- list(ans$beta, ans$alpha, ans$gamma)
  infls <- .Call("infls_rap",  y, x, offset, soffset, weights,
               linkwaves, zsca, zcor, corp,
               clusz, geestr, corr, param, control, PACKAGE = "geepack")
  rownames(infls) <- c(paste("beta", names(ans$beta), sep="_"),
                       if (length(ans$gamma) > 0) paste("gamma", names(ans$gamma), sep="_") else NULL,
                       if (length(ans$alpha) > 0) paste("alpha", names(ans$alpha), sep="_") else NULL)

  ans <- c(ans,
           list(infls=infls,
                clusz=clusz, control=control,
                model=list(mean.link=mean.link,
                  variance=variance, sca.link=sca.link,
                  cor.link=cor.link, corstr=corstr, scale.fix=scale.fix)))
  ans
}

geese.control <- function (epsilon = 1e-04, maxit = 25, trace = FALSE,
                           scale.fix = FALSE, jack = FALSE,
                           j1s = FALSE, fij = FALSE) {
  if (!is.numeric(epsilon) || epsilon <= 0) 
    stop("value of epsilon must be > 0")
  if (!is.numeric(maxit) || maxit <= 0) 
    stop("maximum number of iterations must be > 0")
  list(trace = as.integer(trace),
       jack = as.integer(jack), j1s = as.integer(j1s), fij = as.integer(fij),
       maxit = as.integer(maxit), epsilon = epsilon)
}


## compare coefficients
compCoef <- function(fit0, fit1) {
  v0 <- names(fit0$beta)
  v1 <- names(fit1$beta)
  v0idx <- (1:length(v0))[v0 %in% v1]
  v1idx <- (1:length(v1))[v1 %in% v0]
  delta <- fit0$beta[v0idx] - fit1$beta[v1idx]
  infls <- fit0$infls[v0idx,] - fit1$infls[v1idx,]
  robvar <- infls %*% t(infls)
  list(delta = delta, variance = robvar)  
}
