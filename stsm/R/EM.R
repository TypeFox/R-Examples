
##FIXME see other interface for parallel in windows, Sys.info()["sysname"]=="Windows"

maxlik.em <- function(model, type = c("standard", "modified", "mix"), 
  tol = 0.001, maxiter = 300, kfconv = c(0, 10, 1),  
  ur.maxiter = 1000, r.interval = c(0.001, var(model@y)*0.75, 20),
  mod.steps = seq(3, maxiter, 10), 
  parallel = FALSE, num.cores = NULL)
{
  type <- match.arg(type)[1]

  switch(type, 
    "standard" = res <- EM.v1(model, tol, maxiter, kfconv),
    "modified" = res <- EM.v2(model, tol, maxiter, ur.maxiter, 
      r.interval, kfconv, parallel, num.cores),
    "mix" = {
      mix <- list(mod.steps = mod.steps, nj = 1)
      res <- EM.v4(model, tol, maxiter, mix, ur.maxiter, 
        r.interval, kfconv, parallel, num.cores) })

  res
}

# the functions below are extracted from the 'stsm.em' package
# in this package they are not part of the namespace
# the wrapper 'maxlik.em' defined above must be used instead

elogdens.grad.v2 <- function(x = NULL, model, id = 'var1', kfconv = c(0.001, 10, 1))
{
  if (!is.null(x))
  {
      model@pars[id] <- x
  }

  if (!(id %in% names(model@pars))) {
    stop("'id' is not the name of a non-fixed parameter of the model.")
  }

  n <- length(model@y)
  pars <- c(model@pars, model@nopars)
  ss <- char2numeric(model)
  r <- ncol(ss$V)
  ss$ZtZ <- crossprod(ss$Z)

  if (id == 'var1')
  {
    ds <- .C("KFKSDS_steady_v2eps", 
      dim = as.integer(c(n, 1, ncol(ss$Z), r, 1, -1)), 
      sy = as.numeric(model@y), sZ = as.numeric(ss$Z), sZtZ = as.numeric(ss$ZtZ), 
      sT = as.numeric(t(ss$T)), sH = as.numeric(ss$H), sR = as.numeric(t(ss$R)), 
      V = as.numeric(ss$V[1,1]), sQ = as.numeric(ss$Q), 
      sa0 = as.numeric(ss$a0), sP0 = as.numeric(ss$P0),
      tol = as.numeric(kfconv[1]), maxiter = as.integer(kfconv[2]),
      ksconvfactor = as.numeric(kfconv[3]),
      res = double(1),  
      PACKAGE = "stsm")

    dvar <- ds$res[1]
  } else
  {
    id3 <- which(diag(model@ss$V) == id)
    id2 <- id3 - 1
    ds <- .C("KFKSDS_steady_v2eta", 
      dim = as.integer(c(n, 1, ncol(ss$Z), r, 1, -1, id2)),  
      sy = as.numeric(model@y), sZ = as.numeric(ss$Z), sZtZ = as.numeric(ss$ZtZ), 
      sT = as.numeric(t(ss$T)), sH = as.numeric(ss$H), sR = as.numeric(t(ss$R)), 
      V = as.numeric(ss$V[id3,id3]), sQ = as.numeric(ss$Q), 
      sa0 = as.numeric(ss$a0), sP0 = as.numeric(ss$P0),
      tol = as.numeric(kfconv[1]), maxiter = as.integer(kfconv[2]),
      ksconvfactor = as.numeric(kfconv[3]),
      res = double(1),
      PACKAGE = "stsm")

    dvar <- ds$res[1]
  }
  names(dvar) <- id

  dvar
}

step.EM.v1 <- function(model, kfconv)
{
  n <- length(model@y)
  pars <- c(model@pars, model@nopars)
  ss <- char2numeric(model)
  r <- ncol(ss$V)
  nr <- n*r

  ds <- .C("KFKSDS_steady", 
    dim = as.integer(c(n, 1, ncol(ss$Z), r, 1, -1, 0)), 
    sy = as.numeric(model@y), sZ = as.numeric(ss$Z), sT = as.numeric(t(ss$T)), 
    sH = as.numeric(ss$H), sR = as.numeric(t(ss$R)), sV = as.numeric(ss$V), 
    sQ = as.numeric(ss$Q), sa0 = as.numeric(ss$a0), sP0 = as.numeric(ss$P0),
    tol = as.numeric(kfconv[1]), maxiter = as.integer(kfconv[2]), 
    ksconvfactor = as.numeric(kfconv[3]),
    mll = double(1), epshat = double(n), vareps = double(n),
    etahat = double(nr), vareta = double(nr),
    sumepsmisc = double(1), sumetamisc = double(r),
    PACKAGE = "stsm")

  names(ds$sumepsmisc) <- 'var1'
  vareps <- ds$sumepsmisc / n
  vareps <- vareps[which(names(vareps) %in% names(model@pars))]

  names(ds$sumetamisc) <- paste('var', seq(2,r+1), sep = "")
  varseta <- ds$sumetamisc / (n-1)
  varseta <- varseta[which(names(varseta) %in% names(model@pars))]

  newpars <- c(vareps, varseta)

  newpars
}

EM.v1 <- function(model, tol = 0.001, maxiter = 300, 
  kfconv = c(0.001, 10, 1))
{
  Mpars <- rbind(model@pars)
  iter <- 0
  cond <- TRUE

  while (cond && iter < maxiter)
  {
    oldpars <- model@pars

    newpars <- step.EM.v1(model, kfconv)
    
    Mpars <- rbind(Mpars, newpars)
    model <- set.pars(model, newpars, FALSE)

    cond <- sqrt(sum((newpars - oldpars)^2)) > tol
    iter <- iter + 1
  }

  list(Mpars = Mpars, pars = newpars, iter = iter)
}

EM.v2 <- function(model, tol = 0.001, maxiter = 300, ur.maxiter = 1000, 
   r.interval = c(0.001, var(model@y)*0.75, 20), kfconv = c(0.001, 10, 1),
   parallel = FALSE, num.cores = NULL)
{
  newpars <- model@pars
  Mpars <- rbind(newpars)
  calls.v1 <- newpars
  calls.v1[] <- 0
  crep1 <- rep(1, length(calls.v1))
  iter <- 0 
  cond <- TRUE

  iv1 <- rep(r.interval[1], length(newpars))
  iv2 <- rep(r.interval[2], length(newpars))
  iv <-  r.interval[3]

  tmp <- as.list(model@pars)
  nms <- names(tmp)

  if (parallel)
  {
    names(iv2) <- names(iv1) <- nms

    search.root <- function(id)
    {
      #id <- nms[i]
      res <- try(uniroot(f = elogdens.grad.v2,
        interval = c(iv1[id], iv2[id]),
        model = model, id = id, kfconv = kfconv,
        maxiter = ur.maxiter)$root, silent = TRUE)

      if (inherits(res, "numeric")) {
        return(res)
      } else
      if (inherits(res, "try-error")) {
        return(NA)  # unlist(list(a=1,b=NULL,c=3)) would remove the second element
      } else stop("unexpected class of 'tmp'.")
    }

    if (is.null(num.cores))
      num.cores <- detectCores()
  }

  while (cond && iter < maxiter)
  {
    oldpars <- model@pars

    if (iter > iv)
    {
      ivaux <- Mpars[seq(iter-iv, iter+1),]
      iv1 <- apply(ivaux, 2, min)
      iv2 <- apply(ivaux, 2, max)
    }

    if (!parallel)
    {
      for (i in seq(along = tmp))
      {
        id <- nms[i]

        tmp[[id]] <- try(uniroot(f = elogdens.grad.v2,
          interval = c(iv1[i], iv2[i]),
          model = model, id = id, kfconv = kfconv,
          maxiter = ur.maxiter)$root, 
          silent = TRUE)

        if (class(tmp[[id]]) == "numeric") {
          newpars[id] <- tmp[[id]]
        } else
        if (class(tmp[[id]]) == "try-error") { # do nothing here
        } else stop("unexpected class of 'tmp'.")
      }

      ref <- !unlist(lapply(tmp, is.numeric))      

    } else {
      tmp <- mclapply(X = nms, FUN = search.root, mc.cores = num.cores)
      tmp <- unlist(tmp)
      ref <- is.na(tmp)
      id <- !ref
      if (any(id))
        newpars[id] <- tmp[id]
    }

    # after updating the parameters where a root was found,
    # swith to step.EM.v1 for those parameters where a root was not found

#     ref <- unlist(lapply(tmp, is.numeric))
#     ref <- !ref
    if (any(ref))
    {
      newpars[ref] <- step.EM.v1(model, kfconv)[ref]
      calls.v1[ref] <- calls.v1[ref] + crep1[ref]
    }

    Mpars <- rbind(Mpars, newpars)
    model <- set.pars(model, newpars, FALSE)

    cond <- sqrt(sum((newpars - oldpars)^2)) > tol
    iter <- iter + 1
  }

  list(Mpars = Mpars, pars = newpars, iter = iter, calls.v1 = calls.v1)
}

EM.v4 <- function(model, tol = 0.001, maxiter = 300,   
  opt = list(mod.steps = seq(3, maxiter, 10), nj = 1),
  ur.maxiter = 1000, r.interval = c(0.001, var(model@y)*0.75, 20),
  kfconv = c(0.001, 10, 1),
  parallel = FALSE, num.cores = NULL)
{
  Mpars <- rbind(model@pars)
  calls.v1 <- model@pars
  calls.v1[] <- 0
  iter <- 0
  cond <- TRUE

  ref <- opt$mod.steps
  nj <- opt$nj

  tmp <- as.list(model@pars)

  while (cond && iter < maxiter)
  {
    oldpars <- model@pars

    if (iter %in% ref)
    {
      sout <- EM.v2(model, tol, nj, ur.maxiter, r.interval, kfconv, parallel, num.cores)
      newpars <- sout$pars
      calls.v1 <- calls.v1 + sout$calls.v1
      iter <- iter + sout$iter
    } else {
      newpars <- step.EM.v1(model, kfconv)
      iter <- iter + 1
    }

    Mpars <- rbind(Mpars, newpars)
    model <- set.pars(model, newpars, FALSE)

    cond <- sqrt(sum((newpars - oldpars)^2)) > tol
  }

  list(Mpars = Mpars, pars = newpars, iter = iter, calls.v1 = calls.v1)
}
