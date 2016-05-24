
setClassUnion("OptionalNumeric", c("numeric", "NULL"))
setClassUnion("OptionalInteger", c("integer", "NULL"))
setClassUnion("OptionalCharacterORListORFunction", c("character", "list", "function", "NULL"))
setClassUnion("OptionalMatrix", c("matrix", "NULL"))
setClassUnion("OptionalMatrixNumeric", c("matrix", "numeric", "NULL"))

setClass(Class = "stsm", 
representation = representation(
  call = "language",
  model = "character",
  y = "ts",
  diffy = "ts",
  xreg = "OptionalMatrixNumeric",
  fdiff = "function",
  ss = "list",
  pars = "numeric",
  nopars = "OptionalNumeric",
  cpar = "OptionalNumeric",
  lower = "numeric",
  upper = "numeric",
  transPars = "OptionalCharacterORListORFunction",
  ssd = "OptionalNumeric", #sample spectral density (periodogram)
  sgfc = "OptionalMatrix") # constant terms in the spectral generating function
)

stsm.model <- function(model = c("local-level", "local-trend", "BSM", 
  "llm+seas", "trend+ar2"),
  y, pars = NULL, nopars = NULL, cpar = NULL, xreg = NULL,
  lower = NULL, upper = NULL, transPars = NULL, 
  ssd = FALSE, sgfc = FALSE)
{
  if (missing(y) || !is.ts(y))
    stop("A time series object must be provided in argument 'y'.")
  
  if (any(duplicated(c(names(pars), names(nopars), names(cpar)))))
    stop("duplicated parameters passed through 'pars', 'nopars' and 'cpar'")

  if (!is.null(xreg))
  {
    if (NROW(xreg) != length(y))
      stop("lengths of ", sQuote("y"), " and ", sQuote("xreg"), " do not match.")
  }

  model <- match.arg(model)
  type <- "var"
  y1 <- if (!is.na(y[1])) y[1] else mean(y, na.rm = TRUE)
  vy <- var(y, na.rm = TRUE) / 100
  vyl <- 1e6 * vy

  switch(model,

    "local-level" = 
    {
      diffy <- diff(y)
##FIXME see a way to get "s" from the object that calls the function
      fdiff <- function(x, s) diff(x)
      Z <- rbind(1)
      mT <- rbind(1)
      H <- paste(type, 1, sep = "")
      Q <- matrix(paste(type, 2, sep = ""), nrow = 1, ncol = 1)
      R <- rbind(1)
      V <- matrix(paste(type, 2, sep = ""), nrow = 1, ncol = 1)
      a0 <- cbind("a01")
      P0 <- matrix("P01", nrow = 1, ncol = 1)

      pars0 <- c("var1" = 1, "var2" = 1)
      nopars0 <- c("a01" = y1, "P01" = vyl)
    },

    "local-trend" = 
    {
      diffy <- diff(diff(y))
      fdiff <- function(x, s) diff(diff(x))
      Z <- rbind(c(1, 0))
      mT <- rbind(c(1, 1), c(0, 1))
      H <- paste(type, 1, sep = "")
      Q <- diag(2)
      diag(Q) <- paste(type, 2:3, sep = "")
      R <- diag(2)
#http://r.789695.n4.nabble.com/Behaviors-of-diag-with-character-vector-in-R-3-0-0-td4663735.html
      #V <- diag(paste(type, 2:3, sep = ""))
      V <- diag(2)
      diag(V) <- paste(type, 2:3, sep = "")
      a0 <- paste("a0", seq(1, 2), sep = "")
      P0 <- diag(2)
      diag(P0) <- paste("P0", seq(1, 2), sep = "")

      pars0 <- c("var1" = 1, "var2" = 1, "var3" = 1)
      nopars0 <- c("a01" = y1, "a02" = 0, "P01" = vyl, "P02" = vyl)
    },

    "BSM" = 
    {
      s <- frequency(y)
      diffy <- diff(diff(y, s))
      fdiff <- function(x, s) diff(diff(x, s))
      Z <- rbind(c(1, 0, 1, rep(0, s - 2)))
      mT <- cbind(
        cbind(c(1, rep(0, s)), c(1, 1, rep(0, s - 1))),
        rbind(0, 0, -1, diag(s - 2)), c(0, 0, -1, rep(0, s - 2)))
      H <- paste(type, 1, sep = "")
      #Q <- diag(c(paste(type, 2:4, sep = ""), rep(0, s - 2)))
      Q <- diag(s + 1)
      diag(Q) <- c(paste(type, 2:4, sep = ""), rep(0, s - 2))
      R <- rbind(diag(3), matrix(0, nrow = s - 2, ncol = 3))
      #V <- diag(c(paste(type, 2:4, sep = "")))
      V <- diag(3)
      diag(V) <- c(paste(type, 2:4, sep = ""))
      a0 <- paste("a0", seq(1, s + 1), sep = "")
      #P0 <- diag(paste("P0", seq(1, s + 1), sep = ""))
      P0 <- diag(s + 1)
      diag(P0) <- paste("P0", seq(1, s + 1), sep = "")

      pars0 <- c("var1" = 1, "var2" = 1, "var3" = 1, "var4" = 1)
      nopars0 <- c(rep(0, s), rep(vyl, s + 1))
      names(nopars0) <- c(paste("a0", seq(2, s + 1), sep = ""), 
        paste("P0", seq(1, s + 1), sep = ""))
      nopars0 <- c("a01" = y1, nopars0)
    },

    "llm+seas" = 
    {
      s <- frequency(y)
      diffy <- diff(y, s)
      fdiff <- function(x, s) diff(x, s)
      Z <- rbind(c(1, 1, rep(0, s - 2)))
      mT <- cbind(c(1, rep(0, s-1)),
        rbind(0, rep(-1, s-1), cbind(diag(s-2), 0)))

      H <- paste(type, 1, sep = "")
      #Q <- diag(c(paste(type, 2:3, sep = ""), rep(0, s - 2)))
      Q <- diag(s)
      diag(Q) <- c(paste(type, 2:3, sep = ""), rep(0, s - 2))
      R <- rbind(diag(2), matrix(0, nrow = s - 2, ncol = 2))
      #V <- diag(c(paste(type, 2:3, sep = "")))
      V <- diag(2)
      diag(V) <- c(paste(type, 2:3, sep = ""))
      a0 <- paste("a0", seq(1, s), sep = "")
      #P0 <- diag(paste("P0", seq(1, s), sep = ""))
      P0 <- diag(s)
      diag(P0) <- paste("P0", seq(1, s), sep = "")

      pars0 <- c("var1" = 1, "var2" = 1, "var3" = 1)
      nopars0 <- c(rep(0, s-1), rep(vyl, s))
      names(nopars0) <- c(paste("a0", seq(2, s), sep = ""), 
        paste("P0", seq(1, s), sep = ""))
      nopars0 <- c("a01" = y1, nopars0)
    },

    "trend+ar2" = # Clark's:87 model 
    {
      diffy <- diff(diff(y))
      fdiff <- function(x, s) diff(diff(x))
      p <- 2
      Z <- rbind(c(1, 0, 1, rep(0, p - 1)))
      mT <- rbind(c(1, 1, rep(0, p)), c(0, 1, rep(0, p)),
        c(0, 0, paste("phi", seq(1, p), sep = "")), 
        cbind(0, 0, diag(p - 1), 0))
      H <- paste(type, 1, sep = "")
      #Q <- diag(c(paste(type, 2:4, sep = ""), 0))
      Q <- diag(4)
      diag(Q) <- c(paste(type, 2:4, sep = ""), 0)
      R <- rbind(diag(3), 0)
      #V <- diag(c(paste(type, 2:4, sep = "")))
      V <- diag(3)
      diag(V) <- c(paste(type, 2:4, sep = ""))
      a0 <- paste("a0", seq(1, 2 + p), sep = "")
      #P0 <- diag(paste("P0", seq(1, 2 + p), sep = ""))
      P0 <- diag(2 + p)
      diag(P0) <- paste("P0", seq(1, 2 + p), sep = "")

      pars0 <- c("var1" = 1, "var2" = 1, "var3" = 1, "var4" = 1,
        "phi1" = 0.5, "phi2" = 0.5)
      nopars0 <- c("a01" = y1, "a02" = 0, "a03" = 0, "a04" = 0, 
        "P01" = vy, "P02" = vy, "P03" = vy, "P04" = vy)
    }
  )

  ss <- list(Z = Z, T = mT, H = H, R = R, V = V, Q = Q, a0 = a0, P0 = P0)

  allpars0 <- c(pars0, nopars0)

  if (!is.null(pars))
  {
    # na.omit necessary in the presence of "xreg" cofficients
    ref <- match(names(pars), names(allpars0))
    id <- which(is.na(ref))
    lid <- length(id) 
    if (lid > 0) {
      xregcoefs <- pars[id]
      pars <- pars[-id]
      if (length(pars) == 0)
        pars <- NULL
    } else
      xregcoefs <- NULL
    ref <- na.omit(ref)
    if (length(ref) > 0)
      allpars0 <- allpars0[-ref]
  } else 
    xregcoefs <- NULL

  # do not use "else" because it depends on the previous "if" statement
  #else { 
  if (is.null(pars))
  {
    ref <- match(names(nopars), names(pars0))
    if (length(ref) > 0) {
      pars <- pars0[-ref]
    } else 
      pars <- pars0
    ref <- match(names(pars), names(allpars0))
    if (length(ref) > 0)
      allpars0 <- allpars0[-ref]
  }

  if (!is.null(nopars))
  {
    ref <- match(names(nopars), names(allpars0))
    if (length(ref) > 0)
      allpars0[ref] <- nopars
  } #else
    #nopars <- allpars0
    #nopars <- nopars0
  nopars <- allpars0

  if (!is.null(cpar))
  {
    nms <- names(cpar)
    nmspars <- names(pars)
    nmsnopars <- names(nopars)
    if (nms %in% nmspars){
      pars <- pars[-which(nmspars == nms)]
    } else if (nms %in% nmsnopars){
      nopars <- nopars[-which(nmsnopars == nms)]
    } else stop("bug e1") #'nms' should have been assigned to 'pars' or 'nopars'
  }

  ub <- rep(Inf, length(pars))
  names(ub) <- names(pars)
  lb <- -ub
  if (!is.null(lower)) {
    lb[names(lower)] <- lower
  } else if (is.null(transPars) || (!is.function(transPars) && transPars == "StructTS"))
    lb[which(substr(names(pars), 1, nchar(type)) == type)] <- 0
  if (!is.null(upper))
    ub[names(upper)] <- upper

  # sample spectral density (periodogram)

  if (ssd)
  {
    ssd <- as.vector(Mod(fft(diffy))^2 / (2 * pi * length(diffy)))
  } else ssd <- NULL

  # create "stsm" object

  x <- new("stsm", call = match.call(), model = model, y = y, 
    diffy = diffy, xreg = NULL, fdiff = fdiff,
    ss = ss, pars = pars, nopars = nopars, cpar = cpar, 
    lower = lb, upper = ub,
    transPars = transPars, ssd = ssd, sgfc = NULL)

  # constant terms in the spectral generating function

  if (sgfc)
    x@sgfc <- stsm.sgf(x, FALSE, FALSE, FALSE)$constants

  # external regressors

  if (!is.null(xreg))
    x <- set.xreg(x, xreg, xregcoefs)

  x
}
