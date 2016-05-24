
KalmanFilter <- function(y, ss, 
  KF.version = c("KFKSDS", "StructTS", "KFAS", "FKF", "dlm", "dse"),
  KF.args = list(), check.args = TRUE, debug = FALSE)
{
  KF.version <- match.arg(KF.version)[1]
  a <- if (check.args) {
    make.KF.args(ss, KF.version, KF.args) 
  } else KF.args

  P0cov <- if (is.null(KF.args$P0cov)) FALSE else KF.args$P0cov
  if (P0cov)
    ss$P0[] <- ss$P0[1]

  switch(KF.version,

    "KFKSDS" = 
    {
      res <- list(mloglik = 
        KF.C(y, ss, convergence = a$convergence, t0 = a$t0))
    },

    "StructTS" = 
    {
      mod <- list(Z = ss$Z, a = ss$a0, P = ss$P0, T = ss$T, 
        V = ss$Q, h = ss$H, Pn = ss$P0)
      res <- stats::KalmanRun(y, mod, -1, TRUE)
      res[[1]][1] <- 2 * res[[1]][1] - log(res[[1]][2])

      names(res) <- c("values", "residuals", "states")
      res$mloglik <- 0.5 * sum(res[[1]])  #0.5 * sum(res[[1]]) / length(y)
	
      res[[1]] <- NULL
      id <- match(diag(ss$V), diag(ss$Q))
      res[[2]] <- res[[2]][,id]
      attributes(res[[1]]) <- attributes(y)
      res[[2]] <- ts(res[[2]])
      tsp(res[[2]]) <- tsp(y)
    },

    "KFAS" = 
    {
      #require("KFAS")
    
      #mkfas <- KFAS::SSModel(y = y, Z = ss$Z, H = ss$H, T = ss$T, 
      #  R = ss$R, Q = ss$V, a1 = a$a0, P1 = ss$P0, P1inf = a$P1inf,
      #  distribution = "Gaussian", transform = "none", tolF = a$tolF)

      #ss <- char2numeric(m0, FALSE)
      #do not use "KFAS::" inside KFASS::SSmodel
      tmp <- KFAS::SSMcustom(ss$Z, ss$T, ss$R, ss$V, a$a0, ss$P0, a$P1inf, 1, length(y))
      mkfas <- KFAS::SSModel(y ~ -1 + tmp, H = ss$H, 
        data = data.frame(y = y), distribution = "gaussian", tol = a$tolF)

      res <- KFAS::KFS(mkfas, smoothing = "none", simplify = TRUE)
      res$mloglik <- -res$logLik
    },

    "FKF" = 
    {
      #require("FKF")

      res <- FKF::fkf(a0 = ss$a0, P0 = ss$P0, 
        dt = matrix(0, nrow = ncol(ss$T)), ct = a$ct, 
        Tt = ss$T, Zt = ss$Z, HHt = ss$Q, GGt = as.matrix(ss$H), 
        yt = rbind(y), check.input = a$check.input)
      res$mloglik <- -res$logLik
    },

    #"sspir" = 
    #{
    #  require("sspir")
    #
    #  yaux <- ts(matrix(y))
    #  tsp(yaux) <- tsp(y)
    #  msspir <- sspir::SS(y = yaux,
    #    Fmat = function(tt,x,phi) { return(matrix(phi$Z)) },
    #    Gmat = function(tt,x,phi) { return(phi$T) },
    #    Vmat = function(tt,x,phi) { return(matrix(phi$H)) },
    #    Wmat = function(tt,x,phi) { return(phi$Q) },
    #    m0 = rbind(ss$a0), C0 = ss$P0, phi = ss)
    #
    #  res <- sspir::kfilter(msspir)
    #  res$mloglik <- drop(-res$loglik)
    #},

    "dlm" = 
    {
      #require("dlm")

      mdlm <- list(m0 = ss$a0, C0 = ss$P0,
        FF = rbind(ss$Z), V = matrix(ss$H), GG = ss$T, W = ss$Q)

      res1 <- dlm::dlmFilter(y = as.vector(y), mod = mdlm, 
        debug = FALSE, simplify = TRUE)
      res2 <- dlm::dlmLL(y = as.vector(y), mod = mdlm, debug = FALSE)
      res <- c(res1, mloglik = res2)
    },

    "dse" = 
    {
      #require("dse")

      mdse <- dse::SS(F. = ss$T, G = NULL, H = ss$Z, K = NULL, Q = ss$Q, 
        R = matrix(ss$H), z0 = rbind(ss$a0), P0 = ss$P0, rootP0 = NULL,
        constants = NULL, description = NULL, names = NULL, 
        input.names = NULL, output.names = NULL)
      ##NOTE requires loading "dse"
      TSdata <- dse::TSdata
      res <- dse::l(mdse, TSdata(output = y), result = "like")
      res <- list(mloglik = -res)
    }
  )

  res
}

make.KF.args <- function(ss, KF.version, KF.args = list())
{
#NOTE only ss$a0 and ss$P0 are used in argument 'ss'
#keep the entire list 'ss' as argument for possible extensions
  check.KF.args <- function(list1, list2, label)
  {
    nms1 <- names(list1)
    nms2 <- names(list2)
    notinlist2 <- !(nms1 %in% nms2)
    wno <- which(notinlist2)

    if (any(notinlist2))
      warning("the following elements in argument 'KF.args'\n",
      "are omitted with 'KF.version = ", label, "': ", 
      paste(nms1[wno], collapse = ", "), ".", call. = FALSE)

    if (length(list1) > 0)
    {
      if (any(notinlist2)) {
        list2[nms1[which(!notinlist2)]] <- list1[-wno] 
      } else
        list2[nms1] <- list1
    }

    list2
  }

#  KF.version <- match.arg(KF.version, 
#    eval(formals(KFKSDS::KalmanFilter)$KF.version))

  P0cov <- if (is.null(KF.args$P0cov)) FALSE else KF.args$P0cov
  if (P0cov && !is.null(ss))
    ss$P0[] <- ss$P0[1]

  #NOTE 'inf' is currently not used
  ldef <- list(P0cov = FALSE, inf = 99999)

  ldef <- # list of default parameters for each interface
  switch(KF.version,

    "KFKSDS" = 
    {
      c(ldef, list(t0 = 1, convergence =  c(0, 9999))) #sUP = 0
    },

    "stats" = { ldef },

    "StructTS" = { ldef },

    "KFAS" = 
    {
      a0 <- ss$a0
      P1inf <- matrix(0, nrow(ss$P0), ncol(ss$P0))
      if (all(ss$P0 == 0))
      {
        a0[] <- 0
              diag(P1inf) <- 1  
      }
    
      c(ldef, list(#yt = y, Zt = ss$Z, Tt = ss$T, Rt = ss$R, 
        #Ht = ss$H, Qt = ss$V, a1 = ss$a0, P1 = ss$P0, 
        a0 = a0, P1inf = P1inf, tolF = 1e-08))
    },

    "FKF" = 
    {
      c(ldef, list(
        ct = matrix(0), 
        check.input = TRUE))
    },

    "sspir" = { ldef },

    "dlm" = { ldef },

    "dse" = { ldef }
  )

  args <- check.KF.args(KF.args, ldef, KF.version)

  args
}
