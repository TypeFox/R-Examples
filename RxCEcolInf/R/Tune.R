`Tune` <-  function(fstring, data=NULL, num.runs=12, num.iters=10000,
                    rho.vec=rep(0.05, ntables),
                    kappa=10, nu=(mu.dim+6), psi=mu.dim,
                    mu.vec.0=rep(log((.45/(mu.dim-1))/.55), mu.dim),
                    mu.vec.cu=runif(mu.dim, -3, 0),
                    nolocalmode=50,
                    sr.probs=NULL, sr.reps=NULL, 
                    numscans=1, Diri=100, dof=4, debug=1)
{

  data.mats <- get.mat(fstring, data)
  NNtots <- cbind(data.mats$rowmat, data.mats$colmat)
  nrow.pt <- ncol(data.mats$rowmat)
  ncol.pt <- ncol(data.mats$colmat)
  mu.dim <- (ncol.pt - 1) * nrow.pt
  ntables <- nrow(data.mats$colmat)
  
  # Check there are at least 2 rows/cols:
  if (nrow.pt<2 || ncol.pt<2)
    stop("Table must be at least 2x2")

  # Check the NNtots is valid:
  if (any(apply(NNtots[,1:nrow.pt],1,sum) != apply(NNtots[,(nrow.pt+1):(nrow.pt+ncol.pt)],1,sum)))
    stop("Row and column totals do not sum to the same total total")

  NNbounds <- gq.bounds(nrow.pt,ncol.pt,NNtots)
  storage.mode(NNtots) <- storage.mode(NNbounds) <- "double"
  if (any(is.null(sr.probs)))
    sr.probs <- matrix(1/nrow.pt,nrow=nrow(NNtots),ncol=nrow.pt)
  storage.mode(sr.probs) <- "double"
  if (any(is.null(sr.reps)))
    sr.reps <- matrix(1,nrow=nrow(NNtots),ncol=nrow.pt)
  storage.mode(sr.reps) <- "integer"

  # Error-check Diri:
  Diri <- rep(as.numeric(Diri),nrow(NNtots))
  if (length(Diri)!=nrow(NNtots))
    stop("Diri should be a scalar quantity")
  if (any(Diri<0.0))
    stop("Diri should be non-negative")
  storage.mode(Diri) <- "double"

  args <- list("NNtots" = as.matrix(NNtots),
               "NNbounds" = as.matrix(NNbounds),
               "numrows_pt" = as.integer(nrow.pt),
               "numcols_pt" = as.integer(ncol.pt),
               "sr_probs" = as.matrix(sr.probs),
               "sr_reps" = as.matrix(sr.reps),
               "num_iters" = as.integer(num.iters),
               "psi_0" = as.numeric(psi),
               "kappa_0" = as.numeric(kappa),
               "nu_0" = as.numeric(nu),
               "nolocalmode" = as.numeric(nolocalmode),
               "num_runs" = as.integer(num.runs),
               "dof" = as.numeric(dof),
               "num_runs" = as.integer(num.runs),
               "mu_vec_0" = as.numeric(mu.vec.0),
               "mu_vec_cu" = as.numeric(mu.vec.cu),
               "use_Diri_every_vec" = as.numeric(Diri),
               "rho_vec" = as.numeric(rho.vec),
               "numscans" = as.integer(numscans),
               "dbg" = as.integer(debug))
  .Call("Tune", PACKAGE="RxCEcolInf", args)
}
