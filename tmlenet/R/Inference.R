#----------------------------------------------------------------------------------
# Estimate IC-based Variance (as. Var) and CIs (based on f_W and fY)
# New fast method for as var calculation (matrix vs)
#----------------------------------------------------------------------------------
# use estimates of fWi (hold all W's fixed at once),
# loop over all intersecting friends networks
# calculate R_ij*(fW_i*fW_j) - see page 33 Network paper vdL
#----------------------------------------------------------------------------------
# unified estimator naming used throughout the package:
# c("TMLE", "h_IPTW", "MLE")
#----------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# THIS FUNCTION IS NOT USED, ONLY HERE FOR TESTING SPARSE CONNECTIVITY MATRIX CONVERSION / CREATION
get_sparse_Fiintersectmtx <- function() {
  # ------------------------------------------------------------------------------
  # Simulate some network data
  # ------------------------------------------------------------------------------
  n <- 20000
  # n <- 10000
  # n <- 1000
  # fvec_i <- abs(rnorm(n))
  fvec_i <- rnorm(n)
  # fvec_i <- rep_len(1,n)
  # NetInd_k <- matrix(NA, nrow=n, ncol=3)
  # NetInd_k[1,] <- c(2,3,NA)
  # NetInd_k[4,] <- c(2,NA,NA)
  # nF <- rep.int(0,n)
  # nF[1] <- 2; nF[4] <- 1
  # Kmax <- 150
  Kmax <- 200
  # Kmax <- 15
  NetInd_k <- t(replicate(n, sample(1:n, Kmax, replace = FALSE)))
  nF <- rep.int(Kmax, n)
  # ------------------------------------------------------------------------------

  # ------------------------------------------------------------------------------
  # Using sparse matrix implementation to get conn_ind_mtx_1st_indir
  # ------------------------------------------------------------------------------
  # estvartime <- system.time({
    # NetInd as sparse adjacency matrix (new version returns pattern sparse mat nngCMatrix):
    sparse_mat <- NetInd.to.sparseAdjMat(NetInd_k, nF = nF, add_diag = TRUE)
    # Second pass over columns of connectivity mtx to connect indirect intersections (i and j have a common friend but are not friends themselves):
    sparse_mat_ind <- Matrix::crossprod(sparse_mat) # t(sparse_mat)%*%sparse_mat returns nsCMatrix (only non-zero entries)
    # MOVED TO est.sigma_sparse():
    Dstar_crossprod_new <- sum_crossprod_Fij(sparseAdjMat = sparse_mat_ind, fvec_i = fvec_i)
    Dstar_new <- (1/n) * ((2*Dstar_crossprod_new) + sum(fvec_i^2))
  # })
  # print("estvartime"); print(estvartime)

  # print(object.size(NetInd_k), units="Mb")
  # 3906.4 Kb # 0 Gb
  # sparse_mat as dgCMatrix:
  # print(object.size(sparse_mat), units="Mb")
  # 12,032.1 Kb
  # system.time({
  #   conn_ind_mtx_1st_indir <- get.Fiintersectmtx(n = n)$conn_ind_mtx_1st
  #   Dstar_old <- est.sigma_fsum(get.crossprodmtx(fvec_i), conn_ind_mtx_1st_indir)        
  # })
  # print(all.equal(Dstar_new, Dstar_old))
  # # print(all.equal(Dstar, Dstar_old)); print(Dstar);
  # print(Dstar_new); print(Dstar_old)
  return(list(conn_ind_mtx_1st = sparse_mat_ind))
}

# sparse matrix class from Matrix package:
# dgCMatrix
  # class is a class of sparse numeric matrices in the compressed, sparse, column-oriented format.
# nsparseMatrix-classes
# ?'nsCMatrix-class'
# ngCMatrix, nsCMatrix, and ntCMatrix
  # class of sparse pattern matrices, i.e., binary matrices conceptually with TRUE/FALSE entries. Only the positions of the elements that are TRUE are stored.
  # Objects can be created by calls of the form new("ngCMatrix", ...) and so on. 
  # More frequently objects are created by coercion of a numeric sparse matrix to the pattern form for use in the symbolic 
  # analysis phase of an algorithm involving sparse matrices. Such algorithms often involve two phases: a symbolic phase wherein the 
  # positions of the non-zeros in the result are determined and a numeric phase wherein the actual results are calculated. 
  # During the symbolic phase only the positions of the non-zero elements in any operands are of interest, hence numeric sparse matrices can be treated as sparse pattern matrices.
# lsparseMatrix-classes {Matrix}
# lsCMatrix
# The second letter in the name of these non-virtual classes indicates general, symmetric, or triangular.
# The lsparseMatrix class is a virtual class of sparse matrices with TRUE/FALSE or NA entries. Only the positions of the elements that are TRUE are stored.

# Copied from simcausal, then modified to work as sparse pattern matrix:
# return pattern sparse matrix, no @x is recorded (ngMatrix):
NetInd.to.sparseAdjMat <- function(NetInd_k, nF, add_diag = FALSE) {
  nobs <- nrow(NetInd_k)
  sdims <- c(nobs, nobs)
  nnonzero <- sum(!is.na(NetInd_k))
  sparse_p <- as.integer(c(0, cumsum(nF)))
  sparse_iwNA <- as.vector(t(NetInd_k))
  sparse_i <- sparse_iwNA[!is.na(sparse_iwNA)] - 1
  out_sparseAdjMat <-  Matrix::sparseMatrix(i = sparse_i, p = sparse_p, dims = sdims, index1 = FALSE, giveCsparse = TRUE)
  if (add_diag) diag(out_sparseAdjMat) <- TRUE
  return(out_sparseAdjMat)
}

# For correlated i,j (s.t. i<j), add a term f_vec[i]*f_vec[j] to the cumulative sum
# Helper function to calculate cross product sum of correlated f_Wi (see p.33 of vdL)
sum_crossprod_Fij <- function(sparseAdjMat, fvec_i) {
  # sparseAdjMat:
    # nnzero: number of non-zero elements
    # i: These are the 0-based row numbers for each non-zero element in the matrix.
    # "integer" of length nnzero, 0- based row numbers for each non-zero element in the matrix, i.e., i must be in 0:(nrow(.)-1).
    # p: integer vector for providing pointers, one for each column, to the initial (zero-based) index of elements in the column.
    # .@p is of length ncol(.) + 1, with p[1] == 0 and
    # p[length(p)] == nnzero, such that in fact, diff(.@p) are the number of non-zero elements for each column.  
  assertthat::assert_that(is(sparseAdjMat, "sparseMatrix"))
  # 1) The number of friends for each observation:
  nF <- as.integer(diff(sparseAdjMat@p))
  # 2) Column based cummulative number of non-zero entries (cummulative nF)
  cumFindx <- sparseAdjMat@p
  # 3) All non-zero elements as a vector of 0-based row numbers:
  base0_IDrownums <- sparseAdjMat@i
  # 4) For each observation i that has non-zero nF (friends), add fvec_i[i]*fvec_Fj for each friend Fj of i:  
  non0nF.idx <- which(nF > 1L) # don't care if nF[i]=1 since it means i has 0 actual friends (i itself is included in nF)
  # non0nF.idx <- which(nF > 0L)
  Dstar_crossprod <- 0
  for (idx in non0nF.idx) {
    Fidx_base0 <- (cumFindx[idx]) : (cumFindx[idx + 1] - 1)
    FriendIDs <- base0_IDrownums[Fidx_base0 + 1] + 1
    # remove the diag term (idx,idx) from FriendIDs (always the last entry),
    # since we only need to sum all fvec_i[i]^2 once (outside this fun)
    FriendIDs <- FriendIDs[-length(FriendIDs)]
    Dstar_crossprod <- Dstar_crossprod + sum(fvec_i[idx] * fvec_i[FriendIDs])
  }
  return(Dstar_crossprod)
}

# Sum the cross prod vector over connectivity mtx (prod will only appear if (i,j) entry is 1):
est.sigma_sparse <- function(fvec_i, sparse_connectmtx)  {
  n <- length(fvec_i)
  # sum of fvec_i[i]*fvec[j] for correlated cross-product terms (i,j) s.t. i<j
  Dstar_crossprod <- sum_crossprod_Fij(sparseAdjMat = sparse_connectmtx, fvec_i = fvec_i)
  # double cross prod sum + sum of squares over i=1,...,n
  Dstar <- (1/n) * ((2*Dstar_crossprod) + sum(fvec_i^2))
  return(Dstar)
}

est_sigmas <- function(estnames, n, NetInd_k, nF, obsYvals, ests_mat, QY_mat, wts_mat, fWi_mat) {
# est_sigmas <- function(n, NetInd_k, nF, obsYvals, ests_mat, QY_mat, wts_mat, fWi_mat, onlyTMLE_B) {
  fWi <- fWi_mat[, "fWi_Qinit"]
  QY.init <- QY_mat[, "QY.init"] 
  h_wts <- wts_mat[, "h_wts"]
  var_tmle <- 0

  # NetInd as sparse adjacency matrix (new version returns pattern sparse mat ngCMatrix):
  sparse_mat <- NetInd.to.sparseAdjMat(NetInd_k, nF = nF, add_diag = TRUE)
  # Second pass over columns of connectivity mtx to connect indirect intersections (i and j have a common friend but are not friends themselves):
  connectmtx_1stO <- Matrix::crossprod(sparse_mat) # t(sparse_mat)%*%sparse_mat returns nsCMatrix (only non-zero entries)

  # TMLE inference based on the iid IC:
  iidIC_tmle <- h_wts * (obsYvals - QY.init) + (fWi - ests_mat[rownames(ests_mat)%in%"TMLE",])
  var_tmle <- est.sigma_sparse(iidIC_tmle, connectmtx_1stO)
  # Simple estimator of the iid asymptotic IC-based variance (no adjustment made when two observations i!=j are dependent):
  var_iid.tmle <- mean((iidIC_tmle)^2)
  # IPTW h (based on the mixture density clever covariate (h)):
  iidIC_iptw_h <- h_wts * (obsYvals) - (ests_mat[rownames(ests_mat)%in%"h_IPTW",])
  var_iptw_h <- est.sigma_sparse(iidIC_iptw_h, connectmtx_1stO)

  var.ests <- c(abs(var_tmle), abs(var_iptw_h), NA)
  as.var_mat <- matrix(0, nrow = length(var.ests), ncol = 1)
  as.var_mat[,1] <- var.ests
  rownames(as.var_mat) <- estnames
  colnames(as.var_mat) <- "Var"

  # QY.star <- QY_mat[, "QY.star"]
  # g_wts <- wts_mat[,"g_wts"]  
  # var_tmle_A <- var_tmleiptw_1stO <- var_tmleiptw_2ndO <- var_iptw_1stO <- var_iptw_2ndO <- 0
  # var_tmle_A_Q.init <- var_tmle_B_Q.init <- 0
  # # TMLE A (clever covariate update): Inference based on the iid IC analogy, QY.init := initial Q model predictions, h_wts := h_tilde
  # if (!onlyTMLE_B) {
  #   iidIC_tmle_A <- h_wts * (obsYvals - QY.init) + (fWi - ests_mat[rownames(ests_mat)%in%"tmle_A",])
  #   var_tmle_A <- est.sigma_sparse(iidIC_tmle_A, connectmtx_1stO)
  # }
  # # TMLE B (weighted model update): Inference based on the iid IC:
  # iidIC_tmle_B <- h_wts * (obsYvals - QY.init) + (fWi - ests_mat[rownames(ests_mat)%in%"tmle_B",])
  # var_tmle_B <- est.sigma_sparse(iidIC_tmle_B, connectmtx_1stO)
  # # simple iid estimator of the asymptotic variance (no adjustment made when two observations i!=j are dependent):
  # var_iid.tmle_B <- mean((iidIC_tmle_B)^2)
  # TMLE based on iptw clever covariate (more non-parametric):
  # if (!onlyTMLE_B) {
  #   iidIC_tmleiptw <- g_wts * (obsYvals - QY.init) + (fWi - ests_mat[rownames(ests_mat)%in%"tmle_g_iptw",])
  #   var_tmleiptw_1stO <- est.sigma_sparse(iidIC_tmleiptw, connectmtx_1stO)
  # }
  # # IPTW g:
  # if (!onlyTMLE_B) {
  #   iidIC_iptw_g <- g_wts * (obsYvals) - (ests_mat[rownames(ests_mat)%in%"g_iptw",])
  #   var_iptw_1stO <- est.sigma_sparse(iidIC_iptw_g, connectmtx_1stO)
  # }
  # Inference based on the EIC, with factorization into orthogonal components sigma2_DY and sigma2_W_N
  # sigma2_DY_i are independent (since they are conditioned on W,A)
  # sigma2_W_N_i are dependent => need to take double sum of their crossprod among dependent units
  # if (!onlyTMLE_B) {
  #   D_star_Yi.Qinit <- h_wts * (obsYvals - QY.init) # h*(Y-Q_bar_N):
  #   sigma2_DY <- (1/n) * sum(D_star_Yi.Qinit^2)  # Sum_{i} (D_star_Yi)^2

  #   # fW_A_crossprod <- get.crossprodmtx((fWi - ests_mat[rownames(ests_mat)%in%"tmle_A",]))
  #   # sigma2_W_N_A <- est.sigma_fsum(fW_A_crossprod, connectmtx_1stO)
  #   fW_A_i <- fWi - ests_mat[rownames(ests_mat)%in%"tmle_A",]
  #   sigma2_W_N_A <- est.sigma_sparse(fW_A_i, connectmtx_1stO)
  #   var_tmle_A_Q.init <- sigma2_W_N_A + sigma2_DY

  #   # **NEW** TMLE B (weights model update)
  #   # fW_B_crossprod <- get.crossprodmtx((fWi - ests_mat[rownames(ests_mat)%in%"tmle_B",]))
  #   # sigma2_W_N_B <- est.sigma_fsum(fW_B_crossprod, connectmtx_1stO)
  #   fW_B_i <- fWi - ests_mat[rownames(ests_mat)%in%"tmle_B",]
  #   sigma2_W_N_B <- est.sigma_sparse(fW_B_i, connectmtx_1stO)
  #   var_tmle_B_Q.init <- sigma2_W_N_B + sigma2_DY

  #   # D_star_Yi.Qstar <- h_wts * (obsYvals - QY.star)
  #   # D_star_Yi.Qstar[determ.Q] <- 0
  #   # fDY_crossprod <- get.crossprodmtx(D_star_Yi.Qstar)
  #   # double sum over dependent subjects, Sum_{i,j} R_W(i,j)*D_star_Yi*D_star_Yj
  #   # sigma2_Y_N <- est.sigma_fsum(fDY_crossprod, connectmtx_1stO)
  #   # sigma2_Y_N <- est.sigma_sparse(D_star_Yi.Qstar, connectmtx_1stO)
  #   # var_tmle_Q.init_c <- sigma2_W_N_A + sigma2_Y_N

  #   # #--------
  #   # # conservative estimate of the as. variance from EIC for TMLE A:
  #   # # abs terms double sum over dependent subjects, Sum_{i,j} R_W(i,j)*|D_star_Yi|*|D_star_Yj|:
  #   # fabsDY_crossprod <- get.crossprodmtx(abs(D_star_Yi.Qstar))
  #   # abs_sigma2_Y_N <- est.sigma_fsum(fabsDY_crossprod, connectmtx_1stO)
  #   # abs_sigma2_Y_N <- est.sigma_sparse(abs(D_star_Yi.Qstar), connectmtx_1stO)
  #   # var_tmle_A_Q.star_cons <- sigma2_W_N_A + abs_sigma2_Y_N
  #   # # --------
  # }

  # var.ests <- c(abs(var_tmle_A), abs(var_tmle_B), abs(var_tmleiptw_1stO), abs(var_iptw_h), abs(var_iptw_1stO), 0)
  # estnames <- c( "TMLE_A", "TMLE_B", "TMLE_g_IPTW", "h_IPTW", "g_IPTW", "MLE")
  # other.vars = c(
  #               var_iid.tmle_B = abs(var_iid.tmle_B), # no adjustment for correlations i,j
  #               var_tmleiptw_2ndO = abs(var_tmleiptw_2ndO), # adjusting for 2nd order dependence of i,j
  #               var_iptw_2ndO = abs(var_iptw_2ndO), # adjusting for 2nd order dependence of i,j
  #               var_tmle_A_Q.init = abs(var_tmle_A_Q.init), # using the EIC & Q.init for TMLE A
  #               var_tmle_B_Q.init = abs(var_tmle_B_Q.init)  # using the EIC & Q.init for TMLE B
  #               )

  return(list(as.var_mat = as.var_mat))
  # return(list(as.var_mat = as.var_mat, other.vars = other.vars))
}

# create output object with param ests of EY_gstar, vars and CIs for given gstar (or ATE if two tmle obj are passed)
make_EYg_obj <- function(estnames, estoutnames, alpha, DatNet.ObsP0, tmle_g_out, tmle_g2_out=NULL) {
  nobs <- DatNet.ObsP0$nobs
  NetInd_k <- DatNet.ObsP0$netind_cl$NetInd_k
  nF <- DatNet.ObsP0$netind_cl$nF
  ests_mat <- tmle_g_out$ests_mat
  QY_mat <- tmle_g_out$QY_mat
  fWi_mat <- tmle_g_out$fWi_mat
  wts_mat <- tmle_g_out$wts_mat
  if (!is.null(tmle_g2_out)) {
    ests_mat <- tmle_g_out$ests_mat - tmle_g2_out$ests_mat
    fWi_mat <- tmle_g_out$fWi_mat - tmle_g2_out$fWi_mat
    wts_mat <- tmle_g_out$wts_mat - tmle_g2_out$wts_mat
  }

  # get the iid IC-based asymptotic variance estimates:
  getVar_time <- system.time(
    as.vars_obj <- est_sigmas(estnames = estnames, n = nobs, NetInd_k = NetInd_k, nF = nF, 
                              obsYvals = DatNet.ObsP0$noNA.Ynodevals,
                              ests_mat = ests_mat, QY_mat = QY_mat, 
                              wts_mat = wts_mat, fWi_mat = fWi_mat)
  )
  if (gvars$verbose) {
    print("time to estimate Vars: "); print(getVar_time)  
  }
  
  get_CI <- function(xrow, n) {
    f_est_CI <- function(n, psi, sigma2_N) { # get CI
      z_alpha <- qnorm(1-alpha/2)
      CI_est <- c(psi - z_alpha*sqrt(sigma2_N) / sqrt(n), psi + z_alpha*sqrt(sigma2_N) / sqrt(n))
      return(CI_est)
    }
    psi <- xrow["estimate"];
    sigma2_N <- xrow["Var"];
    return(f_est_CI(n = n, psi = psi, sigma2_N = sigma2_N))
  }

  CIs_mat <- t(apply(cbind(ests_mat, as.vars_obj$as.var_mat), 1, get_CI, n = nobs))
  colnames(CIs_mat) <- c("LBCI_"%+%as.character(alpha/2), "UBCI_"%+%as.character(1-alpha/2))

  # ------------------------------------------------------------------------------------------
  # RENAME ESTIMATORS FOR THE FINAL OUTPUT:
  # ------------------------------------------------------------------------------------------
  # print("ests_mat: "); print(ests_mat)
  # print("ests_vars: "); print(as.vars_obj$as.var_mat)
  # print("ests_CIs: "); print(CIs_mat)
  rownames(ests_mat) <- estoutnames
  rownames(as.vars_obj$as.var_mat) <- estoutnames
  rownames(CIs_mat) <- estoutnames

  EY_g.star <- list(estimates = ests_mat,
                    vars = (as.vars_obj$as.var_mat / nobs),
                    CIs = CIs_mat,
                    other.vars = (as.vars_obj$other.vars / nobs),
                    # other.vars = lapply(as.vars_obj$other.vars, function(var) var / nobs)
                    h_g0_SummariesModel = NULL,
                    h_gstar_SummariesModel = NULL)

  if (is.null(tmle_g2_out)) {
    EY_g.star[["h_g0_SummariesModel"]] <- tmle_g_out$h_g0_SummariesModel
    EY_g.star[["h_gstar_SummariesModel"]] <- tmle_g_out$h_gstar_SummariesModel
  }

  return(EY_g.star)
}