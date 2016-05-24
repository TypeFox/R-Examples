# Find the active set of items for each observation
find_active <- function(modfr) {
  Zt <- modfr$reTrms$Zt
  Lambdat <- modfr$reTrms$Lambdat
  LambdatZt <- Lambdat %*% Zt
  LambdatZt <- as(LambdatZt, "dgTMatrix")
  tapply(LambdatZt@i + 1, LambdatZt@j + 1, c, simplify = FALSE)
}

add_obs_C <- function(C, act, obs_rem) {
  obs_C <- obs_rem[vapply(act[obs_rem],
                          function(items){all(is.element(items, C))},
                          TRUE)]
  list(C = C, obs_C = obs_C)
}

add_obs <- function(cliques, act, obs) {
  obs_rem <- obs
  cliques_ext <- list()
  for(i in seq_along(cliques)) {
    cliques_ext[[i]] <- add_obs_C(cliques[[i]], act, obs_rem)
    obs_rem <- setdiff(obs_rem, cliques_ext[[i]]$obs)
  }
  if(length(obs_rem) != 0) {
    stop("Couldn't allocate some observations to a clique")
  }
  cliques_ext
}

find_local_term <- function(clique_ext, modfr) {
  C <- clique_ext$C
  obs_C <- clique_ext$obs_C
  X_C <- unname(modfr$X[obs_C, , drop = FALSE])
  Zt_C <- unname(as.matrix(modfr$reTrms$Zt[C, obs_C, drop = FALSE]))
  Lambdat_C <- modfr$reTrms$Lambdat[C, C, drop = FALSE]
  Lind_mat <- modfr$reTrms$Lambdat
  Lind_mat@x <- as.numeric(modfr$reTrms$Lind)
  Lind_C <- as.integer(Lind_mat[C, C, drop = FALSE]@x) - 1
  resp <- model.response(modfr$fr)
  if(is.factor(resp)) {
    resp <- (resp != levels(resp)[1L])
  }
  weights <- model.weights(modfr$fr)
  if(is.null(weights)) {
    weights_C <- rep(1, length(obs_C))
  } else {
    weights_C <- weights[obs_C]
  }
  if(NCOL(resp) == 1) {
    resp_C <- resp[obs_C]
    resp_C <- matrix(unname(resp_C), ncol = 1)
  } else {
    #####################
    # FIXME: should convert to vector form
    resp_C <- resp[obs_C, , drop = FALSE]
    #####################
  }
  list(C = C, X = X_C, Zt = Zt_C, Lambdat = Lambdat_C,
       Lind = Lind_C, resp = resp_C, weights = weights_C)
}

find_factorization_terms <- function(modfr) {
  act <- find_active(modfr)
  cliques <- unname(unique(act))
  n_obs <- ncol(modfr$reTrms$Zt)
  n_re <- nrow(modfr$reTrms$Zt)
  cliques_ext <- add_obs(cliques, act, 1:n_obs)
  factorization_terms <- continuous_beliefs()
  for(i in seq_along(cliques_ext)) {
    lmodfr <- find_local_term(cliques_ext[[i]], modfr)
    factorization_terms$append_glmm_belief(items = lmodfr$C - 1,
                                           X = lmodfr$X,
                                           Zt = lmodfr$Zt,
                                           Lambdat = lmodfr$Lambdat,
                                           Lind = lmodfr$Lind,
                                           response = lmodfr$resp,
                                           weights = lmodfr$weights)
  }
  I1 <- matrix(1, nrow = 1, ncol = 1)
  for(i in seq_len(n_re)) {
    factorization_terms$append_normal_belief(items = c(i - 1),
                                             mean = c(0),
                                             precision = I1)
  }
  return(factorization_terms)
}

lmodfr_to_oneline <- function(lmodfr, file = "") {
  size_clique <- length(lmodfr$C)
  nobs <- nrow(lmodfr$X)
  nfixed <- ncol(lmodfr$X)
  Lambdat_triplet <- as(lmodfr$Lambdat, "dgTMatrix")
  nentries <- length(Lambdat_triplet@i)
  Lambdat_print <- as.numeric(rbind(Lambdat_triplet@i, Lambdat_triplet@j,
                                    Lambdat_triplet@x))
  ret <- c("G", size_clique, lmodfr$C - 1, nobs, nfixed, lmodfr$X, lmodfr$Zt,
           nentries, Lambdat_print,
           lmodfr$Lind, lmodfr$resp, lmodfr$weights)
  cat(cat(ret, file = file, append = TRUE),
      "\n", sep = "", file = file, append = TRUE)
}

save_factorization_to_file <- function(modfr, file = "") {
  unlink(file)

  act <- find_active(modfr)
  cliques <- unname(unique(act))
  n_obs <- ncol(modfr$reTrms$Zt)
  n_re <- nrow(modfr$reTrms$Zt)
  cliques_ext <- add_obs(cliques, act, 1:n_obs)

  # GLMM model frame terms
  for(i in seq_along(cliques_ext)) {
    lmodfr <- find_local_term(cliques_ext[[i]], modfr)
    lmodfr_to_oneline(lmodfr, file)
  }

  # normal terms
  for(i in seq_len(n_re)) {
    cat("N", 1, i - 1, 0, 1, "\n", sep = " ", file = file, append = TRUE)
  }
}

save_normal <- function(mean, precision, file = "") {
  unlink(file)

  precision <- as(precision, "dgTMatrix")
  nentries <- length(precision@i)
  precision_print <- c(nentries, as.numeric(rbind(precision@i,
                                                  precision@j,
                                                  precision@x)))
  cat(mean, "\n", file = file, append = TRUE)
  cat(precision_print, "\n", file = file, append = TRUE)
}
