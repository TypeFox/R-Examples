#' Computing simulations results
#'
#'@examples
#'
#'\dontrun{
#'data_sims <- data_sim_voomlike(maxGSsize=300, minGSsize=30)
#'res <- compute_sim_voomlike(counts = data_sims$counts,
#'                            design = data_sims$design,
#'                            gs_keep = data_sims$gs_keep,
#'                            indiv = data_sims$indiv,
#'                            alternative=FALSE,
#'                            RE_indiv_sd=0.1)
#'res_all <- cbind(res$res_voom, res$res_perso, res$res_noweights)
#'colnames(res_all) <- paste0(rep(c("asym", "perm", "camera"), 3),
#'                            rep(c("_voom", "_perso", "_noweights"), each=3))
#'}
#'
#'@keywords internal
#'@importFrom stats rnorm
#'@export
compute_sim_voomlike <- function(counts, design, gs_keep, indiv, alternative=FALSE,
                                 fixed_eff = 0.5, fixed_eff_sd = 0.2,
                                 rand_eff_sd = 0.25, RE_indiv_sd=NULL, eps_sd=0.05){

  logcoutspm <- apply(counts, MARGIN=2, function(v){log2((v + 0.5)/(sum(v) + 1)*10^6)})

  if(alternative){
    n <- ncol(counts)
    nb_g <- nrow(counts)
    beta_time <- stats::rnorm(1, mean=fixed_eff, sd=fixed_eff_sd)
    gamma_genes_time <- matrix(stats::rnorm(nb_g, mean=0, sd=rand_eff_sd), ncol=1)
    logcoutspm_alt <- logcoutspm + (gamma_genes_time + beta_time)%*%design[, "time"]

    # for(gs in gs_keep){
    #   nb_g_gs <- length(gs)
    #   gamma_genes_time <- matrix(stats::rnorm(nb_g_gs, mean = 0, sd = rand_eff_sd), ncol = 1)
    #   logcoutspm_alt[gs,] <- logcoutspm[gs,] + (gamma_genes_time+beta_time)%*%design[, "time"]
    # }

    err <- matrix(stats::rnorm(n*nb_g, mean = 0, sd = eps_sd), nrow = nb_g, ncol = n)
    logcoutspm_alt <- logcoutspm_alt + err

  }else{
    logcoutspm_alt <- logcoutspm
  }

  #Individual Random effects
  if(!is.null(RE_indiv_sd)){
    RE_indiv <- stats::rnorm(n=length(unique(indiv)), mean=0, sd=RE_indiv_sd)
    logcoutspm_alt <- logcoutspm_alt + RE_indiv[as.numeric(as.character(indiv))]
  }


  # Weights estimation ----
  #########################

  #png(width=5.5, height=4.5, units="in", file="voom_ex.png", res=300)
  #w_voom <- voom_weights(x=design, y=t(counts2), doPlot = TRUE, preprocessed = FALSE)
  #dev.off()
  w_voom <- voom_weights(x=design, y=logcoutspm_alt, doPlot = FALSE, preprocessed = TRUE)

  w_perso <- sp_weights(x=design[,-which(colnames(design)=="time")],
                        y=t(logcoutspm_alt),
                        phi = cbind(design[, "time"]),
                        doPlot = FALSE,
                        exact=FALSE,
                        preprocessed = TRUE
  )

  # Testing ----
  ##############
  res_voom <- NULL
  res_perso <- NULL
  res_noweights <- NULL
  gs_ind_list <- limma::ids2indices(gs_keep, identifier=rownames(logcoutspm_alt))
  cor_limma <- limma::duplicateCorrelation(logcoutspm_alt, design, block = indiv)

  for(gs in 1:length(gs_keep)){
    gs_test <- gs_keep[[gs]]
    y_test <- logcoutspm_alt[gs_test,]
    x_test <- design[,-which(colnames(design)=="time")]
    w_test_perso <- w_perso[as.numeric(gs_test), ]
    w_test_voom <- w_voom[as.numeric(gs_test), ]
    phi_test <- cbind(design[, "time"])

    ## Fitting the null model
    n <- ncol(y_test)
    g <- length(gs_test)
    y_T_vect <- as.vector(y_test)
    indiv_vect <- rep(indiv, g)
    g_vect <- as.factor(rep(1:g, n))
    x_vect <-do.call(rbind, replicate(g, x_test, simplify=FALSE))
    phi_vect <- do.call(rbind, replicate(g, phi_test, simplify=FALSE))

    res_voom <- rbind(res_voom, cbind("asym" = vc_test_asym(y = y_test, x=x_test, indiv=indiv, phi=phi_test,
                                                            Sigma_xi = as.matrix(diag(ncol(phi_test))),
                                                            w = w_test_voom)[["pval"]],
                                      "permut" = vc_test_perm(y = y_test, x=x_test, indiv=indiv, phi=phi_test,
                                                              Sigma_xi = as.matrix(diag(ncol(phi_test))),
                                                              w = w_test_voom)[["pval"]],
                                      "camera" = limma::camera(y=logcoutspm_alt, index=gs_ind_list[[gs]],
                                                               design=design, contrast=3,
                                                               weights = w_voom)$PValue,
                                      "roast" = limma::roast(y=logcoutspm_alt, index=gs_ind_list[[gs]],
                                                             design=design, contrast=3,
                                                             block = indiv, correlation = cor_limma$consensus.correlation,
                                                             weights = w_voom)$p.value["Mixed", "P.Value"])
    )
    res_perso <- rbind(res_perso, cbind("asym" = vc_test_asym(y = y_test, x=x_test, indiv=indiv, phi=phi_test,
                                                              Sigma_xi = as.matrix(diag(ncol(phi_test))),
                                                              w = w_test_perso)[["pval"]],
                                        "permut" = vc_test_perm(y = y_test, x=x_test, indiv=indiv, phi=phi_test,
                                                                Sigma_xi = as.matrix(diag(ncol(phi_test))),
                                                                w = w_test_perso)[["pval"]],
                                        "camera" = limma::camera(y=logcoutspm_alt, index=gs_ind_list[[gs]],
                                                                 design=design, contrast=3,
                                                                 weights = w_perso)$PValue,
                                        "roast" = limma::roast(y=logcoutspm_alt, index=gs_ind_list[[gs]],
                                                               design=design, contrast=3,
                                                               block = indiv, correlation = cor_limma$consensus.correlation,
                                                               weights = w_perso)$p.value["Mixed", "P.Value"])
    )
    res_noweights <- rbind(res_noweights, cbind("asym" = vc_test_asym(y = y_test, x=x_test, indiv=indiv, phi=phi_test,
                                                                      Sigma_xi = as.matrix(diag(ncol(phi_test))),
                                                                      w = rep(1, length(gs_test)))[["pval"]],
                                                "permut" = vc_test_perm(y = y_test, x=x_test, indiv=indiv, phi=phi_test,
                                                                        Sigma_xi = as.matrix(diag(ncol(phi_test))),
                                                                        w = rep(1, length(gs_test)))[["pval"]],
                                                "camera" = limma::camera(y=logcoutspm_alt, index=gs_ind_list[[gs]],
                                                                         design=design, contrast=3)$PValue,
                                                "roast" = limma::roast(y=logcoutspm_alt, index=gs_ind_list[[gs]],
                                                                       block = indiv, correlation = cor_limma$consensus.correlation,
                                                                       design=design, contrast=3)$p.value["Mixed", "P.Value"])
    )

  }

  return(list("res_voom"=res_voom, "res_perso"=res_perso, "res_noweights"=res_noweights))
}
