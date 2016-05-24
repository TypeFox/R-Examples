CommT.kfdist <-
function (post_gt_distrs_BEAST,
                         post_gt_distrs_starBEAST,
                         outlier_num = 1,
                         treedist_select = 2) {
  # Descr:  calculate the tree distances
  # Deps:   phangorn::treedist
  #         reshape::melt
  # I/p:    post_gt_distrs_BEAST = object of class phylo; multiple posterior distributions of gene trees from BEAST
  #         post_gt_distrs_starBEAST = object of class phylo; multiple posterior distributions of gene trees from starBEAST
  #         outlier_num
  #         treedist_select = number of the tree distance selected (2 = Kuhner-Felsenstein distance)
  # O/p:    matrix of floats

  n_genes = length(post_gt_distrs_BEAST)
  # 1. Check if valid input
    if (n_genes > 1) {
  # 2. Initialize outmatrix
        out_m = matrix(nrow=length(post_gt_distrs_BEAST[[1]]),
                       ncol=n_genes,
                       NA)
  # 3. Calculate tree distance for each tree pair
        for (h in 1:n_genes) {
            cat("\n\nANALYZING LOCUS '", names(post_gt_distrs_BEAST)[h], "'\n", sep="")
            cat("\nComparing tree pair:\n ")
            for (i in 1:length(post_gt_distrs_BEAST[[h]])) {
                cat(i, " ", sep="")
                out_m[i,h] = phangorn::treedist(post_gt_distrs_BEAST[[h]][[i]],
                                       post_gt_distrs_starBEAST[[h]][[i]],
                                       check.labels = TRUE)[treedist_select]
            }
        }

  # 4. Rename outmatrix
        clmnames = lapply(sprintf("%03d", 1:n_genes), function(x){paste("gene", x, sep="")})
        colnames(out_m) = clmnames

    } else {stop(cat("\n\nERROR: List of posterior distributions of gene trees from BEAST:  not greater than 1.\n"))}

  # 5. Melt outmatrix
    out_df = reshape::melt(out_m, id.vars=c(clmnames))

  # 6. Add grouping variable
    out_df[,ncol(out_df)+1] = rep("regular", length(out_df[,3]))
    outlier_id = paste("gene", sprintf("%03d", as.integer(outlier_num)), sep="")
    out_df[which(out_df[,2]==outlier_id),4] = "outlier"
    colnames(out_df) = c('post_gen', 'gene_id', 'KF_dist', 'grouping_var')
    
  # 7. Return out_df
    return(out_df)
}
