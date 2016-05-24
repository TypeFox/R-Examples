#### iterative_filteration ####

.iterative_filteration <- function(models = list(),
                                   fisher_combined_res,
                                   res_by_marker_set,
                                   eQTL_marker_set,
                                   T.i=5, T.e=10,
                                   reference_data,
                                   expression_data,
                                   genotyping_data) {

  num_of_marker_sets <- length(models)

  num_of_iter <- 0
  num_of_removed_genes <- -1

  eQTL_results_above_thresh <- unique(colnames(eQTL_marker_set)[which(eQTL_marker_set>=T.e, arr.ind=T)])

  marker_info <- list()

  while (num_of_removed_genes != 0) {
    print(paste("num of iteration:", num_of_iter, sep=" "))

    num_of_removed_genes <- 0

    print("applying Refinement")
    max_of_each_row <- apply(fisher_combined_res, 1, max)
    cond <- sapply(max_of_each_row, function(x) x >= T.i)

    places <- apply(fisher_combined_res, 1, which.max)
    iQTL_results_above_thresh <- unique(unlist(places[cond]))

    eQTL_results_above_thresh_combined_with_iQTL <- eQTL_results_above_thresh[
      as.numeric(eQTL_results_above_thresh) %in% iQTL_results_above_thresh]

    genes_connected_to_eQTL <- unique(rownames(
      which(eQTL_marker_set[, as.numeric(eQTL_results_above_thresh_combined_with_iQTL)]>=T.e, arr.ind=T)))

    print("Refining marker sets space")

    for (i in c(1:num_of_marker_sets)) {

      marker_set_name <- colnames(models[[i]][1])
      marker_set_data <- models[[i]][!models[[i]][,1] %in% marker_info[[marker_set_name]],]
      print(paste("working on marker_set", marker_set_name, sep=" "))

      new_marker_set_data <- marker_set_data[!(marker_set_data %in% genes_connected_to_eQTL)]

      marker_info[[marker_set_name]] <- c(marker_info[[marker_set_name]],
                                          marker_set_data[!(marker_set_data %in% new_marker_set_data)])

      print(paste("we removed" , length(marker_set_data[!(marker_set_data %in% new_marker_set_data)]),
                  "genes from", marker_set_name, ":", paste(marker_set_data[!(marker_set_data %in% new_marker_set_data)],
                                                            collapse = " "), sep=" "))

      num_of_removed_genes <- num_of_removed_genes +
        length(marker_set_data[!(marker_set_data %in% new_marker_set_data)])

      if (length(new_marker_set_data) != length(marker_set_data)) {

        new_marker_set_data = as.data.frame(new_marker_set_data)
        names(new_marker_set_data) = marker_set_name

        models[[i]] = new_marker_set_data

        res_by_marker_set[[marker_set_name]] <- .create_iQTL_association_scores(reference_data = reference_data,
                                                                             mix_data = expression_data,
                                                                             marker_set = models[[i]][1],
                                                                             genotype_data = genotyping_data)
      }
    }

    if (num_of_removed_genes == 0) {
      break
    } else {
      fisher_combined_res <- .combine_iQTL_association_scores(res_by_marker_set)
    }
    num_of_iter <- num_of_iter + 1
  }

  final_res <- list("final_association_score" = fisher_combined_res, "marker_info" = marker_info)
  return(final_res)

}

