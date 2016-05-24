new_dat <- change_data(input_dat(), as.factor(rep_names_new()), as.factor(exp_names_new()))
roi <- extract_dpcr(new_dat, input[["array_choice"]])

summs <- summary(roi, print = FALSE)[["summary"]]
summs <- cbind(region = rep("Whole array", nrow(summs)), summs)

if(!is.null(array_val[["selected"]])) {
  slot(roi, ".Data") <- slot(roi, ".Data")[array_val[["selected"]], , drop = FALSE]
  slot(roi, "n") <- sum(array_val[["selected"]])
  summs <- rbind(summs, cbind(region = rep("Selected region", nrow(summs)), 
                              summary(roi, print = FALSE)[["summary"]]))
}

colnames(summs) <- c("Region", "Experiment name", "Replicate ID", "Assay", "Method", "&lambda;", 
                     "&lambda; (lower CI)", "&lambda; (upper CI)", "m", 
                     "m (lower CI)", "m (upper CI)", "k", "n")
summs
