#new_dat <- input_dat()
res <- summary(new_dat, print = FALSE)[["summary"]]
res <- cbind(run = paste0(res[["experiment"]], ".", res[["replicate"]]), res)
colnames(res) <- c("Run", "Experiment name", "Replicate ID", "Assay", "Method", "&lambda;", "&lambda; (lower CI)",
                   "&lambda; (upper CI)", "m", "m (lower CI)", "m (upper CI)", "k", "n")
