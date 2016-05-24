object <- test_counts_dat()
signif_stars <- symnum(slot(object, "test_res")[, "p_value"], corr = FALSE, na = FALSE, 
                       cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                       symbols = c("***", "**", "*", ".", " "))
res <- data.frame(runs = rownames(slot(object, "test_res")), 
                  slot(object, "test_res"), 
                  signif = as.vector(signif_stars))
colnames(res) <- c("Compared pair of runs", "p-value", "Significance")
