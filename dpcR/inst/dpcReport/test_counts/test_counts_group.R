dat <- coef(test_counts_dat())
dat[["run"]] <- as.factor(rownames(dat))
rownames(dat) <- NULL
dat
