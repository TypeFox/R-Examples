mom_dat <- moments(single_run)[, -c(1L:3)]
mom_tab <- cbind(mom_dat[1L:4, ], mom_dat[5L:8, 2])
colnames(mom_tab) <- c("Moment", "Theoretical", "Empirical")
