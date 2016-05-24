kn <- unlist(summary(single_run, print = FALSE)[["summary"]][1, c("k", "n")])

conf <- dpcr_density(k = kn["k"], n = kn["n"], average = input[["density_plot_avg"]], 
                     methods = input[["density_plot_methods"]], 
                     conf.level = input[["density_plot_cil"]], plot = FALSE)

dens <- data.frame(dpcR:::dpcr_calculator(kn["k"], kn["n"], average = input[["density_plot_avg"]]))
colnames(dens) <- c("x", "y")

dens[["conf_low"]] <- dens[["x"]] <= conf[["lower"]] 
dens[["conf_up"]] <- dens[["x"]] >= conf[["upper"]]