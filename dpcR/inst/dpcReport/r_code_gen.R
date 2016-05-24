#Generate R code for the report

#separates biggers chunks of the code
separator <- "\n    \n    ###############"

setup_l <- c("# Load packages",
             'library(dpcR)',
             if(input[["data_summary_scatter_rep"]] || input[["data_summary_test_counts"]] || input[["plot_panel"]] || input[["poisson_distr"]])
               c('library(ggplot2) # ggplot2 library for nice plots',
                 "# Define theme for plots",
                 'cool_theme <- theme(plot.background=element_rect(fill = "transparent", colour = "transparent"), panel.grid.major = element_line(colour="lightgrey", linetype = "dashed"), panel.background = element_rect(fill = "white", colour = "black"), legend.background = element_rect(fill="NA"), legend.position = "bottom", axis.text = element_text(size = 14), axis.title.x = element_text(size=17, vjust = -0.1), axis.title.y = element_text(size = 17, vjust = 1), strip.text = element_text(size = 17, face = "bold"), strip.background = element_rect(fill = "#9ecae1", colour = "black"), legend.text = element_text(size=14), legend.title = element_text(size = 17), plot.title = element_text(size = 22), legend.key = element_rect(fill = "white", colour = "black", linetype = "dashed", size = 0.5))'))


#read data line
read_data_l <- c("# Read and adjust data", if(is.null(input[["input_file"]])) {
  'input_data <- six_panels'
} else {
  c("# The input file is assumed to be in the current R working directory",
    sub("file_name", input[["input_file"]][["name"]], 
        switch(input[["input_type"]],
               raw_adpcr = 'input_data <- read_dpcr("file_name", format = "raw", adpcr = TRUE)',
               raw_ddpcr = 'input_data <- read_dpcr("file_name", format = "raw", adpcr = FALSE)',
               QX100 = 'input_data <- read_dpcr("file_name", format = "QX100")',
               BioMark_det = 'input_data <- read_dpcr("file_name", format = "BioMark", detailed = TRUE)',
               BioMark_sum = 'input_data <- read_dpcr("file_name", format = "BioMark", detailed = FALSE)')))
})

change_replicate_l <- if(any(as.character(slot(input_dat(), "replicate")) != rep_names_new())) {
  paste0('slot(input_data, "replicate") <- factor(', 
         paste0(capture.output(dput(rep_names_new())), collapse = ""), ')')
} else {
  ""
}

change_exper_l <- if(any(as.character(slot(input_dat(), "exper")) != exp_names_new())) {
  paste0('slot(input_data, "exper") <- factor(', 
         paste0(capture.output(dput(exp_names_new())), collapse = ""), ')')
} else {
  ""
}

data_summary_table_lc <- if(input[["data_summary_table_rep"]]) {
  c(separator, '# Print only table from summary.dpcr function',
    'summary(input_data, print = FALSE)[["summary"]]')
} else {
  ""
}


data_summary_scatter_lc <- if(input[["data_summary_scatter_rep"]]) {
  c(separator, '# Prepare data for plots',
    'plot_data <- summary(input_data, print = FALSE)[["summary"]]',
    paste0('plot_data <- plot_data[plot_data[["method"]] == "', input[["CI_method"]], '", ]'),
    # Plot boxplot',
    paste0('ggplot(plot_data, aes(x = experiment, y = lambda, ymax = lambda.up, ymin = lambda.low)) + geom_point(size = 4, alpha = 0.6, lty = 2, colour = "blue") + cool_theme + geom_boxplot(outlier.colour = NA, fill = adjustcolor("lightgrey", alpha.f = 0.25), shape = 15) + ggtitle(paste0("Experiment boxplot\\nCI method: ", "', input[["CI_method"]], '")) + scale_x_discrete("Experiment name") + scale_y_continuous(expression(lambda))'),
    '# Add new column, unique for every experiment/replicate combination',
    'plot_data[["exprep"]] <- factor(paste0(plot_data[["experiment"]], "\\n", plot_data[["replicate"]]))',
    '\n    # Plot stripchart',
    paste0('ggplot(plot_data, aes(y = exprep, x = lambda, colour = experiment, ymin = exprep, ymax = exprep)) + geom_point(size = 4) + cool_theme + ggtitle(paste0("Experiment/replicate scatter chart\\nCI method: ", "', input[["CI_method"]], '")) + scale_y_discrete("Replicate id", labels = plot_data[["replicate"]]) + scale_x_continuous(expression(lambda)) + scale_color_discrete("Experiment name") + geom_errorbarh(aes(x = lambda, xmin = lambda.low, xmax = lambda.up), size = 1.2, heigth = nlevels(plot_data[["exprep"]])/160)'))
} else {
  ""
}

data_summary_test_counts_lc <- if(input[["data_summary_test_counts"]]) {
  c(separator, '# Compare individual runs',
    'test_res <- test_counts(input_data, model = "ratio")',
    '# Results of the test with significance stars', 
    'test_res',
    '# Summary of results', 
    'summary(test_res)',
    '# Coefficients of runs',
    'run_coefs <- coef(test_res)',
    '# Add "run" column to prepare data for plot',
    'run_coefs[["run"]] <- as.factor(rownames(run_coefs))',
    '\n    # Plot coefficients',
    'ggplot(run_coefs, aes(y = run, x = lambda, colour = experiment, label = group)) + geom_point(size = 4) + cool_theme + geom_text(aes(x = lambda.up, y = run), show_guide = FALSE, hjust = -0.25, vjust = 0) + ggtitle("Grouped experiments") + scale_y_discrete("Replicate id", labels = run_coefs[["replicate"]] ) + scale_x_continuous(expression(lambda)) + coord_cartesian(xlim = c(ifelse(min(run_coefs[["lambda.low"]]) > 0, min(run_coefs[["lambda.low"]]) * 0.9, min(run_coefs[["lambda.low"]]) * 1.1), ifelse(max(run_coefs[["lambda.up"]]) < 0,  max(run_coefs[["lambda.up"]]) * 0.9, max(run_coefs[["lambda.up"]]) * 1.1))) + scale_size_discrete(guide = FALSE, range = c(5, 7)) + scale_color_discrete("Experiment name") + geom_errorbarh(aes(x = lambda, xmin = lambda.low, xmax = lambda.up), size = 1.2, heigth = nlevels(run_coefs[["run"]])/160)')
} else {
  ""
}

plot_panel_lc <- if(input[["plot_panel"]]) {
  c(separator, '# Test all panels, use "panels_test[[1]]" to see the first result and so on',
    paste0('panels_test <- test_panel(input_data, nx = ', 
           ifelse(is.null(input[["nx"]]), 5, input[["nx"]]), ', ',
           'ny = ', ifelse(is.null(input[["ny"]]), 5, input[["ny"]]), ')'),
    '# Plot panels, use "panels_plot[[1]]" to see the first array and so on',
    'panels_plot <- lapply(adpcr2panel(input_data), function(single_array) { ggplot(calc_coordinates(single_array, half = "none")[["ggplot_coords"]], aes(x = col, y = row , fill = as.factor(value))) + geom_tile(colour = "black", linetype = 2) + cool_theme  + scale_x_discrete("Column") + scale_y_discrete("Row") + scale_fill_discrete("Value") + theme(panel.border = element_blank(), panel.background = element_blank())})')
} else {
  ""
}

poisson_distr_lc <- if(input[["poisson_distr"]]) {
  plot_line <- 'lapply(1L:length(dens), function(i) ggplot(dens[[i]], aes(x = x, y = y)) + geom_line(colour = "lightskyblue1", size = 1.2) + geom_area(aes(fill = conf_up)) + geom_area(aes(fill = conf_low)) + scale_fill_manual(values = c("FALSE" = NA, "TRUE" = adjustcolor("cyan4", alpha.f = 0.5)), guide = FALSE) + cool_theme +  scale_y_continuous("Density") + ggtitle(names(dens)[i])'
  plot_line <- if(input[["density_plot_avg"]]) {
    paste0(plot_line, ' + scale_x_continuous(expression(lambda))')
  } else {
    paste0(plot_line, ' + scale_x_continuous("k")')
  }
  
  if(input[["density_plot_bars"]])
    plot_line <- paste0(plot_line, 
                        '+ geom_bar(stat = "identity", fill = adjustcolor("lightskyblue1", alpha.f = 0.5))')
  plot_line <- paste0(plot_line, ')')
  
  c(separator, '# Compute moments for all runs',
    'moments(input_data)',
    '# Plot distribution for all runs',
    paste0('dens <- dpcr_density_table(input_data, average = ', input[["density_plot_avg"]], ', ',
           'methods = "', input[["density_plot_methods"]], '", ',
           'conf.level = ', input[["density_plot_cil"]], ')'),
    plot_line)
} else {
  ""
}

all_lines <- c("\n    ", setup_l,
               read_data_l, change_replicate_l, change_exper_l, #input
               data_summary_table_lc, #summary table
               data_summary_scatter_lc, #summary plots
               data_summary_test_counts_lc,
               plot_panel_lc,
               poisson_distr_lc) 
all_lines <- all_lines[all_lines != ""]