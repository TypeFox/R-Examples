# Produces a plot illustrating teh results of the cLHS sampling.
# - a plot of the objective function
# - histograms or density plots of the sampled attributes
#
plot.cLHS_result <- function(
  x,
  modes = "obj",
  ...
  ){

  # Hack to avoid compilation error on ggplot args
  variable <- level <- value <- percent <- id <- NULL

  # Number of canvas to init
  n_views <- length(modes)
  iter <- 1:length(x$obj)
  obj <- x$obj
  cost <- x$cost
  if (is.null(cost)) df_obj <- data.frame(iter, obj)
  else df_obj <- data.frame(iter, obj, cost)

  pl <- list()

  ## Back-compatibility snippet code for ggplot2
  ver <- as.numeric_version(packageVersion('ggplot2'))

  # Objective function plot
  if ("obj" %in% modes) {

    ## Back-compatibility snippet code for ggplot2
    if (ver >= as.numeric_version("0.9.2")) {
      objective_plot <- ggplot(df_obj) + geom_line(aes(x = iter, y = obj)) + labs(title = "Evolution of the objective function", x = "Iteration", y = "Objective function") + theme_bw()
    } else {
      objective_plot <- ggplot(df_obj) + geom_line(aes(x = iter, y = obj)) + labs(x = "Iteration", y = "Objective function") + theme_bw() + ggtitle("Evolution of the objective function")
    }

    pl[[length(pl) + 1]] <- objective_plot
    names(pl)[length(pl)] <- 'obj'
  }

  # Cost function plot
  if ("cost" %in% modes) {

    ## Back-compatibility snippet code for ggplot2
    if (ver >= as.numeric_version("0.9.2")) {
      cost_plot <- ggplot(df_obj) + geom_line(aes(x = iter, y = cost)) + labs(title = "Evolution of the cost function", x = "Iteration", y = "Cost function") + theme_bw()
    } else {
      cost_plot <- ggplot(df_obj) + geom_line(aes(x = iter, y = cost)) + labs(x = "Iteration", y = "Cost function") + theme_bw()  + ggtitle("Evolution of the cost function")
    }
    
    pl[[length(pl) + 1]] <- cost_plot
    names(pl)[length(pl)] <- 'cost'
  }

  # Histogram/density plot
  if (any(c("hist", "dens", "box") %in% modes)) {
    
    if (length(which(c('den', 'box', 'hist') %in% modes)) > 1) stop('"hist", "dens", and "box" modes are mutually exclusive.')

    init <- x$initial_object
    spl <- x$sampled_data
    
    if (.is.spatial(init)) {
      if (inherits(init, "Spatial")) {
        init <- init@data
        spl <- spl@data
      }
      if (inherits(init, "Raster")) {
        init <- as.data.frame(init)
        spl <- as.data.frame(spl)
      }
    }

    # Separate continuous from factor variables
    i_factor <- which(!sapply(init, is.numeric))
    i_continuous <- setdiff(1:ncol(init), i_factor)
    n_factor <- length(i_factor)
    n_continuous <- length(i_continuous)

    init_continuous <- init[, i_continuous, drop = FALSE]
    spl_continuous <- spl[, i_continuous, drop = FALSE]
    init_factor <- init[, i_factor, drop = FALSE]
    spl_factor <- spl[, i_factor, drop = FALSE]

    # initiate an "id" column
    idcolname <- .create_unique_colname(id = "id", nm = names(init))
    
    if (n_factor > 0) {
      init_factor[[idcolname]] <- "init"
      spl_factor[[idcolname]] <- "spl"
      
      # merge df
      df_hist_factor <- melt(rbind(init_factor, spl_factor), idcolname)
      
      vars <- unique(df_hist_factor$variable)
      nvars <- length(vars)
      lvs <- dlply(df_hist_factor, "variable", function(x) unique(x$value))
      lst_prop_table <- lapply(1:nvars, function(x) {
        cur_var <- vars[x]
        cur_df <- df_hist_factor[which(df_hist_factor$variable == cur_var), , drop = FALSE]
        res <- ddply(cur_df, "id", function(y) {
          vect_vals <- factor(y$value, levels = lvs[[x]])
          prop.table(table(vect_vals))
        })
        res
      })
      names(lst_prop_table) <- vars
      
      lst_prop_table_melt <- lapply(1:length(lst_prop_table), function(x) data.frame(
            variable = names(lst_prop_table)[x], 
            melt(lst_prop_table[[x]], id = "id")
          )
        )

      df_prop_table <- do.call("rbind", lst_prop_table_melt)
      names(df_prop_table)[3] <- 'level'
 
      # Plot for factors (bar counts)
      ## Back-compatibility snippet code for ggplot2
      if (ver >= as.numeric_version("0.9.2")) {
        distrib_factor <- ggplot(df_prop_table) + 
          geom_point(aes(x = level, y = value, colour = id)) + 
          facet_wrap(~ variable, scales = "free") + 
          theme_bw() + labs(title = "Discrete variables") 
      } else {
        distrib_factor <- ggplot(df_prop_table) + 
          geom_point(aes(x = level, y = value, colour = id)) + 
          facet_wrap(~ variable, scales = "free") + 
          theme_bw() + ggtitle("Discrete variables") 
      }
            
      pl[[length(pl) + 1]] <- distrib_factor + 
        scale_y_continuous(name = "Relative  Frequency", labels = percent) + scale_x_discrete(name = "Level") +
        scale_colour_discrete(name = "")

      names(pl)[length(pl)] <- 'dens_factor'
    }
    
    if (n_continuous > 0) {
      init_continuous[[idcolname]] <- "init"
      spl_continuous[[idcolname]] <- "spl"
  
      # merge df
      df_hist_continuous <- melt(rbind(init_continuous, spl_continuous), idcolname)
      
      # Plot for continuous
      if ('dens' %in% modes) {
        ## Back-compatibility snippet code for ggplot2
        if (ver >= as.numeric_version("0.9.2")) {
          distrib_continuous <- ggplot(df_hist_continuous) + 
            geom_density(aes_string(x = 'value', fill = idcolname), alpha = 0.5) + 
            facet_wrap( ~ variable, scales = "free") + 
            theme_bw()  + 
            labs(title = "Continuous variables") + 
            scale_fill_discrete(name = "") + 
            scale_x_continuous(name = "Value") + 
            scale_y_continuous(name = "Density")
        } else {
          distrib_continuous <- ggplot(df_hist_continuous) + 
            geom_density(aes_string(x = 'value', fill = idcolname), alpha = 0.5) + 
            facet_wrap( ~ variable, scales = "free") + 
            theme_bw()  + 
            ggtitle("Continuous variables") + 
            scale_fill_discrete(name = "") + 
            scale_x_continuous(name = "Value") + 
            scale_y_continuous(name = "Density")
        }
      }
      if ('hist' %in% modes) {
        ## Back-compatibility snippet code for ggplot2
        if (ver >= as.numeric_version("0.9.2")) {
          distrib_continuous <- ggplot(df_hist_continuous) + 
            geom_histogram(aes_string(x = 'value', fill = idcolname), position = 'dodge') + 
            facet_wrap(~ variable, scales = "free") + 
            theme_bw() + 
            labs(title = "Continuous variables") + 
            scale_fill_discrete(name = "") + 
            scale_x_continuous(name = "Value") + 
            scale_y_continuous(name = "Count")
        } else {
          distrib_continuous <- ggplot(df_hist_continuous) + 
            geom_histogram(aes_string(x = 'value', fill = idcolname), position = 'dodge') + 
            facet_wrap(~ variable, scales = "free") + 
            theme_bw() + 
            ggtitle("Continuous variables") + 
            scale_fill_discrete(name = "") + 
            scale_x_continuous(name = "Value") + 
            scale_y_continuous(name = "Count")
        }
      }
      if ('box' %in% modes) {
        ## Back-compatibility snippet code for ggplot2
        if (ver >= as.numeric_version("0.9.2")) {
          distrib_continuous <- ggplot(df_hist_continuous) + 
            geom_boxplot(aes_string(x = idcolname, y = 'value')) + 
            facet_wrap( ~ variable, scales = "free") + 
            theme_bw() + 
            labs(title = "Continuous variables") + 
            scale_x_discrete(name = "") + 
            scale_y_continuous(name = "Value")
        } else {
          distrib_continuous <- ggplot(df_hist_continuous) + 
            geom_boxplot(aes_string(x = idcolname, y = 'value')) + 
            facet_wrap( ~ variable, scales = "free") + 
            theme_bw() + 
            ggtitle("Continuous variables") + 
            scale_x_discrete(name = "") + 
            scale_y_continuous(name = "Value")
        }
      }

      pl[[length(pl) + 1]] <- distrib_continuous 

      names(pl)[length(pl)] <- 'dens_continuous'
    }
  }

  if (length(pl) > 1) {
    grid.newpage()
    vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

    # If there's density plots
    if ("dens_continuous" %in% names(pl) & "dens_factor" %in% names(pl)) {
      pushViewport(viewport(layout = grid.layout(2, length(pl) - 1)))
      
      k_row <- 1
      k_col <- 1

      for (i_pl in 1:length(pl)) {
        
        if (names(pl)[i_pl] == "dens_continuous" | names(pl)[i_pl] == "dens_factor") {

          print(pl[[i_pl]], vp = vplayout(k_row, k_col))
          k_row <- k_row + 1

        } else {
          print(pl[[i_pl]], vp = vplayout(1:2, k_col))
          k_col <- k_col + 1
        }
      }
    } else {
      pushViewport(viewport(layout = grid.layout(1, length(pl))))
      k_col <- 1

      for (i_pl in 1:length(pl)) {
        print(pl[[i_pl]], vp = vplayout(1, k_col))
        k_col <- k_col + 1
      }
    }
  }
  else {
    print(pl[[1]])
  }

}
