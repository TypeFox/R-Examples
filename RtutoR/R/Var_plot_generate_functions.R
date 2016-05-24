applicable_plotlist <- function(x,y){

  plots <- vector()

  if (is.null(y)) {
    if (is.numeric(x)){
      plots <- c("geom_density","geom_histogram")
    } else {
      plots <- c("geom_bar")
    }

  } else if (is.factor(x) & is.numeric(y)) {
    plots <- c("geom_boxplot","geom_violin")
  } else if (is.numeric(x) & is.numeric(y)) {
    plots <- c("geom_point","geom_smooth","geom_line","geom_area","geom_text")
  } else {
    plots <- c("geom_bar")
  }
  plots

}


######### Generate Variable Dropdown ###########################333

generateVarDropdown <- function(selected_dataset,number_layers) {
  selected_dataset <- eval(parse(text = selected_dataset))
  vars <- colnames(selected_dataset)
  layerCount <- as.integer(number_layers)

  # Dynamically generate the variable list
  lapply(1:layerCount,function(i) {
    div(
      id = "var_control",
      selectInput(
        inputId = paste("xaxis_var",i),paste("Layer ", i,  "- Select x-axis"),choices = c("",vars)
      ),
      selectInput(
        inputId = paste("yaxis_var",i),paste("Layer ", i, "- Select y-axis"),choices = c("",vars)
      )
    )
  })

}

############# Generate plot list ##############################333

generatePlotList <- function(selected_dataset,x_var,y_var,show_all_plots_flag,layer_count) {

    plot_list <- applicable_plotlist(selected_dataset[[x_var]],selected_dataset[[y_var]])

    if (show_all_plots_flag == F) {
      selectInput(inputId = paste("plot_type",layer_count),paste("Layer ", layer_count, " - Select Plot"), choices = plot_list)
    } else {
      selectInput(inputId = paste("plot_type",layer_count),paste("Layer ", layer_count, " - Select Plot"), choices = all_plots)
    }

}




