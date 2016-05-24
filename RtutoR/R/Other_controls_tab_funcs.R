
########### Axis slider function #############################33

genAxisSlider <- function(selected_dataset,plot_out) {

  y_var <- as.character(plot_out$labels[2])
  x_var <- as.character(plot_out$labels[1])

  if(is.numeric(selected_dataset[[y_var]])) {
    min_val_y = min(selected_dataset[[y_var]], na.rm = T)
    max_val_y = max(selected_dataset[[y_var]], na.rm = T)
  }

  if(is.numeric(selected_dataset[[x_var]])) {
    min_val_x = min(selected_dataset[[x_var]], na.rm = T)
    max_val_x = max(selected_dataset[[x_var]], na.rm = T)
  }


  if (is.numeric(selected_dataset[[x_var]]) & is.numeric(selected_dataset[[y_var]])) {
    div(id ="axis_range",
        sliderInput("y_slider","Change y axis range", min = min_val_y, max = max_val_y, value = c(min_val_y,max_val_y)),
        sliderInput("x_slider","Change x axis range", min = min_val_x, max = max_val_x, value = c(min_val_x,max_val_x))
    )
  } else if(!is.numeric(selected_dataset[[x_var]]) & is.numeric(selected_dataset[[y_var]])) {
    div(id="axis_range",
        sliderInput("y_slider","Change y axis range", min = min_val_y, max = max_val_y, value = c(min_val_y,max_val_y))
    )
  }

}

