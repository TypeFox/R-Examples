scale_weights <- function(weights,scale) {
  if(is.null(scale)) {
    return(weights)
  } else if(length(unique(weights)) == 1) {
    weights <- rep(mean(scale), length(weights))
    return(weights)
  }  else {
    #Scale the values so that they lie in the interval (-1,1)
    weights <- (weights-mean(weights)) / (max(weights)-min(weights)) 
    #Then scale them into the region given by the scale parameter
    #This parameter should be a vector of length 2 with the upper and lower
    #boundaries of the scale
    weights <- (weights * (max(scale)-min(scale)) + mean(scale))
    return(weights)
  }
}
from_coords_to_plot <- function(x,y,domain,number_of_dimensions) {
  #normalize x and y so that they are in [0,1]
  xnorm <- (x-domain$left) / (domain$right-domain$left)
  #The y coordinates are "special" in that we want to show the first plot at the top,
  #but the coordinates start at the bottom, so we have to "invert" y coordinate first
  yinv <- domain$top - y
  ynorm <- (yinv-domain$bottom) / (domain$top-domain$bottom)
  x_plot_number <- as.integer(xnorm / (1/number_of_dimensions)) + 1
  y_plot_number <- as.integer(ynorm / (1/number_of_dimensions)) + 1
  return(c(x_plot_number,y_plot_number))
}
dataframe_row_to_html <- function(row) {
  if(is.null(row)) {
    return("")
  }
  value_strings <- sapply(names(row),function(name){return(paste0(name,": ",row[[name]],"<br>"))})
  return(paste0(value_strings,collapse=""))
}
#Problem: the ggplot hover output only gives us the x and y value of the data 
#point that we are hovering over and does not tell us whether we are hovering 
#over an actual data point or just over some empty space in the plot. We would 
#like to have the full dataframe row of the closest point to the cursor if there
#is a point within some area around the cursor. Solution: Take out the correct x
#and y column from the plot. Calculate the distances w.r.t. that point. Find out
#the smallest value, then return the corresponding row if the value is smaller
#than area_around_cursor.
row_from_two_values <- function(dataframe,hover_list,area_around_cursor) {
  if(is.null(hover_list) | is.null(dataframe) | is.null(hover_list$mapping)) {
    return(NULL)
  }
  dims <- c(hover_list$mapping[[1]],hover_list$mapping[[2]])
  dim_vals <- c(hover_list$x,hover_list$y)
  dataframe_columns <- dataframe[dims]
  #Now we need to normalize the x and y axis values. The reasoning for this is
  #that if we compute distances without normalizing and the space of possible values
  #for the x and y axis is very different (e.g. x is between 0 and 1, y is between 1 and 10000),
  #the axes are not taken into account equally, which they really should.
  #To do so, we first subtract the smallest values of the data spaces for x and y
  dim_vals <- c(dim_vals[1]-hover_list$domain$left,dim_vals[2]-hover_list$domain$bottom)
  dataframe_columns[dims[1]] <- dataframe_columns[dims[1]] - hover_list$domain$left
  dataframe_columns[dims[2]] <- dataframe_columns[dims[2]] - hover_list$domain$bottom
  #Then we divide by the width of the x and the y dataspace
  x_width <- hover_list$domain$right - hover_list$domain$left
  y_width <- hover_list$domain$top - hover_list$domain$bottom
  dim_vals <- c(dim_vals[1]/x_width,dim_vals[2]/y_width)
  dataframe_columns[dims[1]] <- dataframe_columns[dims[1]] / x_width
  dataframe_columns[dims[2]] <- dataframe_columns[dims[2]] / y_width
  
  best_index <- which.min(fields::rdist(dataframe_columns,matrix(dim_vals,nrow=1)))
  best_value <- fields::rdist(dataframe_columns[best_index,],matrix(dim_vals,nrow=1))
  if(best_value > area_around_cursor) {
    return(NULL)
  }else {
    return(dataframe[best_index,])
  }
}