#helper function returning y values of confidence intervals
y_val_conf <- function(conf, data, side) {
  side_id <- ifelse(side == "left", 2, 3)
  id <- which(sort(c(data[ ,1], conf[[side_id]])) == conf[[side_id]])
  if (id != 1 && id <= nrow(data)) {
    y_l <- data[id - 1, 2]
    y_r <- data[id, 2]
    x_l <- data[id - 1, 1]
    x_r <- data[id, 1]
    #calculate y value for conf
    c(id + side_id - 3, (y_r - y_l)/(x_r - x_l)*(conf[[side_id]] - x_l) + y_l,
      conf[[side_id]])
  } 
}