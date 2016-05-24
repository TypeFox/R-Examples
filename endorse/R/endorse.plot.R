endorse.plot <- function(Y,
                         data,
                         scale,
                         dk = 98,
                         ra = 99,
                         yaxis = NULL,
                         col.seq = NA
                         ) {

  v.names <- paste("data$", Y, sep = "")
  plot.data <- as.data.frame(matrix(NA, nrow = nrow(data), ncol = length(v.names)))
  
  for (i in 1:length(v.names)) {
    plot.data[, i] <- eval(parse(text = paste("factor(", v.names[i],
                                   ", levels = c(1:",
                                   scale, ", ", dk, ",", ra, "))", sep = "")))
  }

  plot.table <- matrix(prop.table(table(plot.data[, 1])[c(scale:1, scale + 1, scale +2)]),
                       nrow = 1)
  for (i in 2:length(v.names)) {
    plot.table <- rbind(plot.table,
                        matrix(prop.table(table(plot.data[, i])[c(scale:1, scale + 1, scale +2)]),
                               nrow = 1))
  }

  if (is.na(col.seq[1])) {
    col.seq <- c(grey(seq(0.1, 0.9, length = scale)), rep(rgb(0,0,0,0), 2))
  }

  barplot(t(plot.table),
          horiz = TRUE, xaxt = "n", names.arg = yaxis,
          density = c(rep(-1, 5), rep(25, 2)),
          angle = c(rep(-45, 6), 45), lwd = 1.25, col = rgb(.31, .31, .31)) 

  barplot(t(plot.table),
          horiz = TRUE, xaxt = "n", names.arg = yaxis,
          density = c(rep(-1, 5), rep(25, 2)),
          angle = c(rep(-45, 6), 45), lwd = 1.25, col = col.seq, add = TRUE) 

}
