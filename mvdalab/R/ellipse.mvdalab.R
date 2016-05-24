ellipse.mvdalab <- function (data, center = c(0, 0), radius = "chi", scale = TRUE, 
                             segments = 51, level = c(0.95, 0.99), plot.points = FALSE, ...) 
{
  dat <- my.dummy.df(data)
  names(dat) <- c("X1", "X2")
  n <- nrow(dat)
  p <- ncol(dat)
  if (!(is.vector(center) && 2 == length(center))) 
    stop("center must be a vector of length 2")
  if (scale == TRUE) {
    data <- cor(dat)
  } else {
    data <- cov(dat)
  }
  if (radius == "chi") {
    stat.mult <- sapply(level, function(x) sqrt(qchisq(x, 2)))
  } else {
    stat.mult <- sapply(level, function(x) sqrt((((n - 1) * p)/(n - p)) * qf(x, p, n - p)))
  }
  angles <- (0:segments) * 2 * pi/segments
  unit.circle <- cbind(cos(angles), sin(angles))
  DVt <- (diag(sqrt(svd(data)$d)) %*% t(svd(data)$v))[1:2, ]
  Ellipse <- llply(stat.mult, function(x) {
    Out <- data.frame(t(center + t(x * unit.circle %*% DVt)))
    names(Out) <- c("Axis 1", "Axis 2")
    Out[, 1:2]
  })
  "Axis 1" <- NULL
  "Axis 2" <- NULL
  "X1" <- NULL
  "X2" <- NULL
  output <- with(Ellipse[[1]], ggplot(Ellipse[[1]], aes(`Axis 1`, `Axis 2`)) + 
                   theme_bw() + 
                   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
                   geom_path() + 
                   theme(legend.position = "none") + 
                   ggtitle("Confidence Ellipse") + 
                   geom_hline(yintercept = center[2], col = "lightgrey") + geom_vline(xintercept = center[1], 
                                                                                      col = "lightgrey") + 
                   theme(plot.title = element_text(size = 20)) + 
                   theme(axis.title.x = element_text(size = 20)) + 
                   theme(strip.text.x = element_text(size = 20, 
                                                     colour = "black", face = "bold")) + 
                   theme(strip.text.y = element_text(size = 20, 
                                                     colour = "black", face = "bold")) + 
                   theme(axis.title.y = element_text(size = 20, angle = 90)) + 
                   theme(axis.text.x = element_text(size = 15, angle = 90, vjust = 0.5, face = "bold")) + 
                   theme(axis.text.y = element_text(size = 15, angle = 0, face = "bold")))
  if (length(level) > 1) {
    output <- output + sapply(2:length(Ellipse), function(x) geom_path(data = Ellipse[[x]], 
                                                                       aes(`Axis 1`, `Axis 2`), col = x))
  }
  if (plot.points) {
    output <- output + geom_point(data = dat, aes(X1, X2), col = "blue", pch = 1, alpha = 0.5)
  }
  print(output)
  return(Ellipse)
}