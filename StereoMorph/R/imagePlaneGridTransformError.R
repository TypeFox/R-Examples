imagePlaneGridTransformError <- function(p, nx, ny, grid) mean(sqrt(rowSums((grid - imagePlaneGridTransform(p, nx, ny))^2)))
