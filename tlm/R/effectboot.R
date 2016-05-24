effectboot <-
function(data, i, object, x1, x2, level)
 {
  d <- data[i, ]
  objecti <- object
  objecti$model <- update(object$model, data = d)
  M1i <- getM(object = objecti, x = x1, untransform = T, level = level)["Estimate"]
  M2i <- getM(object = objecti, x = x2, untransform = T, level = level)["Estimate"]
  c(M2i - M1i, 100 * (M2i / M1i - 1))
 }
