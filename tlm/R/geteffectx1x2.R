geteffectx1x2 <-
function(object, x1, x2, level, nboot)
 {
  # Point estimate:
  M1 <- getM(object = object, x = x1, untransform = T, level = level)["Estimate"]
  M2 <- getM(object = object, x = x2, untransform = T, level = level)["Estimate"]
  DeltaM <- M2 - M1
  DeltaM100 <- 100 * (M2 / M1 - 1) 
  # Bootstrap CI:
  dd <- object$model$model
  effectboot <- boot(data = dd, statistic = effectboot, R = nboot, stype = "i", object = object, x1 = x1, x2 = x2, level = level)
  CIDM <- boot.ci(effectboot, conf = level, type = "perc", index = 1)[[4]][4:5]
  CIDM100 <- boot.ci(effectboot, conf = level, type = "perc", index = 2)[[4]][4:5]
  DM <- c(DeltaM, CIDM, DeltaM100, CIDM100)
  aux <- paste(c("lower", "upper"), 100 * level,  "%", sep = "")
  names(DM) <- c("DM", aux, "DM100", aux)
  DM
 }
