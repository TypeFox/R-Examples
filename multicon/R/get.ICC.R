get.ICC <-
function(x) {
  if(sum(is.na(x)) > 0) {warning("Missing observations in x. Remove missing data or interpret results with caution.")}
  GM <- mean(x, na.rm=T)
  ObsMeans <- rowMeans(x, na.rm=T)
  ObsNs <- apply(x, 1, function(x) length(x) - sum(is.na(x)))
  SSobs <- sum(ObsNs*(ObsMeans - GM)^2)
  ItemMeans <- colMeans(x, na.rm=T)
  ItemNs <- apply(x, 2, function(x) length(x) - sum(is.na(x)))
  SSitems <- sum(ItemNs*(ItemMeans - GM)^2)
  SSTot <- sum((x - GM)^2, na.rm=T)
  SSres <- SSTot - SSobs - SSitems
  DFobs <- nrow(x) - 1
  DFitems <- ncol(x) - 1
  DFtot <- nrow(x)*ncol(x) - 1
  DFres <- DFtot - DFobs - DFitems
  MSobs <- SSobs / DFobs
  MSitems <- SSitems / DFitems
  MSres <- SSres / DFres
  MSwith <- (SSitems + SSres) / (DFitems + DFres)
  ICC1 <- (MSobs - MSwith) / (MSobs + DFitems*MSwith)
  ICC1k <- (MSobs - MSwith) / MSobs
  ICC2 <- (MSobs - MSres)/(MSobs + DFitems * MSres + ncol(x) * (MSitems - MSres)/nrow(x))
  ICC2k <- (MSobs - MSres)/(MSobs + (MSitems - MSres)/nrow(x))
  ICC3 <- (MSobs - MSres) / (MSobs + DFitems*MSres)
  ICC3k <- (MSobs - MSres) / MSobs
  out <- cbind(ICC1, ICC1k, ICC2, ICC2k, ICC3, ICC3k)
  return(out)
}
