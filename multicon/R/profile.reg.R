Profile.reg <-
function(x.set, y.set, center="group", std=FALSE, nomiss=.8) {
  if(sum(dim(x.set) != dim(y.set))>0) {stop("x.set and y.set must have the same dimensions.")}
  if(center=="group" & std==F) {
    cent.x <- data.frame(scale2(x.set, scale=F))
    cent.y <- data.frame(y.set)
  }
  if(center=="group" & std==T) {
    cent.x <- data.frame(scale2(x.set))
    cent.y <- data.frame(scale2(y.set))
  }
  if(center=="grand" & std==F) {
    cent.x <- data.frame(x.set - mean(colMeans(x.set, na.rm=T)))
    cent.y <- data.frame(y.set)
  }
  if(center=="grand" & std==T) {
    cent.x <- data.frame((x.set - mean(colMeans(x.set, na.rm=T))) / popsd(unlist(x.set)))
    cent.y <- data.frame((y.set - mean(colMeans(y.set, na.rm=T))) / popsd(unlist(y.set)))
  }
  if(center=="none" & std==F) {
    cent.x <- data.frame(x.set)
    cent.y <- data.frame(y.set)
  }
  if(center=="none" & std==T) {
    stop("Standardizing without centering not allowed. Choose different options.")
  }
  out <- data.frame(t(mapply(lin.coef, as.data.frame(t(cent.x)), as.data.frame(t(cent.y)), nomiss=nomiss)))
  rownames(out) <- seq(1:nrow(out))
  colnames(out) <- c("Intercepts", "Slopes")
  return(out)
}
