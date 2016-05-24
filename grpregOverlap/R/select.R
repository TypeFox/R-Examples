## function: select, experimental
# -------------------------------------------------------------------------------
select.grpregOverlap <- function(obj, criterion=c("BIC","AIC",
                                                   "GCV","AICc","EBIC"), 
                                  df.method=c("default","active"), 
                                  smooth=FALSE, ...) {
  obj.new <- obj
  obj.new$beta <- obj$beta.latent
  obj.new$group <- obj$grp.vec
  class(obj.new) <- 'grpreg'
  res <- select(obj.new, criterion=criterion, df.method=df.method, 
                smooth=smooth, ...)
  i <- which(obj.new$lambda == res$lambda) # back trace the index of optimal lambda
  res$beta.latent <- res$beta
  res$beta <- obj$beta[, i]
  res
}
# ------------------------------------------------------------------------------
