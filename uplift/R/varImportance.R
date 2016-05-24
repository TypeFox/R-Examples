### Define as generic

varImportance <- function(x, ...)  UseMethod("varImportance")

varImportance.default <- function(x, ...)
  stop("uplift: No method implemented for this class of object")


######################################################################
# Variable importance upliftRF
######################################################################

varImportance.upliftRF <- function(x, n.trees = x$ntree, plotit = TRUE, normalize = TRUE, ...) {
  
  if (!inherits(x, "upliftRF"))
    stop("uplift: x is not of class upliftRF")
  if (n.trees > x$ntree)
    stop("uplift: n.trees cannot be greater than the number of fitted trees")
  if (n.trees < 1)
    stop("uplift: n.trees must be greater than 0")
  
  lx <- n.trees 
  imp.temp <- vector("list", lx)
  # Add values of split criteria for each variable for each tree
  for (i in 1:lx) {   
    imp.temp[[i]] <- tapply(x$trees[[i]]$s_bs.s.value, x$trees[[i]]$s_bs.var, sum)  
  }
  # Average the values over all trees
  imp.val <- unlist(imp.temp)
  cnames <- names(imp.val)
  imp <- sort(tapply(imp.val, cnames, sum) / lx, decreasing = TRUE)
  # Relative importance of a variable can be <0 if eucli.lr < eucli.node (for "ED" criteria)
  imp <- ifelse(imp < 0, 0, imp) 
  
  if (normalize) imp <- 100 * imp / sum(imp)
  var_names <- x$var.names[as.numeric(names(imp))]
  res <- data.frame(var = var_names, rel.imp = imp, row.names = NULL)
  
  if (plotit) {
    barplot(res$rel.imp,
            horiz = TRUE,
            col= "red",
            names.arg = res$var,
            xlab = "Relative importance",...)
  }
  
  return(res)
  
}

### END FUN


######################################################################
# Variable importance ccif
######################################################################


varImportance.ccif <- function(x, n.trees = x$ntree, plotit = TRUE, normalize = TRUE, ...) {
  
  if (!inherits(x, "ccif"))
    stop("uplift: x is not of class ccif")
  if (n.trees > x$ntree)
    stop("uplift: n.trees cannot be greater than the number of fitted trees")
  if (n.trees < 1)
    stop("uplift: n.trees must be greater than 0")
  
  lx <- n.trees 
  imp.temp <- vector("list", lx)
  # Add values of split criteria for each variable for each tree
  for (i in 1:lx) {   
    imp.temp[[i]] <- tapply(x$trees[[i]]$s_bs.s.value, x$trees[[i]]$s_bs.var, sum)  
  }
  # Average the values over all trees
  imp.val <- unlist(imp.temp)
  cnames <- names(imp.val)
  imp <- sort(tapply(imp.val, cnames, sum) / lx, decreasing = TRUE)
  # Relative importance of a variable can be <0 if eucli.lr < eucli.node (for "ED" criteria)
  imp <- ifelse(imp < 0, 0, imp) 
  
  if (normalize) imp <- 100 * imp / sum(imp)
  var_names <- x$var.names[as.numeric(names(imp))]
  res <- data.frame(var = var_names, rel.imp = imp, row.names = NULL)
  
  if (plotit) {
    barplot(res$rel.imp,
            horiz = TRUE,
            col= "red",
            names.arg = res$var,
            xlab = "Relative importance",...)
  }
  
  return(res)
  
}

### END FUN